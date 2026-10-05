# Compute the tasmax-only annual climate stress metrics (see metrics.py) for
# CMIP5 MPI-ESM-LR, historical + RCP8.5.
#
# Source data: gs://clim_data_reg_useast1/cmip5/daily/tasmax/*.nc
# NOTE: despite the .nc extension, these are zip archives (CDS zips daily
# CMIP5 downloads for this dataset) each containing one CMOR-named NetCDF
# file, e.g. tasmax_day_MPI-ESM-LR_historical_r1i1p1_19700101-19791231.nc.
#
# Source files are downloaded, unzipped, read, and deleted one at a time
# (sequentially) -- each is a modest ~100-270MB, so there's no need to hold
# more than one in memory/on disk at once. Within each file, the per-year
# metric computation (the actual CPU-bound work) is parallelized across
# calendar years via mirai. The (small, ~130 annual layers per metric)
# results are then assembled into one NetCDF per metric and uploaded to:
#   gs://clim_data_reg_useast1/cmip5/annual_aggregates/{metric}/{metric}_{model}_CMIP5_rcp85_yr.nc

library(tidyverse)
library(stars)
library(units)
library(fs)
library(mirai)

# Provides gcs_list_objects() / gcs_upload(), both backed by the `gcloud` CLI.
source("gcs_io.R")

# ─── config ────────────────────────────────────────────────────────────────

VARIABLE <- "tasmax"
MODEL <- "MPI-ESM-LR"
SCRATCH_DIR <- path("/mnt/pd-pf")
DAILY_PREFIX <- str_glue("gs://clim_data_reg_useast1/cmip5/daily/{VARIABLE}/")
OUTPUT_PREFIX <- "gs://clim_data_reg_useast1/cmip5/annual_aggregates"
KELVIN_TO_CELSIUS <- 273.15
N_WORKERS <- 10 # years-within-a-file are processed independently; tune to available cores/RAM

# ─── metric registry (tasmax-only entries from metrics.py's METRIC_REGISTRY) ─
# Each function reduces one pixel's vector of daily tasmax values (deg C, for
# a single calendar year) to a single annual value.

tasmax_metrics <- list(
  "days-above-32c" = function(x) sum(x >= 32, na.rm = TRUE),
  "days-above-35c" = function(x) sum(x >= 35, na.rm = TRUE),
  "days-above-38c" = function(x) sum(x >= 38, na.rm = TRUE),
  "days-above-45c" = function(x) sum(x >= 45, na.rm = TRUE),
  "days-above-50c" = function(x) sum(x >= 50, na.rm = TRUE),
  "average-daytime-temperature" = function(x) mean(x, na.rm = TRUE),
  "ten-hottest-days" = function(x) {
    valid <- x[!is.na(x)]
    if (length(valid) == 0) {
      return(NA_real_)
    }
    k <- min(10, length(valid))
    mean(sort(valid, decreasing = TRUE)[seq_len(k)])
  },
  "freezing-days" = function(x) sum(x < 0, na.rm = TRUE)
)

# ─── per-file worker ─────────────────────────────────────────────────────────
#
# Runs sequentially in the main process (not shipped to a mirai daemon), so it
# relies on packages already attached at the top of this script. Downloads
# gs_uri's zip archive to scratch_dir, extracts the single NetCDF it contains,
# computes every metric's annual layer for every calendar year present in the
# file (parallelized across years via mirai_map(), using the daemons started
# by the caller), then deletes both the zip and the extracted NetCDF.
#
# Returns a named list keyed by calendar year (as character), each element a
# named list of stars 2D (lon, lat) layers keyed by metric name.
process_daily_file <- function(gs_uri, scratch_dir, metrics) {
  zip_dest <- path(scratch_dir, path_file(gs_uri))
  dl_status <- system2(
    "gcloud",
    c("storage", "cp", shQuote(gs_uri), shQuote(as.character(zip_dest))),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!is.null(attr(dl_status, "status")) && attr(dl_status, "status") != 0) {
    stop(str_glue(
      "gcloud storage cp failed for {gs_uri}:\n{str_flatten(dl_status, '\n')}"
    ))
  }

  inner_name <- utils::unzip(as.character(zip_dest), list = TRUE)$Name
  stopifnot(length(inner_name) == 1)
  extract_dir <- path(
    scratch_dir,
    str_glue("extract_{path_ext_remove(path_file(gs_uri))}")
  )
  dir_create(extract_dir)
  utils::unzip(
    as.character(zip_dest),
    files = inner_name,
    exdir = as.character(extract_dir)
  )
  nc_path <- path(extract_dir, inner_name)

  daily <- read_mdim(as.character(nc_path))
  # tasmax arrives in Kelvin; every threshold/metric above is defined in deg C.
  daily[[1]] <- units::drop_units(daily[[1]]) - 273.15

  time_vals <- st_get_dimension_values(daily, "time")
  years <- sort(unique(year(time_vals)))

  message(str_glue("processing {first(years)} - {last(years)}"))

  sliced_daily <- map(years, \(yr) dplyr::filter(daily, year(time) == yr))

  out <- mirai_map(
    sliced_daily,
    \(sl) {
      purrr::map(m, \(fn) stars::st_apply(sl, c("lon", "lat"), fn))
    },
    m = metrics
  )[]

  names(out) <- as.character(years)

  year_failures <- keep(out, ~ inherits(.x, "miraiError"))
  if (length(year_failures) > 0) {
    stop(str_glue(
      "{length(year_failures)} year(s) failed to process for {gs_uri}:\n",
      "{str_flatten(map_chr(year_failures, conditionMessage), '\n')}"
    ))
  }

  file_delete(zip_dest)
  dir_delete(extract_dir)

  out
}

# ─── list source files (historical before rcp85, chronological within each) ──
# Filenames embed zero-padded start dates and "historical" sorts before
# "rcp85" alphabetically, so a plain sort already yields chronological order.
daily_files <- sort(gcs_list_objects(DAILY_PREFIX))
stopifnot(length(daily_files) > 0)

# ─── process files sequentially; parallelize per-year metrics within each ────
# One file lives on disk/in memory at a time. Daemons are started once and
# reused across files, since each file's process_daily_file() call parallelizes
# its own per-year metric computation via mirai_map().
daemons(N_WORKERS)
per_file_results <- map(
  daily_files,
  process_daily_file,
  scratch_dir = SCRATCH_DIR,
  metrics = tasmax_metrics
)
daemons(0)

# ─── flatten to one entry per calendar year, in chronological order ──────────
# Source periods are calendar-year-aligned and non-overlapping, so every year
# should appear in exactly one file.
by_year <- list_flatten(per_file_results)
years <- as.integer(names(by_year))
stopifnot(!any(duplicated(years)))
by_year <- by_year[order(years)]
years <- sort(years)

# ─── assemble + write one annual NetCDF per metric ───────────────────────────
for (metric_name in names(tasmax_metrics)) {
  message(str_glue("Assembling {metric_name} ({length(years)} years)..."))

  layers <- map(by_year, pluck, metric_name)
  annual_dates <- as.Date(str_glue("{years}-01-01"))
  combined <- do.call(c, c(layers, list(along = list(time = annual_dates))))
  st_crs(combined) <- 4326

  local_out <- path(
    SCRATCH_DIR,
    str_glue("{metric_name}_{MODEL}_CMIP5_rcp85_yr.nc")
  )
  write_mdim(combined, as.character(local_out)) |>
    suppressWarnings()

  gs_dest <- str_glue(
    "{OUTPUT_PREFIX}/{metric_name}/{metric_name}_{MODEL}_CMIP5_rcp85_yr.nc"
  )
  message(str_glue("Uploading {gs_dest} ..."))
  gcs_upload(as.character(local_out), gs_dest)
  file_delete(local_out)
}

message("Done.")

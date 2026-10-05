# Compute warming-level aggregates for each annual climate stress metric (see
# scripts/tasmax_metrics_cmip5.R), for CMIP5 MPI-ESM-LR, historical + RCP8.5.
#
# For each metric and each warming level of interest (0.5 to 3 deg C, every
# 0.5 deg C), find the 21-year window centered on the year at which the model
# crosses that warming level (per the GCM warming-level lookup table below),
# and compute the per-cell min, 5th/50th/95th percentile, max, and mean across
# that window. The 0.5 deg C level isn't in the lookup table, so it instead
# uses the model's baseline 1971-2000 (30-year) period.
#
# Source: gs://clim_data_reg_useast1/cmip5/annual_aggregates/{metric}/{metric}_{model}_CMIP5_rcp85_yr.nc
# Output: gs://clim_data_reg_useast1/cmip5/warming_level_aggregates/{metric}/{metric}_{model}_CMIP5_rcp85_wls.nc
#
# Annual aggregate files are small (global, coarse ~1.875deg grid, ~130
# timesteps), so unlike the daily-source scripts in this project, each is read
# directly from GCS via GDAL's /vsigs/ virtual filesystem (no local
# download/scratch file needed).

library(tidyverse)
library(stars)
library(fs)

# Provides gcs_upload(), backed by the `gcloud` CLI.
source("gcs_io.R")

# ─── config ────────────────────────────────────────────────────────────────

MODEL <- "MPI-ESM-LR"
SCRATCH_DIR <- path("/mnt/pd-pf")
INPUT_PREFIX <- "gs://clim_data_reg_useast1/cmip5/annual_aggregates"
OUTPUT_PREFIX <- "gs://clim_data_reg_useast1/cmip5/warming_level_aggregates"
WARMING_LEVELS <- seq(0.5, 3, by = 0.5)
WINDOW_HALF_WIDTH <- 10 # 21-year window = central year +/- 10
BASELINE_YEARS <- 1971:2000 # used for the 0.5 deg C level (no central year in the table)

# tasmax-only annual metrics already produced by tasmax_metrics_cmip5.R
METRICS <- c(
  "days-above-32c",
  "days-above-35c",
  "days-above-38c",
  "days-above-45c",
  "days-above-50c",
  "average-daytime-temperature",
  "ten-hottest-days",
  "freezing-days"
)

gs_to_vsigs <- function(gs_uri) str_replace(gs_uri, "^gs://", "/vsigs/")

# ─── warming-level central years for this model ───────────────────────────
# Probable Futures' table of the central year at which each CMIP5 GCM's
# running-mean global temperature first crosses each warming level (relative
# to pre-industrial), under RCP8.5.
wl_table <- read_tsv(
  "https://raw.githubusercontent.com/Probable-Futures/map-data-processing/refs/heads/main/cmip5_model_temp_thresholds.csv",
  show_col_types = FALSE
) |>
  filter(Model == MODEL)
stopifnot(nrow(wl_table) == 1)

# Named integer vector, e.g. c(`1` = 2006, `1.5` = 2025, `2` = 2038, ...);
# the table has no column for 0.5, handled separately via BASELINE_YEARS.
central_years <- WARMING_LEVELS |>
  discard(~ .x == 0.5) |>
  set_names() |>
  map_int(\(wl) wl_table[[str_glue("Yr{wl}")]])

message(
  "Central years for ",
  MODEL,
  ": ",
  str_flatten_comma(str_glue("{names(central_years)}={central_years}"))
)

# ─── per-cell summary stats across a window of annual layers ──────────────

window_stats <- function(v) {
  c(
    min = min(v, na.rm = TRUE),
    p05 = unname(quantile(v, 0.05, na.rm = TRUE)),
    p50 = unname(quantile(v, 0.50, na.rm = TRUE)),
    p95 = unname(quantile(v, 0.95, na.rm = TRUE)),
    max = max(v, na.rm = TRUE),
    mean = mean(v, na.rm = TRUE)
  )
}

# Slice `annual` (a lon/lat/time stars object) to `years`, then reduce the
# time dimension to per-cell summary stats. Returns a lon/lat stars object
# with one attribute per stat (min, p05, p50, p95, max, mean).
wl_summary <- function(annual, years) {
  annual |>
    filter(year(time) %in% years) |>
    st_apply(c("lon", "lat"), window_stats) |>
    split("window_stats")
}

years_for_wl <- function(wl) {
  if (wl == 0.5) {
    BASELINE_YEARS
  } else {
    central_year <- central_years[[as.character(wl)]]
    (central_year - WINDOW_HALF_WIDTH):(central_year + WINDOW_HALF_WIDTH)
  }
}

# ─── process each metric ────────────────────────────────────────────────────

for (metric_name in METRICS) {
  gs_dest <- str_glue(
    "{OUTPUT_PREFIX}/{metric_name}/{metric_name}_{MODEL}_CMIP5_rcp85_wls.nc"
  )

  # Idempotency: skip metrics that have already been processed.
  if (length(gcs_list_objects(gs_dest)) > 0) {
    message(str_glue("Skipping {metric_name}: {gs_dest} already exists."))
    next
  }

  message(str_glue("Processing {metric_name}..."))

  annual_uri <- str_glue(
    "{INPUT_PREFIX}/{metric_name}/{metric_name}_{MODEL}_CMIP5_rcp85_yr.nc"
  )
  annual <- read_mdim(gs_to_vsigs(annual_uri))

  available_years <- year(st_get_dimension_values(annual, "time"))

  wl_layers <- map(WARMING_LEVELS, \(wl) {
    years <- years_for_wl(wl)
    missing <- setdiff(years, available_years)
    if (length(missing) > 0) {
      stop(str_glue(
        "{metric_name}: warming level {wl} window is missing years ",
        "{str_flatten_comma(missing)} from the source data."
      ))
    }
    wl_summary(annual, years)
  })

  combined <- do.call(c, c(wl_layers, list(along = list(wl = WARMING_LEVELS))))
  st_crs(combined) <- st_crs(annual)

  local_out <- path(
    SCRATCH_DIR,
    str_glue("{metric_name}_{MODEL}_CMIP5_rcp85_wls.nc")
  )
  write_mdim(combined, as.character(local_out)) |>
    suppressWarnings()

  message(str_glue("Uploading {gs_dest} ..."))
  gcs_upload(as.character(local_out), gs_dest)
  file_delete(local_out)
}

message("Done.")

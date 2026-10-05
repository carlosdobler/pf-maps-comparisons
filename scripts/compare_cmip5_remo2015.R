# Compare CMIP5 (MPI-ESM-LR, global ~1.875deg) warming-level aggregates
# against REMO2015-downscaled (MPI-M-MPI-ESM-LR, ~0.2deg regional) warming-
# level aggregates, for the 8 tasmax-only metrics already computed for CMIP5
# (see tasmax_metrics_cmip5.R / warming_level_aggregates_cmip5.R).
#
# REMO2015 warming-level aggregates don't exist as a product -- they're built
# here from REMO2015's annual per-metric files (one per domain) using the
# *same* breaching-year table/windows used for CMIP5 (same underlying GCM,
# MPI-ESM-LR, so the same warming-level-crossing years apply). These
# REMO2015 aggregates are never written to disk/bucket; they only exist
# in-memory for the comparison below.
#
# For each metric x domain:
#   1. Load REMO2015's annual file for that domain (gs_to_vsigs, no local
#      download), fix its (sometimes CMOR-broken) time coordinate, convert
#      K -> degC where needed, and reduce to per-cell min/p05/p50/p95/max/mean
#      over the same 21-year (or 1971-2000 baseline, for 0.5 deg C) windows
#      used for CMIP5.
#   2. Load CMIP5's existing warming-level aggregate file, restrict it to the
#      warming levels of interest, convert its 0-360 deg longitude to
#      -180/180 (REMO2015 domains straddle the prime meridian or antimeridian
#      in some cases) and regularize the lon/lat dimensions (both sources
#      have CF dimensions that are technically irregular due to float32
#      round-trip noise, which GDAL's warp machinery rejects), then
#      nearest-neighbor-regrid (st_warp) it onto the REMO2015 domain's grid.
#   3. Difference = regridded CMIP5 - REMO2015, per stat and warming level.
#      Uploaded to:
#      gs://clim_data_reg_useast1/cmip5/comparisons_remo2015/warming_level_aggregates/{metric}/{metric}_{domain}_MPI-ESM-LR_CMIP5-vs-REMO2015_rcp85_wls.nc
#
# Both sources are small (global ~1.875deg or one regional ~0.2deg domain,
# a handful of warming levels), so everything is read directly from GCS via
# GDAL's /vsigs/ virtual filesystem; only the final small diff NetCDF touches
# local disk (briefly, before upload).

library(tidyverse)
library(stars)
library(units)
library(abind)
library(fs)
library(mirai)

# Provides gcs_list_objects() / gcs_upload(), both backed by the `gcloud` CLI.
source("gcs_io.R")

# ─── config ────────────────────────────────────────────────────────────────

MODEL <- "MPI-ESM-LR"
SCRATCH_DIR <- path("/mnt/pd-pf")
CMIP5_WL_PREFIX <- "gs://clim_data_reg_useast1/cmip5/warming_level_aggregates"
REMO_PREFIX <- "gs://clim_data_reg_useast1/results/global_heat_pf/01_derived"
OUTPUT_PREFIX <- "gs://clim_data_reg_useast1/cmip5/comparisons_remo2015/warming_level_aggregates"
WARMING_LEVELS <- c(0.5, 1, 2, 3)
WINDOW_HALF_WIDTH <- 10 # 21-year window = central year +/- 10
BASELINE_YEARS <- 1971:2000 # used for the 0.5 deg C level (no central year in the table)
DOMAINS <- c(
  "AFR",
  "AUS",
  "CAM",
  "CAS",
  "EAS",
  "EUR",
  "NAM",
  "SAM",
  "SEA",
  "WAS"
)

# REMO2015 annual files cover 1970-2100 (131 contiguous calendar years), one
# layer per year, but some (the "ten-hottest-*" family) have a broken CF time
# coordinate ("day as %Y%m%d.%f", a known CMOR/CDO quirk -- see the
# "fix-time" flag in Probable Futures' tb_vars_all variable table). Rather
# than special-case decoding, every REMO2015 annual file's time dimension is
# just overwritten with this fixed 1970:2100 sequence (values are already in
# chronological order within the file).
REMO_YEARS <- 1970:2100

# Map our (CMIP5) metric names to REMO2015's "var_derived" names; see
# https://raw.githubusercontent.com/Probable-Futures/map-data-processing/refs/heads/main/scripts_v3/tb_vars_all.R
# (the table there is missing the 50C entry; confirmed via bucket listing
# that "days-gte-50C-tasmax" exists with a "_v3-01" filename suffix, hence
# the wildcard match in find_remo_uri() below).
METRIC_MAP <- tribble(
  ~metric_name                  , ~var_derived              ,
  "days-above-32c"              , "days-gte-32C-tasmax"     ,
  "days-above-35c"              , "days-gte-35C-tasmax"     ,
  "days-above-38c"              , "days-gte-38C-tasmax"     ,
  "days-above-45c"              , "days-gte-45C-tasmax"     ,
  "days-above-50c"              , "days-gte-50C-tasmax"     ,
  "average-daytime-temperature" , "mean-tasmax"             ,
  "ten-hottest-days"            , "ten-hottest-days-tasmax" ,
  "freezing-days"               , "days-lt-0C-tasmax"
)

# REMO2015 vars whose values are temperatures in Kelvin (need -> degC to
# match CMIP5, which is already in degC); the rest are day counts.
KELVIN_VARS <- c("mean-tasmax", "ten-hottest-days-tasmax")

gs_to_vsigs <- function(gs_uri) str_replace(gs_uri, "^gs://", "/vsigs/")

# ─── warming-level central years for this model (shared by both sources) ───
# Probable Futures' table of the central year at which each CMIP5 GCM's
# running-mean global temperature first crosses each warming level (relative
# to pre-industrial), under RCP8.5. REMO2015 is a dynamical downscaling of
# this same GCM, so the same crossing years apply to it.
wl_table <- read_tsv(
  "https://raw.githubusercontent.com/Probable-Futures/map-data-processing/refs/heads/main/cmip5_model_temp_thresholds.csv",
  show_col_types = FALSE
) |>
  filter(Model == MODEL)
stopifnot(nrow(wl_table) == 1)

central_years <- WARMING_LEVELS |>
  discard(~ .x == 0.5) |>
  set_names() |>
  map_int(\(wl) wl_table[[str_glue("Yr{wl}")]])

years_for_wl <- function(wl) {
  if (wl == 0.5) {
    BASELINE_YEARS
  } else {
    central_year <- central_years[[as.character(wl)]]
    (central_year - WINDOW_HALF_WIDTH):(central_year + WINDOW_HALF_WIDTH)
  }
}

# ─── geometry helpers ────────────────────────────────────────────────────
#
# Both CMIP5's wl-aggregate NetCDFs and REMO2015's annual NetCDFs encode
# lon/lat as explicit per-cell interval bounds rather than a clean
# offset+delta sequence (float32 round-trip noise from the source data).
# GDAL's warp machinery (used by st_warp) requires a true affine (regular)
# grid, so both are "regularized" below before any warp/write.

regularize_dim <- function(x, dim_name) {
  d <- st_dimensions(x)[[dim_name]]
  if (!is.null(d$offset) && !is.na(d$offset)) {
    return(x) # already regular
  }
  starts <- d$values$start
  n <- length(starts)
  delta <- (starts[n] - starts[1]) / (n - 1)
  crs <- st_crs(x)
  x <- st_set_dimensions(x, dim_name, offset = starts[1], delta = delta)
  st_crs(x) <- crs
  x
}

# Drops a length-1 dimension (e.g. a single wl or time slice) from a stars
# object, keeping all attributes. st_warp() rejects a 2D (x,y) raster that
# still carries a spare length-1 third dimension.
drop_singleton <- function(x, dim_name) {
  dims <- st_dimensions(x)
  dpos <- which(names(dims) == dim_name)
  stopifnot(dim(x)[dpos] == 1)
  arrs <- map(x, ~ abind::adrop(.x, drop = dpos))
  crs <- st_crs(x)
  out <- st_as_stars(
    set_names(arrs, names(x)),
    dimensions = dims[setdiff(names(dims), dim_name)]
  )
  st_crs(out) <- crs
  out
}

# CMIP5's longitude runs 0-360; REMO2015 domains use -180/180 and several
# (e.g. AFR) straddle the prime meridian, so a plain 0-360 grid would warp
# incorrectly. Shifts to -180/180, re-sorting columns to keep the grid
# regular (lon spacing divides 360 evenly, so this is just a circular shift).
shift_to_pm180 <- function(x) {
  lon <- st_get_dimension_values(x, "lon")
  new_lon <- ((lon + 180) %% 360) - 180
  ord <- order(new_lon)
  delta <- st_dimensions(x)[["lon"]]$delta
  crs <- st_crs(x)
  x2 <- x[, ord, , ]
  x2 <- st_set_dimensions(
    x2,
    "lon",
    offset = new_lon[ord][1] - delta / 2,
    delta = delta
  )
  st_crs(x2) <- crs
  x2
}

# ─── per-cell summary stats across a window of annual layers ──────────────
# NA-safe: an all-NA window (REMO2015 cells outside the domain / coastline
# fringe) returns NA for every stat rather than +-Inf (what base min()/max()
# return on an empty/all-NA vector, with a warning).
window_stats <- function(v) {
  if (all(is.na(v))) {
    return(c(
      min = NA_real_,
      p05 = NA_real_,
      p50 = NA_real_,
      p95 = NA_real_,
      max = NA_real_,
      mean = NA_real_
    ))
  }
  c(
    min = min(v, na.rm = TRUE),
    p05 = unname(quantile(v, 0.05, na.rm = TRUE)),
    p50 = unname(quantile(v, 0.50, na.rm = TRUE)),
    p95 = unname(quantile(v, 0.95, na.rm = TRUE)),
    max = max(v, na.rm = TRUE),
    mean = mean(v, na.rm = TRUE)
  )
}

wl_summary <- function(s) {
  s |>
    stars::st_apply(c("lon", "lat"), window_stats) |>
    split("window_stats")
}

# ─── REMO2015 side: build warming-level aggregates for one domain/metric ──

find_remo_uri <- function(domain, var_derived) {
  pattern <- str_glue(
    "{REMO_PREFIX}/{domain}_{var_derived}_yr_REMO2015_MPI-M-MPI-ESM-LR*.nc"
  )
  matches <- gcs_list_objects(pattern)
  if (length(matches) != 1) {
    stop(str_glue(
      "Expected exactly 1 REMO2015 file for domain={domain}, var={var_derived}; ",
      "found {length(matches)}."
    ))
  }
  matches
}

build_remo_wl <- function(domain, var_derived) {
  uri <- find_remo_uri(domain, var_derived)
  annual <- read_mdim(gs_to_vsigs(uri)) |> suppressWarnings()

  n_time <- dim(annual)[["time"]]
  stopifnot(n_time == length(REMO_YEARS))
  annual <- st_set_dimensions(
    annual,
    "time",
    values = as.Date(str_glue("{REMO_YEARS}-01-01"))
  )

  if (var_derived %in% KELVIN_VARS) {
    annual[[1]] <- units::drop_units(annual[[1]]) - 273.15
  } else {
    annual[[1]] <- units::drop_units(annual[[1]])
  }

  annual_split <- map(WARMING_LEVELS, \(wl) {
    annual |>
      filter(year(time) %in% years_for_wl(wl))
  })

  daemons(5)
  wl_layers <-
    # NA-only windows near domain edges/coastlines
    suppressWarnings(
      mirai_map(
        annual_split,
        wl_summary,
        window_stats = window_stats
      )[]
    )
  daemons(0)

  combined <- do.call(c, c(wl_layers, list(along = list(wl = WARMING_LEVELS))))
  st_crs(combined) <- 4326

  combined |> regularize_dim("lon") |> regularize_dim("lat")
}

# ─── CMIP5 side: load + reproject-ready an existing wl-aggregate file ─────

load_cmip5_wl <- function(metric_name) {
  uri <- str_glue(
    "{CMIP5_WL_PREFIX}/{metric_name}/{metric_name}_{MODEL}_CMIP5_rcp85_wls.nc"
  )
  x <- read_mdim(gs_to_vsigs(uri))
  x <- x |> filter(wl %in% WARMING_LEVELS)
  st_crs(x) <- 4326
  x |> shift_to_pm180() |> regularize_dim("lon") |> regularize_dim("lat")
}

# Nearest-neighbor regrids every wl slice of `cmip5_wl` onto `template`
# (a 2D lon/lat stars object) and recombines along wl.
# regrid_cmip5 <- function(cmip5_wl, template) {
#   warp_one <- function(wl_idx) {
#     drop_singleton(cmip5_wl[,,, wl_idx], "wl") |>
#       st_warp(dest = template, method = "near")
#   }
#   out <- map(seq_along(WARMING_LEVELS), warp_one)
#   combined <- do.call(c, c(out, list(along = list(wl = WARMING_LEVELS))))
#   st_crs(combined) <- 4326
#   combined
# }

# ─── process every metric x domain ────────────────────────────────────────

for (i in seq_len(nrow(METRIC_MAP))) {
  metric_name <- METRIC_MAP$metric_name[i]
  var_derived <- METRIC_MAP$var_derived[i]

  message(str_glue("=== {metric_name} ({var_derived}) ==="))
  cmip5_wl <- load_cmip5_wl(metric_name)

  for (domain in DOMAINS) {
    gs_dest <- str_glue(
      "{OUTPUT_PREFIX}/{metric_name}/{metric_name}_{domain}_{MODEL}_CMIP5-vs-REMO2015_rcp85_wls.nc"
    )

    if (length(gcs_list_objects(gs_dest)) > 0) {
      message(str_glue("  skipping {domain}: {gs_dest} already exists."))
      next
    }

    message(str_glue("  {domain}..."))

    remo_wl <- build_remo_wl(domain, var_derived)
    template_2d <- remo_wl |> slice(wl, 1) #drop_singleton(remo_wl[names(remo_wl)[1], , , 1], "wl")
    cmip5_regridded <- st_warp(cmip5_wl, template_2d)

    stopifnot(identical(dim(cmip5_regridded), dim(remo_wl)))
    stopifnot(identical(names(cmip5_regridded), names(remo_wl)))

    diff_obj <- cmip5_regridded - remo_wl

    local_out <- path(
      SCRATCH_DIR,
      str_glue("{metric_name}_{domain}_{MODEL}_CMIP5-vs-REMO2015_rcp85_wls.nc")
    )
    write_mdim(diff_obj, as.character(local_out)) |> suppressWarnings()

    gcs_upload(as.character(local_out), gs_dest)
    file_delete(local_out)
  }
}

message("Done.")

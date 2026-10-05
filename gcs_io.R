# Helpers for reading/writing Google Cloud Storage objects used by both dev
# and production. Two access patterns are used in this project:
#
#   1. Direct, disk-free reads of NetCDF via GDAL's /vsigs/ virtual
#      filesystem, using stars::read_mdim() (used for ERA5, and optionally
#      for formatted NMME).
#   2. Upload/listing of locally-written formatted/output NetCDFs via the
#      `gcloud storage` CLI (mirrors the old pipeline's approach; confirmed
#      available and authenticated in this environment).

library(stringr)
library(fs)
library(stars)

#' List objects under a gs:// prefix.
#'
#' @param gs_prefix character, e.g. "gs://water-balance/era5/monthly/total_precipitation/"
#' @return character vector of full gs:// URIs (character(0) if the prefix
#'   doesn't exist or is empty)
gcs_list_objects <- function(gs_prefix) {
  result <- system2(
    "gcloud",
    c("storage", "ls", shQuote(gs_prefix)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    # gcloud storage ls exits non-zero when the prefix matches nothing.
    return(character(0))
  }
  result
}

# #' Read a single NetCDF object directly from GCS without writing to local disk.
# #'
# #' Uses the /vsigs/ GDAL virtual filesystem via stars::read_mdim(), so
# #' multi-dimensional NetCDFs (e.g. ERA5's monthly files) can be read (and,
# #' with proxy = TRUE, lazily subset) without ever touching local disk.
# #'
# #' @param gs_uri character, full gs:// URI to a NetCDF file
# #' @param ... passed through to stars::read_mdim()
# #' @param proxy logical, passed to stars::read_mdim(); TRUE (default) defers
# #'   actually reading data until the result is subset/materialized (e.g. via
# #'   st_as_stars())
# #' @return a stars (or stars_proxy, if proxy = TRUE) object
# read_nc_from_gcs <- function(gs_uri, ..., proxy = TRUE) {
#   stars::read_mdim(gs_to_vsigs(gs_uri), ..., proxy = proxy)
# }

# #' Download a single GCS object to local disk.
# #'
# #' Symmetric counterpart to gcs_upload(); used for small objects (e.g.
# #' frozen EMOS parameter .rds files) where the /vsigs/ NetCDF-oriented
# #' read_nc_from_gcs() path doesn't apply.
# #'
# #' @param gs_uri character, source gs:// URI (a single object, not a prefix)
# #' @param local_path character, destination local file path
# #' @return character, `local_path`, invisibly
# gcs_download <- function(gs_uri, local_path) {
#   status <- system2(
#     "gcloud",
#     c("storage", "cp", shQuote(gs_uri), shQuote(local_path)),
#     stdout = TRUE,
#     stderr = TRUE
#   )
#   exit_status <- attr(status, "status")
#   if (!is.null(exit_status) && exit_status != 0) {
#     stop(str_glue(
#       "gcloud storage cp failed ({gs_uri} -> {local_path}):\n",
#       "{str_flatten(status, '\n')}"
#     ))
#   }
#   invisible(local_path)
# }

#' Upload a local file to a GCS destination.
#'
#' @param local_path character, path to a local file
#' @param gs_uri character, destination gs:// URI (file or prefix)
gcs_upload <- function(local_path, gs_uri) {
  status <- system2(
    "gcloud",
    c("storage", "cp", shQuote(local_path), shQuote(gs_uri)),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && exit_status != 0) {
    stop(str_glue(
      "gcloud storage cp failed ({local_path} -> {gs_uri}):\n",
      "{str_flatten(status, '\n')}"
    ))
  }
  invisible(gs_uri)
}

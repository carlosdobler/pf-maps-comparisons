"""
Annual climate stress metrics.

Each public function accepts a dict of xr.DataArray (keyed by variable name)
and returns an annual xr.DataArray with a 'time' dimension of yearly Jan-1
timestamps.

The METRIC_REGISTRY at the bottom maps each kebab-case metric name to its
function and the list of input variables it requires.
"""

from __future__ import annotations

import dask
import numpy as np
import pandas as pd
import xarray as xr

from config import LAT_COORD, WB_PERCENTILE_FIRST_VALID_YEAR


# ─── private helpers ──────────────────────────────────────────────────────────


def _annual_sum(condition: xr.DataArray) -> xr.DataArray:
    """Count True days per calendar year."""
    return condition.resample(time="YS").sum()


def _top_n_mean(da: xr.DataArray, n: int) -> xr.DataArray:
    """Annual mean of the n highest daily values.

    Each calendar year's top-n mean depends only on that year's values, so
    there's no cross-year dependency to preserve. Rather than rechunking the
    entire (multi-decade) time dimension to a single chunk — which would
    force each dask task to hold the full time series for its spatial tile
    in memory — we slice the data to one calendar year at a time and only
    rechunk that (much smaller) slice to -1 before calling apply_ufunc. Each
    per-year task then needs roughly
    n_days_per_year * lat_chunk * lon_chunk * 4 bytes, about two orders of
    magnitude less than the full-series approach, at the cost of a larger
    dask graph (one task per spatial tile per year instead of one per
    spatial tile). See CHUNKS in config.py.
    """
    unique_years = np.unique(da.time.dt.year.values)

    def _fn(x: np.ndarray) -> float:
        valid = x[~np.isnan(x)]
        if len(valid) == 0:
            return np.nan
        k = min(n, len(valid))
        return np.mean(np.partition(valid, -k)[-k:])

    annual_results = []
    for yr in unique_years:
        da_year = da.sel(time=str(yr)).chunk({"time": -1})
        result_year = xr.apply_ufunc(
            _fn,
            da_year,
            input_core_dims=[["time"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[float],
        )
        annual_results.append(
            result_year.expand_dims(time=[pd.Timestamp(f"{yr}-01-01")])
        )

    return xr.concat(annual_results, dim="time")


# ─── tasmax metrics ───────────────────────────────────────────────────────────


def days_above_32c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmax"] >= 32.0)


def days_above_35c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmax"] >= 35.0)


def days_above_38c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmax"] >= 38.0)


def days_above_45c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmax"] >= 45.0)


def days_above_50c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmax"] >= 50.0)


def average_daytime_temperature(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return ds["tasmax"].resample(time="YS").mean()


def ten_hottest_days(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _top_n_mean(ds["tasmax"], 10)


def freezing_days(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmax"] < 0.0)


# ─── tasmin metrics ───────────────────────────────────────────────────────────


def frost_nights(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmin"] < 0.0)


def nights_above_20c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmin"] >= 20.0)


def nights_above_25c(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["tasmin"] >= 25.0)


def average_nighttime_temperature(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return ds["tasmin"].resample(time="YS").mean()


def ten_hottest_nights(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _top_n_mean(ds["tasmin"], 10)


# ─── tas metrics ─────────────────────────────────────────────────────────────


def average_temperature(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return ds["tas"].resample(time="YS").mean()


def average_winter_temperature(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    """Hemisphere-aware annual winter mean temperature.

    Northern hemisphere (lat >= 0): mean over DJF months of that calendar year
    (i.e. Dec of year Y is averaged with Jan and Feb of year Y, not year Y+1).
    Southern hemisphere (lat < 0): mean over JJA months of that calendar year.
    """
    tas = ds["tas"]
    lat_mask = tas[LAT_COORD] >= 0  # True = northern hemisphere

    north = (
        tas.where(tas.time.dt.month.isin([12, 1, 2]))
        .where(lat_mask)
        .resample(time="YS")
        .mean()
    )
    south = (
        tas.where(tas.time.dt.month.isin([6, 7, 8]))
        .where(~lat_mask)
        .resample(time="YS")
        .mean()
    )
    return xr.where(lat_mask, north, south)


# ─── wet-bulb metrics ─────────────────────────────────────────────────────────


def days_above_26c_wbmax(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["wetbulbmax"] >= 26.0)


def days_above_28c_wbmax(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["wetbulbmax"] >= 28.0)


def days_above_30c_wbmax(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["wetbulbmax"] >= 30.0)


def days_above_32c_wbmax(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _annual_sum(ds["wetbulbmax"] >= 32.0)


def ten_hottest_wbmax_days(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _top_n_mean(ds["wetbulbmax"], 10)


# ─── precipitation metrics ────────────────────────────────────────────────────


def total_annual_precipitation(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return ds["pr"].resample(time="YS").sum()


def wettest_day(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return ds["pr"].resample(time="YS").max()


def wettest_90_days(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    """Maximum non-overlapping 90-day precipitation sum per year.

    The 90-day rolling sum can span calendar year boundaries (e.g. a window
    ending in January may include days from the previous December), and each
    year's selected window must not overlap the previous year's. Rather than
    rechunking the entire multi-decade time dimension into a single chunk
    (which forces each dask task to hold the full time series for its
    spatial tile in memory), this is computed as a lazy chain of small
    per-year tasks: for year Y we slice in only that year plus its trailing
    89-day lookback buffer (~454 days total, versus the full series), and
    thread the previous year's selected-window end date through the chain as
    an ordinary (small, lat/lon-shaped) array argument. The non-overlap
    comparison happens inside the per-pixel numpy function at runtime, so no
    eager `.compute()` is needed to carry that state forward — the whole
    per-year loop stays lazy, at the cost of a longer sequential dask-task
    chain per spatial tile (years can no longer be computed concurrently for
    a given tile, though different tiles remain fully parallel). Each task
    now needs roughly
    454 * lat_chunk * lon_chunk * 4 bytes, about two orders of magnitude
    less than the full-series approach. See CHUNKS in config.py.

    Note: this also fixes a latent bug in the previous full-series
    implementation, where `np.argmax` over a candidate range containing NaN
    rolling-sum values (from not yet having 90 days of history, e.g. in the
    first calendar year of the record) could spuriously return a NaN
    position -- producing a NaN result for that year and, since the bogus
    position was carried forward as the "previous window" state, corrupting
    the following year's non-overlap constraint too. Candidates with a NaN
    rolling sum are now explicitly excluded.

    Checkpointing: although each per-year task is small, `state` still forms
    a single dask-graph dependency chain spanning the entire record (one
    task per year per spatial tile, each depending on the last). Left fully
    lazy, that unbounded chain -- combined with a single `xr.concat` + write
    at the end -- lets completed-but-not-yet-consumed intermediate results
    accumulate across all years and tiles before anything is flushed to
    durable storage, which can exhaust local disk via dask's spill-to-disk
    mechanism on a long run. To bound this, every `CHECKPOINT_EVERY_N_YEARS`
    years we eagerly materialize (`dask.compute()`) *both* `state` and that
    batch's annual sums together in a single call, then re-wrap each as a
    small constant (re-chunked to match the input's spatial chunks for
    `state`). Computing them jointly is essential: `state` and the annual
    sum for a given year share the same upstream `apply_ufunc` task, so
    computing only `state` and discarding the (still-lazy) sum would force
    that entire upstream chain -- including re-reading the source data from
    GCS -- to be redone later when the sum is finally needed, silently
    doubling the work for every checkpointed year. Materializing both
    together lets dask execute the shared graph exactly once.
    """
    precip = ds["pr"]
    unique_years = np.unique(precip.time.dt.year.values)

    CHECKPOINT_EVERY_N_YEARS = 10

    spatial_chunks = {
        dim: precip.chunks[precip.get_axis_num(dim)][0]
        for dim in precip.dims
        if dim != "time"
    }

    # Sentinel meaning "no previous window" -- guarantees the first year's
    # candidates are unconstrained (day_offsets are real calendar day
    # numbers since a fixed epoch, always >> this value).
    SENTINEL = -1_000_000

    state = xr.zeros_like(precip.isel(time=0, drop=True), dtype=np.int64) + SENTINEL

    annual_sums = []
    pending = []  # indices into annual_sums computed since the last checkpoint
    for i, yr in enumerate(unique_years):
        buffer_start = pd.Timestamp(f"{yr - 1}-12-31") - pd.Timedelta(days=89)
        year_end = pd.Timestamp(f"{yr}-12-31")
        da_year = precip.sel(time=slice(buffer_start, year_end)).chunk({"time": -1})

        local_dates = da_year["time"].values
        day_offsets = local_dates.astype("datetime64[D]").astype(np.int64)
        year_mask = pd.DatetimeIndex(local_dates).year.values == yr
        n_local = len(local_dates)

        def _fn(
            x: np.ndarray,
            prev_state: np.int64,
            day_offsets: np.ndarray = day_offsets,
            year_mask: np.ndarray = year_mask,
            n_local: int = n_local,
        ) -> tuple[float, np.int64]:
            if np.all(np.isnan(x)):
                return np.nan, prev_state

            # 90-day backward rolling sum via cumsum (NaN treated as 0)
            filled = np.where(np.isnan(x), 0.0, x)
            cs = np.concatenate([[0.0], np.cumsum(filled)])
            runsum = np.full(n_local, np.nan)
            if n_local >= 90:
                runsum[89:] = cs[90:] - cs[: n_local - 89]

            candidates = (
                year_mask & ~np.isnan(runsum) & (day_offsets >= prev_state + 90)
            )
            if not np.any(candidates):
                return np.nan, prev_state

            idx = np.where(candidates)[0]
            best = idx[np.argmax(runsum[idx])]
            return runsum[best], day_offsets[best]

        result_year, state = xr.apply_ufunc(
            _fn,
            da_year,
            state,
            input_core_dims=[["time"], []],
            output_core_dims=[[], []],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[float, np.int64],
        )
        annual_sums.append(result_year.expand_dims(time=[pd.Timestamp(f"{yr}-01-01")]))
        pending.append(len(annual_sums) - 1)

        if (i + 1) % CHECKPOINT_EVERY_N_YEARS == 0:
            batch = xr.concat([annual_sums[j] for j in pending], dim="time")
            batch, state_values = dask.compute(batch, state)
            for offset, j in enumerate(pending):
                annual_sums[j] = batch.isel(time=[offset])
            state = xr.DataArray(
                state_values, dims=state.dims, coords=state.coords
            ).chunk(spatial_chunks)
            pending = []

    return xr.concat(annual_sums, dim="time")


# ─── multi-variable metrics ───────────────────────────────────────────────────


def snowy_days(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    """Count days with precipitation ≥ 1 mm and mean temperature < 0°C."""
    return _annual_sum((ds["pr"] >= 1.0) & (ds["tas"] < 0.0))


# ─── water balance percentile metrics ─────────────────────────────────────────


def _drop_incomplete_first_year(da: xr.DataArray) -> xr.DataArray:
    """Drop calendar year 1961, which lacks a full w12 lookback window.

    wb_percentile uses a trailing 12-month rolling window, so Jan-Nov 1961
    are NaN (no prior-year data to look back on) and Dec 1961 alone isn't a
    meaningful annual aggregate. See WB_PERCENTILE_FIRST_VALID_YEAR in
    config.py.
    """
    return da.sel(time=slice(str(WB_PERCENTILE_FIRST_VALID_YEAR), None))


def average_water_balance(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _drop_incomplete_first_year(ds["wb_percentile"].resample(time="YS").mean())


def probability_of_drought(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _drop_incomplete_first_year(
        (ds["wb_percentile"] <= 30.0).resample(time="YS").mean()
    )


def probability_of_extreme_drought(ds: dict[str, xr.DataArray]) -> xr.DataArray:
    return _drop_incomplete_first_year(
        (ds["wb_percentile"] <= 5.0).resample(time="YS").mean()
    )


# ─── registry ─────────────────────────────────────────────────────────────────

METRIC_REGISTRY: dict[str, dict] = {
    "days-above-32c": {"fn": days_above_32c, "variables": ["tasmax"]},
    "days-above-35c": {"fn": days_above_35c, "variables": ["tasmax"]},
    "days-above-38c": {"fn": days_above_38c, "variables": ["tasmax"]},
    "days-above-45c": {"fn": days_above_45c, "variables": ["tasmax"]},
    "days-above-50c": {"fn": days_above_50c, "variables": ["tasmax"]},
    "average-temperature": {"fn": average_temperature, "variables": ["tas"]},
    "average-daytime-temperature": {
        "fn": average_daytime_temperature,
        "variables": ["tasmax"],
    },
    "ten-hottest-days": {"fn": ten_hottest_days, "variables": ["tasmax"]},
    "freezing-days": {"fn": freezing_days, "variables": ["tasmax"]},
    "frost-nights": {"fn": frost_nights, "variables": ["tasmin"]},
    "nights-above-20c": {"fn": nights_above_20c, "variables": ["tasmin"]},
    "nights-above-25c": {"fn": nights_above_25c, "variables": ["tasmin"]},
    "average-nighttime-temperature": {
        "fn": average_nighttime_temperature,
        "variables": ["tasmin"],
    },
    "average-winter-temperature": {
        "fn": average_winter_temperature,
        "variables": ["tas"],
    },
    "ten-hottest-nights": {"fn": ten_hottest_nights, "variables": ["tasmin"]},
    "days-above-26c-wbmax": {"fn": days_above_26c_wbmax, "variables": ["wetbulbmax"]},
    "days-above-28c-wbmax": {"fn": days_above_28c_wbmax, "variables": ["wetbulbmax"]},
    "days-above-30c-wbmax": {"fn": days_above_30c_wbmax, "variables": ["wetbulbmax"]},
    "days-above-32c-wbmax": {"fn": days_above_32c_wbmax, "variables": ["wetbulbmax"]},
    "ten-hottest-wbmax-days": {
        "fn": ten_hottest_wbmax_days,
        "variables": ["wetbulbmax"],
    },
    "total-annual-precipitation": {
        "fn": total_annual_precipitation,
        "variables": ["pr"],
    },
    "wettest-90-days": {"fn": wettest_90_days, "variables": ["pr"]},
    "snowy-days": {"fn": snowy_days, "variables": ["pr", "tas"]},
    "wettest-day": {"fn": wettest_day, "variables": ["pr"]},
    "average-water-balance": {
        "fn": average_water_balance,
        "variables": ["wb_percentile"],
    },
    "probability-of-drought": {
        "fn": probability_of_drought,
        "variables": ["wb_percentile"],
    },
    "probability-of-extreme-drought": {
        "fn": probability_of_extreme_drought,
        "variables": ["wb_percentile"],
    },
}

# /// script
# dependencies = [
#     "cdsapi",
#     "google-cloud-storage",
# ]
# ///
"""Download CMIP5 daily variables (MPI-ESM-LR, historical + RCP8.5) from the
Climate Data Store (CDS) and stream them to Google Cloud Storage.

Files are downloaded one at a time to a scratch directory on the attached
persistent disk (never the boot disk), immediately uploaded to
`gs://clim_data_reg_useast1/cmip5/daily/{variable}/`, and then removed from
local disk so nothing accumulates there. Historical periods are processed
before RCP8.5 periods.

Usage:
    uv run scripts/download_cmip5.py tasmax
"""

import argparse
import sys
from pathlib import Path

import cdsapi
from google.cloud import storage

DATASET = "projections-cmip5-daily-single-levels"
MODEL = "mpi_esm_lr"
MODEL_NAME = "MPI-ESM-LR"
ENSEMBLE_MEMBER = "r1i1p1"

SCRATCH_DIR = Path("/mnt/pd-pf")
GCS_BUCKET = "clim_data_reg_useast1"
GCS_PREFIX = "cmip5/daily/{variable}"

# Maps our short variable name to the CDS variable name.
VARIABLES = {
    "tasmax": "maximum_2m_temperature_in_the_last_24_hours",
}

# Maps our filename experiment label to the CDS "experiment" request value,
# along with the periods to request for that experiment.
EXPERIMENTS = {
    "historical": {
        "cds_experiment": "historical",
        "periods": [
            "19700101-19791231",
            "19800101-19891231",
            "19900101-19991231",
            "20000101-20051231",
        ],
    },
    "rcp85": {
        "cds_experiment": "rcp_8_5",
        "periods": [
            "20060101-20091231",
            "20100101-20191231",
            "20200101-20291231",
            "20300101-20391231",
            "20400101-20491231",
            "20500101-20591231",
            "20600101-20691231",
            "20700101-20791231",
            "20800101-20891231",
            "20900101-21001231",
        ],
    },
}


def build_filename(variable: str, experiment: str, period: str) -> str:
    return f"{variable}_{MODEL_NAME}_CMIP6_{experiment}_day_{period}.nc"


def download_one(client: cdsapi.Client, variable: str, experiment_key: str, period: str, dest: Path) -> None:
    cds_experiment = EXPERIMENTS[experiment_key]["cds_experiment"]
    request = {
        "experiment": cds_experiment,
        "variable": [VARIABLES[variable]],
        "model": MODEL,
        "ensemble_member": ENSEMBLE_MEMBER,
        "period": [period],
    }
    client.retrieve(DATASET, request).download(str(dest))


def upload_and_remove(bucket: storage.Bucket, local_path: Path, variable: str) -> None:
    blob_name = f"{GCS_PREFIX.format(variable=variable)}/{local_path.name}"
    blob = bucket.blob(blob_name)
    blob.upload_from_filename(str(local_path))
    local_path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download CMIP5 daily data (MPI-ESM-LR) and upload it to GCS."
    )
    parser.add_argument(
        "variable",
        choices=sorted(VARIABLES),
        help="Short name of the variable to download (e.g. 'tasmax').",
    )
    args = parser.parse_args()
    variable = args.variable

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

    client = cdsapi.Client()
    storage_client = storage.Client()
    bucket = storage_client.bucket(GCS_BUCKET)

    # Historical before RCP8.5, one file at a time end-to-end.
    for experiment_key in ("historical", "rcp85"):
        periods = EXPERIMENTS[experiment_key]["periods"]
        for period in periods:
            filename = build_filename(variable, experiment_key, period)
            dest = SCRATCH_DIR / filename
            print(f"Downloading {filename} ...", file=sys.stderr)
            download_one(client, variable, experiment_key, period, dest)
            print(f"Uploading {filename} to gs://{GCS_BUCKET}/{GCS_PREFIX.format(variable=variable)}/ ...", file=sys.stderr)
            upload_and_remove(bucket, dest, variable)
            print(f"Done with {filename}.", file=sys.stderr)


if __name__ == "__main__":
    main()

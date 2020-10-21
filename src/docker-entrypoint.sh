#!/bin/bash
set -e

python3 -u /data-ingest/src/pre_run_test.py
Rscript /data-ingest/src/1-preproc.r
python3 -u /data-ingest/src/2-run_scrublet.py
Rscript /data-ingest/src/3-preproc.r
python3 -u /data-ingest/src/4-create-adata.py
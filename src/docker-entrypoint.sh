#!/bin/bash
set -e

./pre_run_test.py
Rscript /data-ingest/src/1-preproc.r
/data-ingest/src/2-run_scrublet.py
Rscript /data-ingest/src/3-preproc.r
/data-ingest/src/4-create-adata.py
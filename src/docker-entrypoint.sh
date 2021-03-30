#!/bin/bash
set -e

#python3 -u /data-ingest/src/pre_run_test.py
#Rscript /data-ingest/src/1_Preproc.r ${INPUT_TYPE-10x}
#Rscript /data-ingest/src/2-1_Compute-metris_emptyDrops.r
#python3 -u /data-ingest/src/2-2_Compute-metrics_scrublet.py
#Rscript /data-ingest/src/3_Seurat.r
Rscript /data-ingest/src/4_Prepare_experiment.r
python3 -u /data-ingest/src/5_upload-to-aws.py
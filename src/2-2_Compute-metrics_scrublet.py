#!/usr/bin/python3

################################################
## 2-2_Compute-metrics_scrublet.r
##  - Compute emptyDrops per sample
##  - Save an rds with emptyDrops result per sample
################################################

import pandas
import scrublet as scr
import json
import os
import copy

with open("/input/meta.json") as f:
    config = json.load(f)

# Check which samples have been selected. Otherwiser we are going to use all of them. 
samples = config["samples"]
if len(samples) == 0:
    samples = [
        name
        for name in os.listdir("/input")
        if os.path.isdir(os.path.join("/input", name))
    ]

print("Now running scrublet...")

for sample in samples:
    # Reading pre-filtered matrix. From 1_Preproc.r we have stored one pre-doublet-matrix per sample
    df = pandas.read_csv("/output/pre-doublet-matrix-" + sample + ".csv", index_col=0)
    print("Running sample " + sample)
    # Computing scrublet
    scrub = scr.Scrublet(df)
    # Checking when filtering inside `scrub_doublets` reduces the matrix less than 30 genes or 30 cells
    # 30 is the hardcode value for PCA analysis in scrublet
    scrub_check = copy.deepcopy(scrub)
    # This pipeline belongs to the pipeline inside `scrub_doublets` (https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py)
    scr.pipeline_normalize(scrub_check)
    scr.pipeline_get_gene_filter(scrub_check, min_counts=3, min_cells=3, min_gene_variability_pctl=85)
    scr.pipeline_apply_gene_filter(scrub_check)
    # The dimension must be strictly less than min(n_samples, n_features), so we substract 1
    min_dimension = min(scrub_check._E_obs.shape)-1 
    n_prin_comps = min(30, min_dimension)
    print("Number of principal components :" + str(n_prin_comps))
    doublet_scores_sample, _ = scrub.scrub_doublets(n_prin_comps=n_prin_comps)
    # We will store the doublet-scores in the following format
    ####  Barcode   score
    with open("/output/doublet-scores-" + sample + ".csv", "w") as f:
        f.writelines(
            (
                df.index[i] + "\t" + f"{doublet_scores_sample[i]}\n"
                for i in range(len(df.index))
            )
        )

print("Step 2-2 completed.")

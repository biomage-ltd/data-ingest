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
    doublet_scores_sample, _ = scrub.scrub_doublets()
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

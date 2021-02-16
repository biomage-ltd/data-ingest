#!/usr/bin/python3
import pandas
import numpy as np
import scrublet as scr
import json

print("Now running scrublet...")
df = pandas.read_csv('/output/pre-doublet-matrix.csv', index_col=0)

with open('/input/meta.json') as f:
    config = json.load(f)

# Checking multisample
if config['samples']['multisample'] == 'TRUE':
  samples = config['samples']['samples_info']['type']
  doublet_scores = []
  # Running the rest
  for i in range(0, len(samples)):
    id_sample = samples[i]
    print("Running sample " + id_sample)
    sample = df.loc[df.index.str.contains(id_sample)]
    scrub = scr.Scrublet(sample)
    doublet_scores_sample, _ = scrub.scrub_doublets()
    doublet_scores = np.append(doublet_scores, doublet_scores_sample)
else:
  scrub = scr.Scrublet(df)
  doublet_scores, _ = scrub.scrub_doublets()

print("Exporting doublet scores...")
with open('/output/doublet-scores.csv', 'w') as f:
    f.writelines((f"{score}\n" for score in doublet_scores))

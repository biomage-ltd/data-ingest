#!/usr/bin/python3
import pandas
import scrublet as scr

print("Now running scrublet...")
df = pandas.read_csv('/output/pre-doublet-matrix.csv')

scrub = scr.Scrublet(df)
doublet_scores, _ = scrub.scrub_doublets()

print("Exporting doublet scores...")
with open('/output/doublet-scores.csv', 'w') as f:
    f.writelines((f"{score}\n" for score in doublet_scores))

#!/bin/bash

# run data-ingest for each folder in user_data
for d in ./user_data/*; do
 f="$(basename -- $d)"
 echo "Working on $f ..."
 rm -rf input
 cp -R $d/. input
 docker-compose up --build

 # add to exp_meta.csv file
 Rscript ./src/wrapper/add_exp_meta.R

 # can't remove output for next as from docker root user
 # can rename
 rename output output_$f output
done


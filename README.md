data-ingest
===========

An application for ingesting data into Biomage

Set-up
------

Make sure you have Docker, docker-compose installed.
Make sure your input folder contains the required files: `barcodes.tsv`, `genes.tsv`, `matrix.mtx`. Make sure your output folder is empty. These conditions will be checked by the script before it runs.

Then run:

    EXPERIMENT_NAME="My Experiment" AWS_ACCESS_KEY_ID=key AWS_SECRET_ACCESS_KEY=key  docker-compose up --build

to start the ingest.
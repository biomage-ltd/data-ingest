data-ingest
===========

An application for ingesting data into Biomage

How to run the script to load user's data into the platform.
------

1. Make sure you have Docker, docker-compose installed.

2. Make sure you have an input and output folders at the root of this project.

3. Go to S3 to `biomage-originals-production` bucket and navigate to the relevant folder with analysis files you want to load in the platform. The folder should contain files `barcodes.tsv`, `genes.tsv`, `matrix.mtx`.

4. Download the `barcodes.tsv`, `genes.tsv`, `matrix.mtx` files, located in S3 from the previous step. Make sure that they are saved inside the `input` folder you created in step 2.

5. Make sure your output folder is empty. This folder will be used by the script to populate all relevant output files from the preprocessing steps.

6. Run:

    EXPERIMENT_NAME="My Experiment" AWS_ACCESS_KEY_ID=key AWS_SECRET_ACCESS_KEY=key  docker-compose up --build

to start the ingest.
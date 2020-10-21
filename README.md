data-ingest
===========

An application for ingesting data into Biomage. 

This repo contains a dockerised script that does data processing and filtering on experiment result files and loads the results into the Biomage single cell platform. This consists for the following steps:

How to run the script to load user's data into the platform.
------

1. Make sure you have Docker, docker-compose installed.

2. Make sure you have an input and output folders at the root of this project. Make sure that these folders are empty. Otherwise their contents can pollute the results after running the script.

3. Go to S3 to `biomage-originals-production` bucket. You should see several folders. Each of them contains unprocessed experiment results for a user. Navigate to the correct folder with analysis files you want to preprocess and load in the platform. The folder should contain files `barcodes.tsv`, `genes.tsv`, `matrix.mtx`.

4. Download the `barcodes.tsv`, `genes.tsv`, `matrix.mtx` files, located in S3 from the previous step. Make sure that they are saved inside the `input` folder you created in step 2.

5. Run:

    EXPERIMENT_NAME="My Experiment" AWS_ACCESS_KEY_ID=key AWS_SECRET_ACCESS_KEY=key  docker-compose up --build

to start the ingest.
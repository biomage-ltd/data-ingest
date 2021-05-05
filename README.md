data-ingest
===========

An application for ingesting data into Biomage

Re-ingest everything
------

Run from `data-ingest` folder:

```bash
# sync from s3 to local
aws s3 sync s3://biomage-originals-production ./user_data/

CLUSTER_ENV=all ./src/wrapper/ingest-all.sh
```

This will:

1. sync all data in `biomage-originals-producting` to `./user_data/`
2. run data-ingest using each folder in `./user_data/*` as input
3. create/append experiment meta-data to a `exp_meta.csv` file

Set-up
------

This repo contains a dockerised script that does data processing and filtering on
experiment result files and loads the results into the Biomage single cell platform.
This consists for the following steps:

1. Make sure you have Docker, docker-compose installed. Make sure your AWS credentials
are correct in `~/.aws`.

2. You might need to increase the amount of RAM being used by
the virtual machine. Go to the Docker icon in your macOS
bar, and click Preferences. Click Resources, and increase the memory slider.

3. Make sure you have an `input/` and `output/` folders at the root of this project.
Make sure that these folders are empty. Otherwise their contents can pollute the
results after running the script. The only exception is the meta.json file that
**must** be present in the `input/` folder.

4. Go to S3 to `biomage-originals-production` bucket. You should see several folders.
Each of them contains unprocessed experiment results for a user. Navigate to the correct
folder with analysis files you want to preprocess and load in the platform. The folder
should contain files `barcodes.tsv`, `genes.tsv`, `matrix.mtx`.

5. Download the `barcodes.tsv`, `genes.tsv`, `matrix.mtx` files, located in S3 from the
previous step. Make sure that they are saved inside the `input` folder you created in step 2.

6. Fill out the `meta.json` file. This is the configuration that is being used to save
and process the files appropriately.

`name` should be the name of the experiment. If you don't know this, you can use the name
of the folder in the S3 bucket.

`organism` must match the organism of the data set. The appropriate organism ID can be found
[here](https://biit.cs.ut.ee/gprofiler/page/organism-list). For example, Human is `hsapiens`,
Chicken is `ggallus`,  Zebrafish is `drerio`, etc.

`input` should not be modified for 10x data sets.

5. Run: `CLUSTER_ENV=production docker-compose up --build` to upload to production. The default will
upload to staging. You can also specify `CLUSTER_ENV=all` to upload to both staging and production.

This process should have given you an Experiment-ID (EID) as an output and also uploaded 3 things:

1. The experiment configuration in `json` format into **DynamoDB** table `experiments-{CLUSTER_ENV}` with the EID as `experimentId` partion key
2. The sample configuration in `json` format into **DynamoDB** table `samples-{CLUSTER_ENV}` with the EID as `experimentId` partion key
3. An `r.rds` object into **S3** this bucket:  `biomage-source-{CLUSTER_ENV}/{experiment_id}/r.rds`

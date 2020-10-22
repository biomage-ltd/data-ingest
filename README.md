data-ingest
===========

An application for ingesting data into Biomage

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

5. Run:

    docker-compose up --build

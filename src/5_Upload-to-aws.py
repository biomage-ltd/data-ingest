#!/usr/bin/python3

################################################
## 5_upload-to-aws.py
##  - Getting ready samples-table and experiments-table
##  - Upload tables to DynamoDB
##  - Upload experiment.rds to S3
################################################

import hashlib
import os
import pandas
from scipy.io import mmread
import matplotlib.pyplot as plt
import boto3
import json
from decimal import Decimal
from datetime import datetime
import uuid

COLOR_POOL = []

# Environment to deploy to, by default staging
ENVIRONMENT = os.getenv("ENVIRONMENT", "staging")

WARN_TXT_COL = "\033[93m"
RESET_TXT_COL = '\033[0m'
print(f"{WARN_TXT_COL}Deploying to {ENVIRONMENT}{RESET_TXT_COL}")

with open("/data-ingest/src/color_pool.json") as f:
    COLOR_POOL = json.load(f)


def calculate_checksum(filenames):
    hash = hashlib.md5()
    for fn in filenames:
        if os.path.isfile(fn):
            hash.update(open(fn, "rb").read())
    return hash.hexdigest()

# This function crate the table information for samples. As input it requires the experiment id and the config.  
def create_samples_table(config, experiment_id):
    # In samples_table we are going to add the core of the information
    samples_table = {}

    # Getting flag_filtered information
    df_prefilered = pandas.read_csv(
        "/output/df_flag_filtered.txt",
        sep="\t",
        na_values=["None"],
    )

    # Firstly, we identify the samples name. To do that we fetch the names of the folders (we suppose that the name
    # of the folders corresponds with the samples name) or direclty get them from the config
    if len(config["samples"])>1:
        samples = config["samples"]
    else:
        samples = [ name for name in os.listdir("/input") if os.path.isdir(os.path.join("/input", name)) ]
    
    samples_table["ids"] = [ "sample-"+sample for sample in samples]

    # For the current datasets it could happen that they are not in the gz format, so we leave the alternative tsv format. 
    mime_options = {"tsv": "application/tsv", "gz" : "application/gzip", "mtx" :  "application/mtx"}

    for sample in samples:
        
        # flag filtered
        preFiltered = df_prefilered.loc[df_prefilered.samples==sample, "flag_filtered"].tolist()[0]=="Filtered"

        # Identify datetime
        createdDate = datetime.now()
        lastModified = datetime.now()
        fileNames = {}
        # Look for the file that are not hidden (the hidden files start with .hidden.tsv)
        sample_files = [sample+"/"+f for f in os.listdir("/input/"+sample) if not f.startswith('.')]

        # Iterate over each file to create the slot
        for sample_file in sample_files:
            fileNames[sample_file] = {
                "objectKey" : '', 
                "name" : sample_file, 
                "size" : os.stat("/input/"+sample_file).st_size, 
                "mime": mime_options[sample_file.split(".")[-1]],
                "success": True, 
                "error": False
            }

        # Add the whole information to each sample
        samples_table["sample-"+sample] = {
            "name" : sample, 
            "uuid" : str(uuid.uuid4()), 
            "species": config["organism"],
            "type": config["input"]["type"],
            "createdDate": createdDate.isoformat(),
            "lastModified": lastModified.isoformat(),
            "complete": True, 
            "error": False,
            "fileNames" : sample_files, 
            "files" : fileNames,
            "preFiltered" : preFiltered
        }

    return {"experimentId": experiment_id, "samples" : samples_table}

# cell_sets fn for seurat samples name
def samples_sets():
    # construct new cell set group

    samples_annotations = pandas.read_csv(
        "/output/samples-cells.csv",
        sep="\t",
        names=["Cells_ID", "Value"],
        na_values=["None"],
    )

    cell_set = {
        "key": "sample",
        "name": "Samples",
        "rootNode": True,
        "children": [],
        "type": "metadataCategorical",
    }

    for sample in samples_annotations["Value"].unique():
        view = samples_annotations[samples_annotations.Value == sample]["Cells_ID"]
        cell_set["children"].append(
            {
                "key": f"sample-{sample}",
                "name": f"{sample}",
                "color": COLOR_POOL.pop(0),
                "cellIds": [int(d) for d in view.tolist()],
            }
        )

    return cell_set


# cell_sets fn for seurat metadata information
def meta_sets():
    meta_annotations = pandas.read_csv(
        "/output/metadata-cells.csv",
        sep="\t",
        header=0,
    )
    
    cell_set_list = list()

    # The first column is the cells_id, the rest is the metadata information
    for i in range(1, len(meta_annotations.columns)):

        # keep as key and name the name of the column
        key = meta_annotations.columns[i]
        name = meta_annotations.columns[i]

        cell_set = {
            "key": key,
            "name": name,
            "rootNode": True,
            "children": [],
            "type": "metadataCategorical",
        }

        for value in meta_annotations.iloc[:, i].unique():
            view = meta_annotations[meta_annotations.iloc[:, i] == value]["cells_id"]
            cell_set["children"].append(
                {
                    "key":   key+f"-{value}",
                    "name": f"{value}",
                    "color": COLOR_POOL.pop(0),
                    "cellIds": [int(d) for d in view.tolist()],
                }
            )

        cell_set_list.append(cell_set)
    return cell_set_list

def main():
    experiment_id = calculate_checksum(
        [
            "/output/r-out-raw.mtx",
            "/output/r-out-normalized.mtx",
            "/output/r-out-cells.tsv",
        ]
    )

    config = None
    with open("/input/meta.json", "r") as f:
        config = json.load(f)

    # read config related with QC pipeline
    config_dataProcessing = None
    with open("/output/config_dataProcessing.json", "r") as f:
        config_dataProcessing = json.load(f)


    # Design cell_set scratchpad for DynamoDB
    scratchpad = {
                "key": "scratchpad",
                "name": "Scratchpad",
                "rootNode": True,
                "children": [],
                "type": "cellSets",
            }

    samples_data = create_samples_table(config, experiment_id)
    samples_set = samples_sets()

    if "metadata" in config.keys():
        # Design cell_set meta_data for DynamoDB
        meta_set = meta_sets()
        cellSets = [samples_set, scratchpad] + meta_set

    else:
        # Design cell_set meta_data for DynamoDB
        cellSets = [scratchpad, samples_set]

    print("Experiment name is", config["name"])

    FILE_NAME = f"biomage-source-{ENVIRONMENT}/{experiment_id}/r.rds"

    experiment_data = {
        "apiVersion": "2.0.0-data-ingest-seurat-rds-automated",
        "experimentId": experiment_id,
        "experimentName": config["name"],
        "meta": {
            "organism": config["organism"],
            "type": config["input"]["type"],
        },
        "processingConfig": config_dataProcessing, 
    }

    cell_sets_data = json.dumps(cellSets)
    
    # Conver to float all decimals
    experiment_data = json.loads(json.dumps(experiment_data), parse_float=Decimal)
    samples_data = json.loads(json.dumps(samples_data), parse_float=Decimal)
    
    access_key = os.getenv("AWS_ACCESS_KEY_ID")
    secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")

    experiments_key = f"experiments-{ENVIRONMENT}"
    samples_key = f"samples-{ENVIRONMENT}"
    cell_sets_bucket = f"cell-sets-{ENVIRONMENT}"
    r_object_bucket, r_object_key = FILE_NAME.split("/", 1)

    print(f"uploading to dynamodb experiments table {experiments_key}...")
    dynamo = boto3.resource(
        "dynamodb",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    ).Table(experiments_key)
    dynamo.put_item(Item=experiment_data)

    print(f"uploading to dynamodb samples table {samples_key}...")
    dynamo = boto3.resource(
        "dynamodb",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    ).Table(samples_key)
    dynamo.put_item(Item=samples_data)

    print(f"uploading cellsets table to s3 bucket {cell_sets_bucket}...")
    s3 = boto3.client(
        "s3",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    )

    s3.put_object(Body=cell_sets_data, Bucket=cell_sets_bucket, Key=experiment_id)

    print(f"uploading R object to s3 bucket {r_object_bucket}...")
    s3 = boto3.client(
        "s3",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    )
    
    with open("/output/experiment.rds", "rb") as f:
        s3.put_object(Body=f, Bucket=r_object_bucket, Key=r_object_key)
    
    print("successful. experiment is now accessible at:")

    url_to_env = ENVIRONMENT == "https://scp.biomage.net" if "production" else "staging_deployment_url"
    print(f"{url_to_env}/experiments/{experiment_id}/data-exploration")



main()

print("Step 5 completed.")


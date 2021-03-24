#!/usr/bin/python3
import hashlib
import os
import pandas
from scipy.io import mmread
import matplotlib.pyplot as plt
import boto3
import json
from decimal import Decimal
from datetime import datetime

COLOR_POOL = []

with open("/data-ingest/src/color_pool.json") as f:
    COLOR_POOL = json.load(f)


def calculate_checksum(filenames):
    hash = hashlib.md5()
    for fn in filenames:
        if os.path.isfile(fn):
            hash.update(open(fn, "rb").read())
    return hash.hexdigest()


# This function crate the table information for samples. As input it requires the experiment id and the user uuid.  
def create_samples_table_multisample(config, experiment_id, uuid):
    # In samples_table we are going to add the core of the information
    samples_table = {}
    # Firstly, we identify the samples name. To do that we fetch the names of the folders (we suppose that the name
    # of the folders corresponds with the samples name)
    samples = [ name for name in os.listdir("/input") if os.path.isdir(os.path.join("/input", name)) ]
    samples_table["ids"] = samples
    # We crate a cnt variable to track the number of samples, since the key value of the is: sample1, sample2, ... samplen 
    # (with n as the total number of samples)
    cnt = 1
    # For the current datasets it could happen that they are not in the gz format, so we leave the alternative tsv format. 
    mime_options = {"tsv": "application/tsv", "gz" : "application/gzip", "mtx" :  "application/mtx"}

    for sample in samples:
        
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
        samples_table["sample"+str(cnt)] = {
            "name" : sample, 
            "uuid" : uuid, 
            "species": config["organism"],
            "type": config["input"]["type"],
            "createdDate": createdDate.isoformat(),
            "lastModified": lastModified.isoformat(),
            "complete": True, 
            "error": False,
            "fileNames" : sample_files, 
            "files" : fileNames
        }
        # Incremen the iterator for the sampleX
        cnt = cnt+1

    return {"experiment_id": experiment_id, "samples" : samples_table}


# This function crate the table information for samples. As input it requires the experiment id and the user uuid.  
def create_samples_table_unisample(config, experiment_id, uuid):
    # In samples_table we are going to add the core of the information
    samples_table = {}
    samples_table["ids"] = "sample1"

    # For the current datasets it could happen that they are not in the gz format, so we leave the alternative tsv format. 
    mime_options = {"tsv": "application/tsv", "gz" : "application/gzip", "mtx" :  "application/mtx"}

    createdDate = datetime.now()
    lastModified = datetime.now()
    fileNames = {}
    
    # We identify the files name. To do that we fetch the names of the files excep the meta.json one
    sample_files = [ f for f in os.listdir("/input") if f!="meta.json" and not f.startswith('.') ]

    for sample_file in sample_files:
        fileNames[sample_file] = {
            "objectKey" : '', 
            "name" : sample_file, 
            "size" : os.stat("/input/"+sample_file).st_size, 
            "mime": mime_options[sample_file.split(".")[-1]],
            "success": True, 
            "error": False
        }
    samples_table["sample1"] = {
        "name" : "sample1", 
        "uuid" : uuid, 
        "species": config["organism"],
        "type": config["input"]["type"],
        "createdDate": createdDate.isoformat(),
        "lastModified": lastModified.isoformat(),
        "complete": True, 
        "error": False,
        "fileNames" : sample_files, 
        "files" : fileNames
    }

    return {"experiment_id": experiment_id, "samples" : samples_table}


# cell_sets fn for seurat clusters
def cell_sets_seurat():
    # construct new cell set group
    cell_set = {
        "key": "louvain",
        "name": "Louvain clusters",
        "rootNode": True,
        "children": [],
        "type": "cellSets",
    }

    cluster_annotations = pandas.read_csv(
        "/output/cluster-cells.csv",
        sep="\t",
        names=["Cells_ID", "Clusters"],
        na_values=["None"],
    )
    
    for cluster in sorted(cluster_annotations["Clusters"].unique()):
        view = cluster_annotations[cluster_annotations.Clusters == cluster]["Cells_ID"]
        cell_set["children"].append(
            {
                "key": f"louvain-{cluster}",
                "name": f"Cluster {cluster}",
                "color": COLOR_POOL.pop(0),
                "cellIds": [int(d) for d in view.tolist()],
            }
        )

    return cell_set

# cell_sets fn for seurat samples name
def meta_sets():
    # construct new cell set group
    cell_set = {
        "key": "sample",
        "name": "Samples",
        "rootNode": True,
        "children": [],
        "type": "metadataCategorical",
    }

    meta_annotations = pandas.read_csv(
        "/output/multisample-cells.csv",
        sep="\t",
        names=["Cells_ID", "type"],
        na_values=["None"],
    )
    
    for sample in meta_annotations["type"].unique():
        view = meta_annotations[meta_annotations.type == sample]["Cells_ID"]
        cell_set["children"].append(
            {
                "key": f"sample-{sample}",
                "name": f"{sample}",
                "color": COLOR_POOL.pop(0),
                "cellIds": [int(d) for d in view.tolist()],
            }
        )

    return cell_set


def main():
    experiment_id = calculate_checksum(
        [
            "/output/r-out-raw.mtx",
            "/output/r-out-normalized.mtx",
            "/output/r-out-cells.tsv",
            "/output/r-out-dispersions.tsv",
        ]
    )

    config = None
    with open("/input/meta.json", "r") as f:
        config = json.load(f)

    # read config related with QC pipeline
    config_dataProcessing = None
    with open("/output/config_dataProcessing.json", "r") as f:
        config_dataProcessing = json.load(f)

    # Design cell_set cluster for DynamoDB
    cell_set = cell_sets_seurat()

    # Design cell_set scratchpad for DynamoDB
    scratchpad = {
                "key": "scratchpad",
                "name": "Scratchpad",
                "rootNode": True,
                "children": [],
                "type": "cellSets",
            }

    if config['samples']['multisample'] == 'TRUE':
        # Design cell_set meta_data for DynamoDB
        meta_set = meta_sets()
        cellSets = [cell_set, meta_set, scratchpad]

        # Design samples-table for DynamoDB
        samples_data = create_samples_table_multisample(config, experiment_id, "a1234")

    else:
        # Design cell_set meta_data for DynamoDB
        cellSets = [cell_set, scratchpad]
        # Design samples-table for DynamoDB
        samples_data = create_samples_table_unisample(config, experiment_id, "a1234")


    print("Experiment name is", config["name"])

    experiment_id = "multisample_unfiltered"

    FILE_NAME = f"biomage-source-production/{experiment_id}/r.rds"

    experiment_data = {
        "apiVersion": "2.0.0-data-ingest-seurat-rds-automated",
        "experimentId": experiment_id,
        "experimentName": config["name"],
        "meta": {
            "organism": config["organism"],
            "type": config["input"]["type"],
        },
        "matrixPath": FILE_NAME,
        "cellSets": cellSets,
        "processingConfig": config_dataProcessing, 
    }

    # Conver to float all decimals
    experiment_data = json.loads(json.dumps(experiment_data), parse_float=Decimal)
    samples_data = json.loads(json.dumps(samples_data), parse_float=Decimal)
    
    access_key = os.getenv("AWS_ACCESS_KEY_ID")
    secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")

    print("uploading to dynamodb experiments table...")
    dynamo = boto3.resource(
        "dynamodb",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    ).Table("experiments-production")
    dynamo.put_item(Item=experiment_data)
    

    print("uploading to dynamodb samples table...")
    dynamo = boto3.resource(
        "dynamodb",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    ).Table("samples-production")
    dynamo.put_item(Item=samples_data)


    print("uploading R object to s3...")
    s3 = boto3.client(
        "s3",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    )
    bucket, key = FILE_NAME.split("/", 1)
    
    with open("/output/experiment.rds", "rb") as f:
        s3.put_object(Body=f, Bucket=bucket, Key=key)
    
    print("successful. experiment is now accessible at:")
    print(f"https://scp.biomage.net/experiments/{experiment_id}/data-exploration")

main()

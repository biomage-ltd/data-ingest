#!/usr/bin/python3
import hashlib
import os
import pandas
from scipy.io import mmread
import matplotlib.pyplot as plt
import boto3
import json

COLOR_POOL = []

with open("/data-ingest/src/color_pool.json") as f:
    COLOR_POOL = json.load(f)


def calculate_checksum(filenames):
    hash = hashlib.md5()
    for fn in filenames:
        if os.path.isfile(fn):
            hash.update(open(fn, "rb").read())
    return hash.hexdigest()

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
    config_qc = None
    with open("/output/config_qc.json", "r") as f:
        config_qc = json.load(f)

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
    else:
        cellSets = [cell_set, scratchpad]

    print("Experiment name is", config["name"])

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
        "processingConfig": config_qc, 
    }

    print(experiment_data)

    access_key = os.getenv("AWS_ACCESS_KEY_ID")
    secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")
    
    print("uploading to dynamodb...")
    dynamo = boto3.resource(
        "dynamodb",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    ).Table("experiments-production")
    dynamo.put_item(Item=experiment_data)
    
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

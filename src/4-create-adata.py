#!/usr/bin/python3
import hashlib
import os
import anndata
import pandas
import scanpy as sc
from scipy.io import mmread
import boto3
import json

COLOR_POOL = []

with open("/data-ingest/src/color_pool.json") as f:
    COLOR_POOL = json.load(f)


def process_cells():
    df = pandas.read_csv("/output/r-out-cells.csv", names=("cell_names",))
    df.reset_index(inplace=True)

    df.set_index("cell_names", inplace=True)
    df.index.names = [None]

    df.rename(columns={"index": "cell_ids"}, inplace=True)

    return df


def process_genes():
    df = pandas.read_csv(
        "/output/r-out-dispersions.csv",
        names=("gene_names", "dispersions"),
    )

    df.set_index("gene_names", inplace=True, drop=False)
    df.index.names = [None]

    return df


def calculate_checksum(filenames):
    hash = hashlib.md5()
    for fn in filenames:
        if os.path.isfile(fn):
            hash.update(open(fn, "rb").read())
    return hash.hexdigest()


def create_file(checksum):

    print("reading mtx files")
    X = mmread("/output/r-out-normalized.mtx").toarray()
    X_raw = mmread("/output/r-out-raw.mtx").tocsr()

    # create cell and gene matrix
    print("creating obs and var dataframes")

    df = {}
    df["obs"] = process_cells()
    df["var"] = process_genes()

    print("initializing with raw alues")
    # initialize with raw values
    adata = anndata.AnnData(X=X_raw, obs=df["obs"], var=df["var"])
    adata.raw = adata

    # overwrite with processed data
    print("overwriting with processed data")
    adata.X = X

    print("running PCA...")
    sc.tl.pca(adata, svd_solver="arpack")

    print("running louvain")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.louvain(adata)

    print("saving file")
    adata.write("/output/experiment.h5ad")

    return adata


def cell_sets(adata):
    # construct new cell set group
    cell_set = {
        "key": "louvain",
        "name": "Louvain clusters",
        "rootNode": True,
        "children": [],
    }

    raw = adata.obs[["louvain", "cell_ids"]]

    for cluster in raw["louvain"].cat.categories:
        view = raw[raw.louvain == cluster]["cell_ids"]
        cell_set["children"].append(
            {
                "key": f"louvain-{cluster}",
                "name": f"Cluster {cluster}",
                "color": COLOR_POOL.pop(0),
                "cellIds": [int(id) for id in view.tolist()],
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

    adata = create_file(experiment_id)
    cell_set = cell_sets(adata)

    name = os.getenv("EXPERIMENT_NAME")

    print("Experiment name is", name)

    FILE_NAME = f"biomage-source-production/{experiment_id}/python.h5ad"

    experiment_data = {
        "experimentId": experiment_id,
        "experimentName": name,
        "cellSets": [
            cell_set,
            {
                "key": "scratchpad",
                "name": "Scratchpad",
                "rootNode": True,
                "children": [],
            },
        ],
    }

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

    print("uploading Python object to s3...")
    s3 = boto3.client(
        "s3",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_access_key,
        region_name="eu-west-1",
    )
    bucket, key = FILE_NAME.split("/", 1)

    with open("/output/experiment.h5ad", "rb") as f:
        s3.put_object(Body=f, Bucket=bucket, Key=key)

    print("uploading R object to s3...")
    with open("/output/experiment.rds", "rb") as f:
        s3.put_object(Body=f, Bucket=bucket, Key=key.replace("python.h5ad", "r.rds"))

    print("successful. experiment is now accessible at:")
    print(f"https://scp.biomage.net/experiments/{experiment_id}/data-exploration")


main()

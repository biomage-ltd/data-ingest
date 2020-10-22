library(pagoda2)
library(magrittr)
library(parallel)
library(Matrix)
require(data.table)
library(RJSONIO)

create_dataframe <- function(config) {
    data_type <- config$input["type"]
    data <- list()

    if (data_type == "10x") {
        message("Loading 10x data set from input folder.")
        data$raw <- pagoda2::read.10x.matrices("/input")

        # column names are barcodes prefixed with `one_`.
        # Remove as we are processing one dataset.
        colnames(data$raw) <- gsub("^one_", "", colnames(data$raw))

        genes <- read.csv(
            "/input/genes.tsv",
            header = F, sep = "\t", stringsAsFactors = F
        )

        colnames(genes) <- c("id", "name")
        rownames(data$raw) <- genes[, "id"]
    }

    if (data_type == "table") {
        path <- config$input["path"]

        message(paste("Loading table-type data set from", path))

        data$raw <- as.matrix(read.table(paste("/input", path, sep = "/")))
    }

    message(
        paste(
            "Found", nrow(data$raw), "genes and", ncol(data$raw), "cells."
        )
    )

    return(data)
}

prepare_scrublet_table <- function(data) {
    table <- data.table(
        as.matrix(
            t(
                data$filtered
            )
        )
    )

    path <- "/output/pre-doublet-matrix.csv"

    file.create(path)
    data.table::fwrite(table, file = path)
}

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")

message("Creating raw dataframe...")
data <- create_dataframe(config)

message("Filtering cells by size...")
data$filtered <- pagoda2::gene.vs.molecule.cell.filter(
    data$raw,
    min.cell.size = 10e1,
    plot = F
)

message("Filtering cells by molecules/gene...")
data$filtered <- data$filtered[
    Matrix::rowSums(data$filtered > 0) > 5,
]

message("Exporting pre-scrublet data...")
prepare_scrublet_table(data)
saveRDS(data, file = "/output/pre-doublet-data.rds", compress = FALSE)
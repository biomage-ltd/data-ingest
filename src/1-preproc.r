library(pagoda2)
library(magrittr)
library(parallel)
library(Matrix)
require(rjson)
require(data.table)

create_dataframe <- function() {
    # load data
    data <- pagoda2::read.10x.matrices("/input")

    # column names are barcodes prefixed with `one_`.
    # Remove as we are processing one dataset.
    colnames(data) <- gsub("^one_", "", colnames(data))

    # read gene table
    genes <- read.csv(
        "/input/genes.tsv",
        header = F, sep = "\t", stringsAsFactors = F
    )

    # change first two column names to be `id` and `name`
    # assign names to the matrix instead of IDs
    colnames(genes) <- c("id", "name")
    rownames(data) <- genes[, "name"]

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

data <- list()

message("Creating raw dataframe...")
data$raw <- create_dataframe()

message("Filtering cells by size...")
data$filtered <- pagoda2::gene.vs.molecule.cell.filter(
    data$raw,
    min.cell.size = 1e3,
    plot = F
)

message("Filtering cells by molecules/gene...")
data$filtered <- data$filtered[
    Matrix::rowSums(data$filtered > 0) > 5,
]

message("Exporting pre-scrublet data...")
prepare_scrublet_table(data)
saveRDS(data, file = "/output/pre-doublet-data.rds", compress = FALSE)
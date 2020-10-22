library(pagoda2)
library(conos)
library(Matrix)
require(data.table)
require(conos)
library(gprofiler2)

get_real_cells <- function(data, score_filter = 0.2) {
    scores <-
        data.table::fread(
            "/output/doublet-scores.csv",
            col.names = c("score")
        )

    scores <- scores[, "barcodes" := colnames(data$filtered)]

    barcodes <- scores[score <= 0.2]$barcodes

    return(barcodes)
}

message("reloading old matrices...")
data <- readRDS("/output/pre-doublet-data.rds")

message("getting scrublet results...")
real_barcodes <- get_real_cells(data)

message("removing doublets...")
data$filtered <- data$filtered[, real_barcodes]

message('finding genome annotations for genes...')
config <- RJSONIO::fromJSON("/input/meta.json")
organism <- config$organism

annotations <- gprofiler2::gconvert(
    query = rownames(data$filtered), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE
)

message("creating pagoda2 object...")
data <- Pagoda2$new(data$filtered)
data$misc$gene_annotations <- annotations

message("variance normalization...")
data$adjustVariance(plot = F, gam.k = 10)

# we're taking 3000 top overdispersed genes
# (a lot more than the number that is significantly overdispersed)
message("computing PCA reduction...")
data$calculatePcaReduction(nPcs = 50, n.odgenes = 3e3)

message("creating kNN graph...")
data$makeKnnGraph(
    k = 20, type = "PCA", center = T,
    distance = "cosine", weight.type = "none"
)

message("computing leiden clusters...")
data$getKnnClusters(method = conos::leiden.community, type = "PCA")

message("Saving gene and cell data...")
write.table(
    data$misc$varinfo[c("qv")],
    file = "/output/r-out-dispersions.csv",
    sep = ",",
    col.names = F
)

write.table(
    rownames(data$counts),
    file = "/output/r-out-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

write.table(
    data$misc$gene_annotations,
    file = "/output/r-out-annotations.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

message("saving R object...")
saveRDS(data, file = "/output/experiment.rds", compress = FALSE)

message("saving normalized matrix...")
Matrix::writeMM(data$counts, file = "/output/r-out-normalized.mtx")

message("saving raw matrix...")
Matrix::writeMM(
    data$misc$rawCounts[
        rownames(data$counts),
        colnames(data$counts)
    ],
    file = "/output/r-out-raw.mtx"
)
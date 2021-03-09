library(Seurat)
library(Matrix)
require(data.table)
library(gprofiler2)

set.seed(123)
options(future.globals.maxSize= 1000 * 1024 ^ 2)
source("src/help.r")


################################################
## LOADING SPARSE MATRIX AND CONFIGURATION
################################################

message("reloading old matrices...")
scdata <- readRDS("/output/pre-doublet-data.rds")

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")
metadata <- check_config(scdata, config)

message("Creating Seurat Object...")
scdata <- Seurat::CreateSeuratObject(scdata$filtered, assay='RNA', min.cells=3, min.features=200, meta.data=metadata)

################################################
## GETTING METADATA AND ANNOTATION
################################################

message('finding genome annotations for genes...')
organism <- config$organism

annotations <- gprofiler2::gconvert(
    query = rownames(scdata), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

if(organism%in%c("hsapiens", "mmusculus")){
  message("Adding MT information...")
  mt.features <-  annotations$input[grep("^mt-", annotations$name, ignore.case = T)]
  scdata <- PercentageFeatureSet(scdata, features=mt.features , col.name = "percent.mt")
}

message("getting scrublet results...")
scores <- get_doublet_score(scdata)

message("Adding doublet scores information...")
idt <- scores$barcodes[scores$barcodes%in%rownames(scdata@meta.data)]
scdata@meta.data[idt, "doublet_scores"] <- scores[idt, "score"]


################################################
## DATA PROCESSING
################################################

#[HARDCODED]

config.cellSizeDistribution <- list()
config.mitochondrialContent <- list()
config.classifier <- list()
config.numGenesVsNumUmis <- list()
config.doubletScores <- list()

#
# Step 1: Cell size distribution filter
#




message("Normalization step...")

####### Default  configuration settings

nfeatures <- 2e3
normalization_method <- "LogNormalize"
normalization_scale_factor<- 10000
pca_nPCs <- 30
neighbors_metric <- "cosine"
clustering_method <- 1 #"Louvain"
clustering_resolution <- 0.5
umap_min_distance <- 0.3
umap_distance_metric <- "euclidean"
tsne_perplexity <- min(30, ncol(scdata)/100)
tsne_learning_rate <- max(200, ncol(scdata)/12)

if(as.logical(config$samples[["multisample"]])){
    
    # Seurat V3 pipeline (see for other methods: https://satijalab.org/seurat/archive/v3.0/integration.html)
    data.split <- SplitObject(scdata, split.by = "type")
    for (i in 1:length(data.split)) {
        data.split[[i]] <- NormalizeData(data.split[[i]], normalization.method = normalization_method, scale.factor = normalization_scale_factor, verbose = F)
        data.split[[i]] <- FindVariableFeatures(data.split[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    }
    
    data.anchors <- FindIntegrationAnchors(object.list = data.split, dims = 1:pca_nPCs, verbose = FALSE)
    scdata <- IntegrateData(anchorset = data.anchors, dims = 1:pca_nPCs)
    DefaultAssay(scdata) <- "integrated"
    scdata <- Seurat::ScaleData(scdata, verbose = F)
    
    scdata <- FindVariableFeatures(scdata, selection.method = "vst", assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
    vars <- HVFInfo(object = scdata, assay = "RNA", selection.method = 'vst')
    
}else{
    scdata <- Seurat::NormalizeData(scdata, normalization.method = normalization_method, scale.factor = normalization_scale_factor, verbose = F)
    scdata <-Seurat::FindVariableFeatures(scdata, selection.method = "vst", nfeatures = nfeatures, verbose = F)
    scdata<-Seurat::ScaleData(scdata, verbose = F)
    
    vars <- HVFInfo(object = scdata, assay = "RNA", selection.method = 'vst')
    # In case of SCTransform
    #vars <- HVFInfo(object = scdata, selection.method = 'sctransform')
    #vars <- vars[, "residual_variance", drop = FALSE]

}

message("computing PCA reduction...")
scdata<-Seurat::RunPCA(scdata, npcs = 50, features = VariableFeatures(object=scdata), verbose=FALSE)

message("creating kNN graph...")
scdata <- FindNeighbors(scdata, k.param = 20, annoy.metric = neighbors_metric, verbose=FALSE) #default method

message("computing louvain clusters...")
scdata <- FindClusters(scdata, resolution=clustering_resolution, verbose = FALSE, algorithm = clustering_method) #default method (Louvain clusters)

message("Running embedding")
scdata <- RunUMAP(scdata, reduction='pca', dims = 1:pca_nPCs, verbose = F, umap.method = "uwot-learn", min.dist = umap_min_distance, metric = umap_distance_metric)
scdata <- RunTSNE(scdata, reduction = 'pca', dims = 1:pca_nPCs, perplexity = tsne_perplexity, learning.rate = tsne_learning_rate)

## Adding more information to misc embedding

scdata@misc$embedding_configuration <- list()

scdata@misc$embedding_configuration[["UMAP"]] <- list()
scdata@misc$embedding_configuration$UMAP["pca_nPCs"] <- pca_nPCs
scdata@misc$embedding_configuration$UMAP["umap_min_distance"] <- umap_min_distance
scdata@misc$embedding_configuration$UMAP["umap_distance_metric"] <- umap_distance_metric

scdata@misc$embedding_configuration[["TSNE"]] <- list()
scdata@misc$embedding_configuration$TSNE["pca_nPCs"] <- pca_nPCs
scdata@misc$embedding_configuration$TSNE["tsne_perplexity"] <- tsne_perplexity
scdata@misc$embedding_configuration$TSNE["tsne_learning_rate"] <- tsne_learning_rate

message("Storing gene annotations...")
scdata@misc[["gene_annotations"]] <- annotations

message("Storing cells id...")
# Keeping old version of ids starting from 0
scdata$cells_id <- 0:(nrow(scdata@meta.data)-1)

message("Storing dispersion...")
# Convert to Gene Symbol
vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
vars$ENSEMBL <- rownames(vars)
scdata@misc[["gene_dispersion"]] <- vars

pdf("/output/umap.pdf")
DimPlot(scdata, reduction = "umap")
dev.off()


if(as.logical(config$samples$multisample)){
    pdf("/output/umap_type.pdf")
    DimPlot(scdata, reduction = "umap", group.by = "type")
    dev.off()
}

message("saving R object...")
saveRDS(scdata, file = "/output/experiment.rds", compress = FALSE)

message("saving cluster info...")
write.table(
    data.frame(Cells_ID = scdata$cells_id[names(scdata@active.ident)], Clusters=scdata@active.ident),
    file = "/output/cluster-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

if(as.logical(config$samples$multisample)){
    message("saving multsiample info...")
    write.table(
        data.frame(Cells_ID = scdata$cells_id[names(scdata@active.ident)], type=scdata$type),
        file = "/output/multisample-cells.csv",
        quote = F, col.names = F, row.names = F,
        sep = "\t"
    )
}

vars <- vars[rownames(scdata), c("ENSEMBL", "variance.standardized")]
message("Saving gene and cell data...")
write.table(
    vars,
    file = "/output/r-out-dispersions.csv",
    sep = ",",
    col.names = F, row.names = F
)

write.table(
    colnames(scdata),
    file = "/output/r-out-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

write.table(
    scdata@misc[["gene_annotations"]][scdata@misc[["gene_annotations"]]$input%in%rownames(scdata), ],
    file = "/output/r-out-annotations.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

message("saving normalized matrix...")
Matrix::writeMM(t(scdata@assays[[scdata@active.assay]]@data
                  ), file = "/output/r-out-normalized.mtx")

message("saving raw matrix...")
Matrix::writeMM(t(
    scdata@assays[["RNA"]]@counts[
        rownames(scdata),
        colnames(scdata)
    ]),
    file = "/output/r-out-raw.mtx", verb
)

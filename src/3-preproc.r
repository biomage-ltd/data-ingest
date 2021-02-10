library(Seurat)
library(Matrix)
require(data.table)
library(gprofiler2)

set.seed(123)
options(future.globals.maxSize= 1000 * 1024 ^ 2)

# get_doublet_score function 
#' @description Get the cells with its doublet scores computed previously through scrublet
#' @param data matrix with barcodes as columns
#' 
#' @export save barcodes and double scores

get_doublet_score <- function(data) {
    scores <-
        data.table::fread(
            "/output/doublet-scores.csv",
            col.names = c("score")
        )

    scores <- as.data.frame(scores[, "barcodes" := colnames(data$filtered)])
    rownames(scores) <- scores$barcodes    
    return(scores)
}

# check_config function 
#' @description Create metadata dataframe from config files
#' @param data matrix with barcodes as columns to assign metadata information
#' @param config config list from meta.json
#' 
#' @export save barcodes to keep

check_config <- function(data, config){
    metadata <- NULL
    
    if("type" %in% names(config$samples$samples_info)){
        metadata <- data.frame(row.names = colnames(data$filtered))
        metadata[colnames(data$filtered), "type"] <- unlist(lapply(strsplit(colnames(data$filtered), "_"), `[`, 1))
        
        rest_metadata <- as.data.frame(config$samples$samples_info)
        for(var in colnames(rest_metadata)[-which(colnames(rest_metadata)%in%"type")]){
            metadata[, var] <- rest_metadata[, var][match(metadata$type, rest_metadata$type)]
        }
    }
    
    return(metadata)
}

message("reloading old matrices...")
data <- readRDS("/output/pre-doublet-data.rds")

message("getting scrublet results...")
scores <- get_doublet_score(data)

# message("removing doublets...")
# data$filtered <- data$filtered[, scores$barcodes[scores$score<=0.2]]

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")
metadata <- check_config(data, config)

message("Creating Seurat Object...")
data <- Seurat::CreateSeuratObject(data$filtered, assay='RNA', min.cells=3, min.features=200, meta.data=metadata)

message('finding genome annotations for genes...')
organism <- config$organism

annotations <- gprofiler2::gconvert(
    query = rownames(data), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

message("Adding MT information...")
data_gene_symbol <- data
rownames(data_gene_symbol@assays$RNA@data) <- annotations$name[match(rownames(data), annotations$input)]
rownames(data_gene_symbol@assays$RNA@counts) <- annotations$name[match(rownames(data), annotations$input)]
data_gene_symbol <- PercentageFeatureSet(data_gene_symbol, pattern = "^MT-", col.name = "percent.mt")
data@meta.data$percent.mt <- data_gene_symbol@meta.data[rownames(data@meta.data), "percent.mt"]
rm(data_gene_symbol)

message("Adding doublet scores information...")
idt <- scores$barcodes[scores$barcodes%in%rownames(data@meta.data)]
data@meta.data[idt, "doublet_scores"] <- scores[idt, "score"]

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
tsne_perplexity <- min(30, ncol(data)/100)
tsne_learning_rate <- max(200, ncol(data)/12)

if(as.logical(config$samples[["multisample"]])){
    
    # Seurat V3 pipeline (see for other methods: https://satijalab.org/seurat/archive/v3.0/integration.html)
    data.split <- SplitObject(data, split.by = "type")
    for (i in 1:length(data.split)) {
        data.split[[i]] <- NormalizeData(data.split[[i]], normalization.method = normalization_method, scale.factor = normalization_scale_factor, verbose = F)
        data.split[[i]] <- FindVariableFeatures(data.split[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    }
    
    data.anchors <- FindIntegrationAnchors(object.list = data.split, dims = 1:pca_nPCs, verbose = FALSE)
    data <- IntegrateData(anchorset = data.anchors, dims = 1:pca_nPCs)
    DefaultAssay(data) <- "integrated"
    data <- Seurat::ScaleData(data, verbose = F)
    
    data <- FindVariableFeatures(data, selection.method = "vst", assay = "RNA", nfeatures = nfeatures, verbose = FALSE)
    vars <- HVFInfo(object = data, assay = "RNA", selection.method = 'vst')
    
}else{
    data <- Seurat::NormalizeData(data, normalization.method = normalization_method, scale.factor = normalization_scale_factor, verbose = F)
    data <-Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = nfeatures, verbose = F)
    data<-Seurat::ScaleData(data, verbose = F)
    
    vars <- HVFInfo(object = data, assay = "RNA", selection.method = 'vst')
}

vars <- vars[, 'variance.standardized', drop=FALSE]
# In case of SCTransform
#vars <- HVFInfo(object = data, selection.method = 'sctransform')
#vars <- vars[, "residual_variance", drop = FALSE]

message("Storing annotations...")
data@misc$gene_annotations <- annotations

message("computing PCA reduction...")
data<-Seurat::RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)

message("creating kNN graph...")
data <- FindNeighbors(data, k.param = 20, annoy.metric = neighbors_metric, verbose=FALSE) #default method

message("computing louvain clusters...")
data <- FindClusters(data, resolution=clustering_resolution, verbose = FALSE, algorithm = clustering_method) #default method (Louvain clusters)

message("Running embedding")
data <- RunUMAP(data, reduction='pca', dims = 1:pca_nPCs, verbose = F, umap.method = "uwot-learn", min.dist = umap_min_distance, metric = umap_distance_metric)
data <- RunTSNE(data, reduction = 'pca', dims = 1:pca_nPCs, perplexity = tsne_perplexity, learning.rate = tsne_learning_rate)

## Adding more information to misc embedding

data@misc[["embedding_configuration"]] <- list()

data@misc$embedding_configuration[["UMAP"]] <- list()
data@misc$embedding_configuration$UMAP["pca_nPCs"] <- pca_nPCs
data@misc$embedding_configuration$UMAP["umap_min_distance"] <- umap_min_distance
data@misc$embedding_configuration$UMAP["umap_distance_metric"] <- umap_distance_metric

data@misc$embedding_configuration[["TSNE"]] <- list()
data@misc$embedding_configuration$TSNE["pca_nPCs"] <- pca_nPCs
data@misc$embedding_configuration$TSNE["tsne_perplexity"] <- tsne_perplexity
data@misc$embedding_configuration$TSNE["tsne_learning_rate"] <- tsne_learning_rate


pdf("/output/umap.pdf")
DimPlot(data, reduction = "umap")
dev.off()


if(as.logical(config$samples$multisample)){
    pdf("/output/umap_type.pdf")
    DimPlot(data, reduction = "umap", group.by = "type")
    dev.off()
}

message("saving R object...")
saveRDS(data, file = "/output/experiment.rds", compress = FALSE)

vars$genes <- rownames(vars)
vars <- vars[, c(2, 1)]
vars <- vars[rownames(data), ]
message("Saving gene and cell data...")
write.table(
    vars,
    file = "/output/r-out-dispersions.csv",
    sep = ",",
    col.names = F, row.names = F
)

write.table(
    colnames(data),
    file = "/output/r-out-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

write.table(
    data@misc[["gene_annotations"]][data@misc[["gene_annotations"]]$input%in%rownames(data), ],
    file = "/output/r-out-annotations.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

message("saving R object...")
saveRDS(data, file = "/output/experiment.rds", compress = FALSE)

message("saving normalized matrix...")
Matrix::writeMM(t(data@assays[[data@active.assay]]@data
                  ), file = "/output/r-out-normalized.mtx")

message("saving raw matrix...")
Matrix::writeMM(t(
    data@assays[["RNA"]]@counts[
        rownames(data),
        colnames(data)
    ]),
    file = "/output/r-out-raw.mtx", verb
)

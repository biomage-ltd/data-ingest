library(Seurat)
library(Matrix)
require(data.table)
library(gprofiler2)

# get_real_cells function 
#' @description Get the cells that pass the filter of doublet scores
#' @param data matrix to filter with barcodes as columns
#' 
#' @export save barcodes to keep

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
real_barcodes <- get_real_cells(data)

message("removing doublets...")
data$filtered <- data$filtered[, real_barcodes]

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")
metadata <- check_config(data, config)

message("Creating Seurat Object...")
data <- Seurat::CreateSeuratObject(data$filtered, assay='RNA', min.cells=3, min.features=200, meta.data=metadata)

#Below are QC plots which are exported as pdf files. Uncomment to look at the qc plots.

#pdf("/output/violin-pl.pdf")
#VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", ncol=2))
#dev.off()
#pdf("/output/scatter-pl.pdf")
#FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#dev.off()

#Using Seurat we can subset data. We can subset based on nUMI, nFeature and mito genes. We have filtered in 1.preproc.r so subset was not used. 

#In order to determine the filter parameters for mito genes we would need to plot them on either a violin plot or a scatter plot. Please refer to the plots above.

#Remove by subsetting... or look under ScaleData below for regressing out mito genes.
#data <- subset(data, subset = percent.mt < 5)

message('finding genome annotations for genes...')
organism <- config$organism

annotations <- gprofiler2::gconvert(
    query = rownames(data), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

message("Adding MT information...")
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")

message("Normalization step...")
if(as.logical(config$samples[["multisample"]])){
    data.split <- SplitObject(data, split.by = "type")
    for (i in 1:length(data.split)) {
        data.split[[i]] <- NormalizeData(data.split[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
        data.split[[i]] <- FindVariableFeatures(data.split[[i]], selection.method = "vst", nfeatures = 3e3, verbose = FALSE)
    }
    
    data.anchors <- FindIntegrationAnchors(object.list = data.split, dims = 1:30,verbose = FALSE)
    data <- IntegrateData(anchorset = data.anchors, dims = 1:30)
    DefaultAssay(data) <- "integrated"
    data <- Seurat::ScaleData(data, verbose = F)
    
    data <- FindVariableFeatures(data, selection.method = "vst", assay = "RNA", nfeatures = 3e3, verbose = FALSE)
    vars <- HVFInfo(object = data, assay = "RNA", selection.method = 'vst')
    
}else{
    data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
    data <-Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 3e3, verbose = F)
    data<-Seurat::ScaleData(data, verbose = F)
    
    vars <- HVFInfo(object = data, assay = "RNA", selection.method = 'vst')
}

vars <- vars[, 'variance.standardized', drop=FALSE]

message("Storing annotations...")
data@misc$gene_annotations <- annotations


#we can additionally regress out the mitochondrial genes
#data <-Seurat::ScaleData(data, vars.to.regress= "percent.mt", verbose=FALSE)

# Or we can use SCTransform (it uses regularized negative binomial regression)
#data <- SCTransform(data, variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = FALSE, do.scale=T)


message("computing PCA reduction...")
data<-Seurat::RunPCA(data, npcs = 50, features = VariableFeatures(object=data), verbose=FALSE)

message("creating kNN graph...")
data <- FindNeighbors(data, k.param = 20, annoy.metric = "cosine", verbose=FALSE) #default method

message("computing louvain clusters...")
data <- FindClusters(data, resolution=0.5, verbose = FALSE) #default method (Louvain clusters)


# In case of SCTransform
#vars <- HVFInfo(object = data, selection.method = 'sctransform')
#vars <- vars[, "residual_variance", drop = FALSE]

message("Running embedding")
data <- RunUMAP(data, reduction='pca', dims = 1:10, verbose = F, umap.method = "uwot-learn")

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

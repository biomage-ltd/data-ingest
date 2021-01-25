library(Seurat)
library(Matrix)
require(data.table)
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

message("Creating Seurat Object...")
data <- Seurat::CreateSeuratObject(data$filtered, assay='RNA', min.cells=3, min.features=200)

#Below are QC plots which are exported as pdf files. Uncomment to look at the qc plots.

#pdf("/output/violin-pl.pdf")
#VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", ncol=2))
#dev.off()
#pdf("/output/scatter-pl.pdf")
#FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#dev.off()

#Using Seurat we can subset data. We can subset based on nUMI, nFeature and mito genes. We have filtered in 1.preproc.r so subset was not used. 

#Uncomment the code below to remove mito genes. 
#message("Subset and remove mito genes...")
#data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

#In order to determine the filter parameters for mito genes we would need to plot them on either a violin plot or a scatter plot. Please refer to the plots above.

#Remove by subsetting... or look under ScaleData below for regressing out mito genes.
#data <- subset(data, subset = percent.mt < 5)

message('finding genome annotations for genes...')
config <- RJSONIO::fromJSON("/input/meta.json")
organism <- config$organism

annotations <- gprofiler2::gconvert(
    query = rownames(data), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

message("Storing annotations...")
#Misc(data, slot="annotations") <- annotations
data@misc$gene_annotations <- annotations

message("Adding MT information...")
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")


message("Normalization step...")
data <- Seurat::NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
data <-Seurat::FindVariableFeatures(data, selection.method = "vst", nfeatures = 3e3, verbose = F)
data<-Seurat::ScaleData(data, verbose = F)
#we can additionally regress out the mitochondrial genes
#data <-Seurat::ScaleData(data, vars.to.regress= "percent.mt", verbose=FALSE)

# Or we can use SCTransform (it uses regularized negative binomial regression)
#data <- SCTransform(data, variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = FALSE, do.scale=T)

vars <- HVFInfo(object = data, assay = "RNA", selection.method = 'vst')
vars <- vars[, 'variance.standardized', drop=FALSE]

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

#We can also create a tsne plot
#message("Running embedding")
#data <- RunTSNE(data, reduction='pca', dims = 1:10)
#pdf("/output/tsne.pdf")
#DimPlot(data, reduction = "tsne")
#dev.off()

message("saving R object...")
saveRDS(data, file = "/output/experiment.rds", compress = FALSE)

message("Saving gene and cell data...")
write.table(
    vars,
    file = "/output/r-out-dispersions.csv",
    sep = ",",
    col.names = F
)

write.table(
    colnames(data),
    file = "/output/r-out-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

write.table(
    data@misc[["gene_annotations"]],
    file = "/output/r-out-annotations.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

message("saving R object...")
saveRDS(data, file = "/output/experiment.rds", compress = FALSE)

message("saving normalized matrix...")
Matrix::writeMM(t(data@assays$RNA@data), file = "/output/r-out-normalized.mtx")

message("saving raw matrix...")
Matrix::writeMM(t(
    data@assays$RNA@counts[
        rownames(data@assays$RNA@data),
        colnames(data@assays$RNA@data)
    ]),
    file = "/output/r-out-raw.mtx", verb
)

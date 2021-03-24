library(Seurat)
library(Matrix)
library(dplyr)
require(data.table)
library(gprofiler2)

set.seed(123)
options(future.globals.maxSize= 1000 * 1024 ^ 2)
source("/data-ingest/src/help.r")
source("/data-ingest/src/QC_helpers/cellSizeDistribution.r")
source("/data-ingest/src/QC_helpers/mitochondrialContent.r")
source("/data-ingest/src/QC_helpers/classifier.r")
source("/data-ingest/src/QC_helpers/numGenesVsNumUmis.r")
source("/data-ingest/src/QC_helpers/doubletScores.r")
source("/data-ingest/src/QC_helpers/dataIntegration.r")
source("/data-ingest/src/QC_helpers/computeEmbedding.r")


################################################
## LOADING SPARSE MATRIX AND CONFIGURATION
################################################

message("reloading old matrices...")
scdata <- readRDS("/output/pre-doublet-data.rds")

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")
metadata <- check_config(scdata, config)

message("Creating Seurat Object...")
seurat_obj <- Seurat::CreateSeuratObject(scdata$filtered, assay='RNA', min.cells=1, min.features=1, meta.data=metadata)

################################################
## GETTING METADATA AND ANNOTATION
################################################

message('finding genome annotations for genes...')
organism <- config$organism

annotations <- gprofiler2::gconvert(
    query = rownames(seurat_obj), organism = organism, target="ENSG", mthreshold = Inf, filter_na = FALSE)

if(organism%in%c("hsapiens", "mmusculus")){
  message("Adding MT information...")
  mt.features <-  annotations$input[grep("^mt-", annotations$name, ignore.case = T)]
  seurat_obj <- PercentageFeatureSet(seurat_obj, features=mt.features , col.name = "percent.mt")
}

message("getting scrublet results...")
scores <- get_doublet_score(seurat_obj)

pdf("/output/doublet_scores_hist.pdf")
hist(scores$score, breaks = 100)
dev.off()

pdf("/output/mitochondrialFractionLogHistogram.pdf")
hist(seurat_obj$fraction.mt, breaks=200)
dev.off()

pdf("/output/mitochondrialFractionLogScatter.pdf")
plot(seurat_obj$nCount_RNA, seurat_obj$fraction.mt)
dev.off()


pdf("/output/UMI_hist.pdf",width = 14, height = 7)
par(mfcol=c(1,2),cex.lab=2)
hist(colSums(seurat_obj@assays$RNA@data), breaks=200)
hist(rowSums(seurat_obj@assays$RNA@data), breaks=200)
dev.off()

genes_umi <- rowSums(seurat_obj@assays$RNA@data)
counts <- colSums(seurat_obj@assays$RNA@data)


pdf("/output/nCount_vs_mito.fr_hist.pdf")
plot(seurat_obj@meta.data$nCount_RNA, seurat_obj@meta.data$fraction.mt)
dev.off()

pdf("/output/GenesVsNumUmis.pdf")
plot(seurat_obj@meta.data$nCount_RNA, seurat_obj@meta.data$nFeature_RNA)
dev.off()

message("Adding doublet scores information...")
idt <- scores$barcodes[scores$barcodes%in%rownames(seurat_obj@meta.data)]
seurat_obj@meta.data[idt, "doublet_scores"] <- scores[idt, "score"]

file_ed <-  "/output/pre-emptydrops-data.rds"
if (file.exists(file_ed)) {
  seurat_obj@tools$flag_filtered <- FALSE
  message("getting emptyDrops results...")
  emptydrops_out <- readRDS(file = "/output/pre-emptydrops-data.rds")

  emptydrops_out_df <- emptydrops_out %>%
    as.data.frame() %>%
    rlang::set_names(~ paste0("emptyDrops_", .)) %>%
    tibble::rownames_to_column("barcode")

  # adding emptydrops data to meta.data for later filtering (using left join)
  meta.data <- seurat_obj@meta.data  %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(emptydrops_out_df)
  rownames(meta.data) <- meta.data$barcode


  # some dignostic plots for empty-drops
  pdf("/output/emptyDrops_Total_hist.pdf",width = 16, height = 14)
  par(mfcol=c(2,2),cex.lab=2)
  # hist(log10(emptydrops_out_df$emptyDrops_Total), breaks = 300)
  hist((emptydrops_out_df$emptyDrops_PValue), breaks = 300)
  hist((emptydrops_out_df$emptyDrops_FDR), breaks = 300)
  plot1_data_1  <- log10(meta.data$emptyDrops_FDR)
  names(plot1_data_1) <- rep("FDR", length(plot1_data_1))
  plot1_data_2  <- log10(meta.data$nCount_RNA)
  names(plot1_data_2) <- rep("log_u", length(plot1_data_2))
  plot(log10(meta.data$nCount_RNA),log10(meta.data$emptyDrops_FDR))
  plot((meta.data$nCount_RNA),(meta.data$emptyDrops_FDR))
  dev.off()


  message("Adding emptyDrops scores information...")
  seurat_obj@meta.data <- meta.data
  # previously (before joining into meta.data), results were just dumped as a additional slot
  # leaving the code here in case bugs arise from above solution
  # seurat_obj@tools$CalculateEmptyDrops <- emptydrops_out
} else {
  # Later on, when creating the config file, the enable will look the value of flag_filtered to deactivate the classifier filter
  message("emptyDrops results not present, skipping...")
  seurat_obj@meta.data$emptyDrops_FDR <- NA
  seurat_obj@tools$flag_filtered <- TRUE
}

sink("/output/routput.Rout")
print("UMI counts of highest expressed genes")
tail(sort(genes_umi), n=30)
print("median rowsums")
median(rowSums(seurat_obj@assays$RNA@data))
print("empty-drops table of FDR threshold categories (# UMIs for a given threshold interval")
table(cut(emptydrops_out_df$emptyDrops_FDR, breaks = c(-Inf,0,0.0001,0.01,0.1,0.5,1)), useNA="ifany")
print("empty-drops table of FDR threshold categories (# UMIs for a given threshold interval")
table(cut(emptydrops_out_df$emptyDrops_PValue, breaks = c(-Inf,0,0.0001,0.01,0.1,0.5,1)), useNA="ifany")
sink()


################################################
## DATA PROCESSING
################################################

#[HARDCODED]
config.cellSizeDistribution <- list(enabled="true", 
    auto="true", 
    filterSettings = list(minCellSize=1080, binStep = 200)
)

config.mitochondrialContent <- list(enabled="true", auto="true", 
    filterSettings = list(method="absolute_threshold", methodSettings = list(
        absolute_threshold=list(maxFraction=0.1, binStep=0.05)
        )
    )
)

config.classifier <- list(enabled=tolower(as.character(!seurat_obj@tools$flag_filtered)) # emptyDrops results not present
    , auto="true", 
    filterSettings = list(FDR=0.1)
)

config.numGenesVsNumUmis <- list(enabled="true", auto="true", 
    filterSettings = list(regressionType = "gam", regressionTypeSettings = list(
        "gam" = list(p.level=0.001)
        )
    )
)

config.doubletScores <- list(enabled="true", auto="true", 
    filterSettings = list(probabilityThreshold = 0.25, binStep = 0.05)
)

# BE CAREFUL! The method is based on config.json. For multisample only seuratv3, for unisample LogNormalize
identified.method <- ifelse(as.logical(config$samples$multisample), "seuratv3", "unisample")
config.dataIntegration <- list(enabled="true", auto="true", 
    dataIntegration = list( method = identified.method , 
                        methodSettings = list(seuratv3=list(numGenes=2000, normalisation="LogNormalize"), 
                                            unisample=list(numGenes=2000, normalisation="LogNormalize"))),
    dimensionalityReduction = list(method = "rpca", numPCs = 30, excludeGeneCategories = c())
)

config.computeEmbedding <- list(enabled="true", auto="true", 
    embeddingSettings = list(method = "umap", methodSettings = list(
                                umap = list(minimumDistance=0.3, distanceMetric="euclidean"), 
                                tsne = list(perplexity=30, learningRate=200)
                            ) 
                        ), 
    clusteringSettings = list(method = "louvain", methodSettings = list(
                            louvain = list(resolution = 0.5)
                            )
                        )
)



#
# Step 1: Cell size distribution filter
#

# message("Filter 1")
# result.step1 <- cellSizeDistribution(seurat_obj, config.cellSizeDistribution)
## Update config
# config.cellSizeDistribution <- result.step1$config
# result.step1$config$filterSettings$minCellSize
## str(result.step1$plotData)
## List of 2
##  $ plot1: Named num [1:11217] 3483 6019 3892 3729 4734 ...
##   ..- attr(*, "names")= chr [1:11217] "u" "u" "u" "u" ...
##  $ plot2:List of 2
##   ..$ : Named num [1:11217] 3483 6019 3892 3729 4734 ...
##   .. ..- attr(*, "names")= chr [1:11217] "u" "u" "u" "u" ...
##   ..$ : Named int [1:11217] 6807 1276 1894 4438 16 6867 887 3494 10873 4161 ...
##   .. ..- attr(*, "names")= chr [1:11217] "rank" "rank" "rank" "rank" ...

#
# Step 2: Mitochondrial content filter
#

# message("Filter 2")
## Be aware that all the currents experiment does not have the slot fracion.mt, but they have the percent.mt. 
## To be consistent, in this new version I have transformed to fracion.mt [Line 46]
# result.step2 <- mitochondrialContent(result.step1$data, config.mitochondrialContent)
## Update config
# config.mitochondrialContent <- result.step2$config
# result.step2$config$filterSettings

## plotData plots
## Plot 1 (histogram with fraction MT)
# h=hist(result.step2$plotData$plot1,plot=FALSE)
# h$density = h$counts/sum(h$counts)
# plot(h,freq=FALSE)
# abline(v=result.step2$config$filterSettings$methodSettings$absolute_threshold$maxFraction, col = "red")
#
## Plot 2 (scatter plot)
# result.step1$data$Valid.cells.MT <- result.step1$data$fraction.MT < result.step2$config$filterSettings$methodSettings$absolute_threshold$maxFraction
# FeatureScatter(result.step1$data, "nCount_RNA", "fraction.mt", group.by = "Valid.cells.MT")
## or direcly plot with the plotData results
# plot(result.step2$plotData$plot2$u, result.step2$plotData$plot2$`MT-content`, xlab = "UMIs", ylab = "MT-fraction")

#
# Step 3: Classifier filter
#

# message("Filter 3")
# Waiting filter 3
# result.step3 <- classifier(result.step2$data, config.classifier)
# str(result.step3$plotData)
## Update config
# config.classifier <- result.step3$config

#
# Step 4: Number of genes vs number of UMIs filter
#

# message("Filter 4")
# result.step4 <- numGenesVsNumUmis(result.step3$data, config.numGenesVsNumUmis)
## Update config
# config.numGenesVsNumUmis <- result.step4$config

## plotData plots
## Plot 1 (scatter plot with bands)
# plot(result.step4$plotData$plot1$log10_UMIs, result.step4$plotData$plot1$log10_genes, ylab = "log10_genes", xlab = "log10_UMIs")
# lines(result.step4$plotData$plot1$log10_UMIs, result.step4$plotData$plot1$upper_cutoff, col = "red")
# lines(result.step4$plotData$plot1$log10_UMIs, result.step4$plotData$plot1$lower_cutoff, col = "red")


#
# Step 5: Doublet scores filter
#

# message("Filter 5")
# result.step5 <- doubletScores(result.step4$data, config.doubletScores)
## Update config
# config.doubletScores <- result.step5$config
## plotData plots
## Plot 1 (histogram with fraction MT)
# h=hist(result.step5$plotData$plot1,plot=FALSE)
# h$density = h$counts/sum(h$counts)
# plot(h,freq=FALSE)


#
# Step 6: Data integration
#

message("Filter 6")
result.step6 <- dataIntegration(seurat_obj, config.dataIntegration)
# Update config
config.dataIntegration <- result.step6$config

#
# Step 7: Compute embedding
#

message("Filter 7")
result.step7 <- computeEmbedding(result.step6$data, config.computeEmbedding)
# Update config
config.computeEmbedding <- result.step7$config


seurat_obj <- result.step7$data

message("Storing gene annotations...")
seurat_obj@misc[["gene_annotations"]] <- annotations

message("Storing cells id...")
# Keeping old version of ids starting from 0
seurat_obj$cells_id <- 0:(nrow(seurat_obj@meta.data)-1)

message("Storing dispersion...")
# Convert to Gene Symbol
# Following the answer in this issue (https://github.com/satijalab/seurat/issues/2778) for Seurat V3, FindVariableFeatures
# does not support for multisample. As a solution we will recompute FindVariables with RNA assays like it is a unisample experiment. 
# HARDCODE: nfeature to 2000 (default value of the function)
nfeautes <- 2000
seurat_obj_dispersion <- FindVariableFeatures(seurat_obj, selection.method = "vst", assay = "RNA", nfeatures = nfeautes, verbose = FALSE)
vars <- HVFInfo(object = seurat_obj_dispersion, assay = "RNA", selection.method = 'vst') # to create vars
vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
vars$ENSEMBL <- rownames(vars)
seurat_obj@misc[["gene_dispersion"]] <- vars


message("Storing color pool...")
# We store the color pool in a slot in order to be able to access it during computeEmbedding
color_pool <- RJSONIO::fromJSON("/data-ingest/src/color_pool.json")
seurat_obj@misc[["color_pool"]] <- color_pool


pdf("/output/umap.pdf")
DimPlot(seurat_obj, reduction = "umap")
dev.off()

pdf("/output/pca.pdf")
DimPlot(seurat_obj, reduction = "pca")
dev.off()


if(as.logical(config$samples$multisample)){
    pdf("/output/umap_type.pdf")
    DimPlot(seurat_obj, reduction = "umap", group.by = "type")
    dev.off()
}

message("saving R object...")
saveRDS(seurat_obj, file = "/output/experiment.rds", compress = FALSE)

message("saving cluster info...")
write.table(
    data.frame(Cells_ID = seurat_obj$cells_id[names(seurat_obj@active.ident)], Clusters=seurat_obj@active.ident),
    file = "/output/cluster-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

if(as.logical(config$samples$multisample)){
    message("saving multsiample info...")
    write.table(
        data.frame(Cells_ID = seurat_obj$cells_id[names(seurat_obj@active.ident)], type=seurat_obj$type),
        file = "/output/multisample-cells.csv",
        quote = F, col.names = F, row.names = F,
        sep = "\t"
    )
}

vars <- vars[rownames(seurat_obj), c("ENSEMBL", "variance.standardized")]
message("Saving gene and cell data...")
write.table(
    vars,
    file = "/output/r-out-dispersions.csv",
    sep = ",",
    col.names = F, row.names = F
)

write.table(
    colnames(seurat_obj),
    file = "/output/r-out-cells.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

write.table(
    seurat_obj@misc[["gene_annotations"]][seurat_obj@misc[["gene_annotations"]]$input%in%rownames(seurat_obj), ],
    file = "/output/r-out-annotations.csv",
    quote = F, col.names = F, row.names = F,
    sep = "\t"
)

message("saving normalized matrix...")
Matrix::writeMM(t(seurat_obj@assays[[seurat_obj@active.assay]]@data
                  ), file = "/output/r-out-normalized.mtx")

message("saving raw matrix...")
Matrix::writeMM(t(
    seurat_obj@assays[["RNA"]]@counts[
        rownames(seurat_obj),
        colnames(seurat_obj)
    ]),
    file = "/output/r-out-raw.mtx", verb
)



################################################
## SAVING CONFIG FILE 
################################################
# 
# We are going to store the final config to config_dataProcessing.json, in order to upload to dynamoDB.
# The unisample experiments does not require any change, but for the multisample experiment we need
# to add the filtering parameter for each sample (only in the steps that is required.)
# We are going to differentiate in samples only in the steps:
# --> cellSizeDistribution
# --> numGenesVsNumUmis
#
# For both of them, we will run again the step fn for each sample (samples names are stored in metadata type)

# Function to recompute the step fn and store the new config of each sample inside the latest config file 
# We need to iterate per sample and compute separately the step fn.
# Example of structure:
# {
# “filterSettings”: {
#   “probabilityThreshold”: 0.2,
#   “binStep”: 0.05
# },
# “sample-KO”: {
#   “filterSettings”: {
#     “probabilityThreshold”: 0.1,
#     “binStep”: 100
#   }
# },
# “sample-WT1": {
#                     “filterSettings”: {
#                         “probabilityThreshold”: 0.1,
#                         “binStep”: 45
#                     }
#                 }
# }

add_config_per_sample <- function(step_fn, config, scdata, samples){
  
  # We upadte the config file, so to be able to access the raw config we create a copy
  config.raw <- config
  
  for(sample in samples){
    # Downsample the seurat object to a unisample experiment
    scdata_sample <- subset(scdata, type %in% sample)
    # Run the step fun with the unisample experiment and keep the config result
    result_config <- step_fn(scdata_sample, config.raw)$config
    # Inside the config of the samples we are not storing the auto and enable settings, so we remove them
    result_config$auto <- NULL
    result_config$enabled <- NULL
    # Update config with the unisample thresholds
    config[[sample]] <- result_config
  }
  
  return(config)
  
}

# Only recompute in multisample case
if (as.logical(config$samples$multisample)){
  
  config.cellSizeDistribution <- add_config_per_sample(cellSizeDistribution, config.cellSizeDistribution, seurat_obj, unique(seurat_obj$type))
  config.numGenesVsNumUmis <- add_config_per_sample(numGenesVsNumUmis, config.numGenesVsNumUmis, seurat_obj, unique(seurat_obj$type))

}

# When we remove the steps from data-ingest we need to change here the default config. 
# Save config for all steps. 
config <- list(
  cellSizeDistribution = config.cellSizeDistribution
  , mitochondrialContent = config.mitochondrialContent
  , classifier = config.classifier
  , numGenesVsNumUmis = config.numGenesVsNumUmis
  , doubletScores = config.doubletScores
  , dataIntegration = config.dataIntegration
  , computeEmbedding = config.computeEmbedding
)

# Export to json
exportJson <- RJSONIO::toJSON(config, pretty = T)
# The RJSONIO library add '' to boolean keys, so we will remove them.
exportJson <- gsub('\"true\"', "true", exportJson)
exportJson <- gsub('\"false\"', "false", exportJson)
# Trnasform null into []
exportJson <- gsub('null', "[]", exportJson)
message("config file...")
write(exportJson, "/output/config_dataProcessing.json")




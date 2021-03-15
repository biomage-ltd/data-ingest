library(Seurat)
library(Matrix)
require(data.table)
library(gprofiler2)
library(dplyr)

set.seed(123)
options(future.globals.maxSize= 1000 * 1024 ^ 2)
source("src/help.r")
source("src/QC_helpers/cellSizeDistribution.r")
source("src/QC_helpers/mitochondrialContent.r")
source("src/QC_helpers/classifier.r")
source("src/QC_helpers/numGenesVsNumUmis.r")
source("src/QC_helpers/doubletScores.r")
source("src/QC_helpers/dataIntegration.r")
source("src/QC_helpers/computeEmbedding.r")


################################################
## LOADING SPARSE MATRIX AND CONFIGURATION
################################################

message("reloading old matrices...")
scdata <- readRDS("/output/pre-doublet-data.rds")

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")
metadata <- check_config(scdata, config)

message("Creating Seurat Object...")
seurat_obj <- Seurat::CreateSeuratObject(scdata$filtered, assay='RNA', min.cells=3, min.features=200, meta.data=metadata)

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

 # To be consistent we will conver to fraction
  seurat_obj$fraction.mt <- seurat_obj$percent.mt/100

}

message("getting scrublet results...")
# Q to Juanlu: which object is this, scdata or seurat_obj?
# [Bug] seb: it seems scdata, but in previous version it got overwritten (everthing was scdata), how could this ever work?
scores <- get_doublet_score(scdata)

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

  message("Adding emptyDrops scores information...")
  seurat_obj@meta.data <- meta.data
  # dump emptyDrops results into tools slot
  # seurat_obj@tools$CalculateEmptyDrops <- emptydrops_out
} else {
  # TODO: or should this be saved in config?
  message("emptyDrops results not present, skipping...")
  seurat_obj@tools$flag_filtered <- TRUE
}


################################################
## DATA PROCESSING
################################################

#[HARDCODED]

config.cellSizeDistribution <- list(enabled="true", auto="true", 
    filterSettings = list(minCellSize=10800, binStep = 200)
)

config.mitochondrialContent <- list(enabled="true", auto="true", 
    filterSettings = list(method="absolute_threshold", methodSettings = list(
        absolute_threshold=list(maxFraction=0.1, binStep=0.05)
        )
    )
)

config.classifier <- list(enabled="true", auto="true", 
    filterSettings = list(minProbability=0.82, bandwidth=-1, filterThreshold=-1)
)

config.numGenesVsNumUmis <- list(enabled="true", auto="true", 
    filterSettings = list(regressionType = "gam", regressionTypeSettings = list(
        "gam" = list(p.level=0.001)
        )
    )
)

config.doubletScores <- list(enabled="true", auto="true", 
    filterSettings = list(probabilityThreshold = 0.2 , binStep = 0.05)
)

config.dataIntegration <- list(enabled="true", auto="true", 
    dataIntegration = list(method = "seuratv3", methodSettings = list(seuratv3=list(numGenes=2000, normalisation="LogNormalize"))),
    dimensionalityReduction = list(method = "rpca", numPCs = 30, excludeGeneCategories = c())
)

config.computeEmbedding <- list(enabled="true", auto="true", 
    embeddingSettings = list(method = "umap", methodSettings = list(
                                umap = list(minimumDistance=0.2, distanceMetric="euclidean"), 
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

# Waiting filter 1
result.step1 <- cellSizeDistribution(seurat_obj, config.cellSizeDistribution)
# result.step1$config$filterSettings$minCellSize
# str(result.step1$plotData)
# List of 2
#  $ plot1: Named num [1:11217] 3483 6019 3892 3729 4734 ...
#   ..- attr(*, "names")= chr [1:11217] "u" "u" "u" "u" ...
#  $ plot2:List of 2
#   ..$ : Named num [1:11217] 3483 6019 3892 3729 4734 ...
#   .. ..- attr(*, "names")= chr [1:11217] "u" "u" "u" "u" ...
#   ..$ : Named int [1:11217] 6807 1276 1894 4438 16 6867 887 3494 10873 4161 ...
#   .. ..- attr(*, "names")= chr [1:11217] "rank" "rank" "rank" "rank" ...
#
# Step 2: Mitochondrial content filter
#

# Be aware that all the currents experiment does not have the slot fracion.mt, but they have the percent.mt. 
# To be consistent, in this new version I have transformed to fracion.mt [Line 46]
result.step2 <- mitochondrialContent(result.step1$data, config.mitochondrialContent)
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

# Waiting filter 3
result.step3 <- classifier(result.step2$data, config.classifier)
# str(result.step3$plotData)
# result.step3 <- result.step2

#
# Step 4: Number of genes vs number of UMIs filter
#

result.step4 <- numGenesVsNumUmis(result.step3$data, config.numGenesVsNumUmis)

## plotData plots
## Plot 1 (scatter plot with bands)
# plot(result.step4$plotData$plot1$log10_UMIs, result.step4$plotData$plot1$log10_genes, ylab = "log10_genes", xlab = "log10_UMIs")
# lines(result.step4$plotData$plot1$log10_UMIs, result.step4$plotData$plot1$upper_cutoff, col = "red")
# lines(result.step4$plotData$plot1$log10_UMIs, result.step4$plotData$plot1$lower_cutoff, col = "red")


#
# Step 5: Doublet scores filter
#

result.step5 <- doubletScores(result.step4$data, config.doubletScores)
## plotData plots
## Plot 1 (histogram with fraction MT)
# h=hist(result.step5$plotData$plot1,plot=FALSE)
# h$density = h$counts/sum(h$counts)
# plot(h,freq=FALSE)


#
# Step 6: Data integration
#

result.step6 <- dataIntegration(result.step5$data, config.dataIntegration)

#
# Step 7: Compute embedding
#

result.step7 <- computeEmbedding(result.step6$data, config.computeEmbedding)


seurat_obj <- result.step7$data

message("Storing gene annotations...")
seurat_obj@misc[["gene_annotations"]] <- annotations

message("Storing cells id...")
# Keeping old version of ids starting from 0
seurat_obj$cells_id <- 0:(nrow(seurat_obj@meta.data)-1)

message("Storing dispersion...")
# Convert to Gene Symbol
# [Bug] seb: 
# For multi-sample: Error: Unable to find highly variable feature information for method 'vst'
# but this works (for each sample at a time)? do we even git multi-sample at this stage?
# HVFInfo(object = data.split[[i]], assay = "RNA", selection.method = 'vst') # to create vars
vars <- HVFInfo(object = seurat_obj, assay = "RNA", selection.method = 'vst') # to create vars
vars$SYMBOL <- annotations$name[match(rownames(vars), annotations$input)]
vars$ENSEMBL <- rownames(vars)
seurat_obj@misc[["gene_dispersion"]] <- vars

pdf("/output/umap.pdf")
DimPlot(seurat_obj, reduction = "umap")
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

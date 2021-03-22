################################################
## STEP 6. Data integration
#################################################
# Data integration step where batch effect is corrected through data integration methods such as "Seurat V3". 
# The data integration step include the normalization and the PCA analysis.

# NEED TO DISCUSS: functio to compute default parameters

############ NEED TO CHECK CONFIG SCHEMA

#   "dataIntegration": {
#       "dataIntegration": {
#           "method": "seuratv4",
#           "methodSettings": {
#               "seuratv4": {
#                   "numGenes": 2000,
#                   "normalisation": "logNormalise"
#               }
#           }
#       },
#       "dimensionalityReduction": {
#           "method": "rpca",  --> NOT USE!
#           "numPCs": 30,
#           "excludeGeneCategories": []
#       }
#   },

dataIntegration <- function(scdata, config){

    # Check wheter the filter is set to true or false
    if (as.logical(toupper(config$enabled))){
        # So far we only support Seurat V3
        scdata.integrated <- run_dataIntegration(scdata, config)
            # Compute explained variance for the plot2
        eigValues = (scdata.integrated@reductions$pca@stdev)^2  ## EigenValues
        varExplained = eigValues / sum(eigValues)

        # As a short solution, we are going to store an intermediate slot for the numPCs, since this parameter is required when performing
        # the computeEmdedding. The main reason to do not have in the config.computeEmbedding is that this parameter does not change in the computeEmbedding step.
        scdata.integrated@misc[["numPCs"]] <- config$dimensionalityReduction$numPCs
    }

    else
        scdata.integrated <- scdata
    

    # For now config is not updated, since there is not new changes
    # config <- ...

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.integrated,
        config = config,
        plotData = list(
            # Plot1 --> embedding plot by UMAP 
            plot1 = list(x=scdata.integrated@reductions$umap@cell.embeddings[, 1], y=scdata.integrated@reductions$umap@cell.embeddings[, 1]),
            plot2 = list(proportion_variance_explained = varExplained, nPCs = 1:50)
        )
    )

    return(result)

}

# This function covers
#   - Integrate the data using the variable "type" (in case of data integration method is selected) and normalize using LogNormalize method.
#   - Compute PCA analysis
#   - To visualize the results of the batch effect, an UMAP with default setting has been made. 
# Seurat V3 pipeline (see for other methods: https://satijalab.org/seurat/archive/v3.0/integration.html)
run_dataIntegration <- function(scdata, config){
    
    method <- config$dataIntegration$method
    nfeatures <- config$dataIntegration$methodSettings[[method]]$numGenes
    normalisation <- config$dataIntegration$methodSettings[[method]]$normalisation
    
    # Q: Should be the numPCs input for the RunPCA? Since we are looking into the explained variance I suggest to RunPCA with 50 and only
    # use this parameter for the data integration purpose.
    numPCs <- config$dimensionalityReduction$numPCs

    #HARDCODE
    # We require just to get an overview of the data-integration
    umap_min_distance <- 0.3
    umap_distance_metric <- "euclidean"

    # Currently, we only support Seurat V3 pipeline for the multisample integration
   if(method=="seuratv3"){
        data.split <- SplitObject(scdata, split.by = "type")
        for (i in 1:length(data.split)) {
            data.split[[i]] <- NormalizeData(data.split[[i]], normalization.method = normalisation, verbose = F)
            data.split[[i]] <- FindVariableFeatures(data.split[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }
        data.anchors <- FindIntegrationAnchors(object.list = data.split, dims = 1:numPCs, verbose = FALSE)
        scdata <- IntegrateData(anchorset = data.anchors, dims = 1:numPCs)
        DefaultAssay(scdata) <- "integrated"
    }else{
        # Else, we are in unisample experiment and we only need to normalize 
        scdata <- Seurat::NormalizeData(scdata, normalization.method = normalisation, verbose = F)
        scdata <-Seurat::FindVariableFeatures(scdata, selection.method = "vst", nfeatures = nfeatures, verbose = F)
    }

    # Scale in order to compute PCA
    scdata <- Seurat::ScaleData(scdata, verbose = F)

    # HARDCODE numPCs to 50
    scdata <- Seurat::RunPCA(scdata, npcs = 50, features = VariableFeatures(object=scdata), verbose=FALSE)

    # Compute embedding with default setting to get an overview of the performance of the bath correction
    scdata <- RunUMAP(scdata, reduction='pca', dims = 1:numPCs, verbose = F, umap.method = "uwot-learn", min.dist = umap_min_distance, metric = umap_distance_metric)

    return(scdata)
}

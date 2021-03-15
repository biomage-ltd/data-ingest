################################################
## STEP 7. Compute embedding
#################################################
# Compute embedding step where we run dimensional reduction technniques such as t-SNE and  UMAP. Moreover, the cluster analysis is
# done also in this step. 

############ NEED TO CHECK CONFIG SCHEMA

# nPCs ?

#"configureEmbedding": {
#    "embeddingSettings": {
#        "method": "umap",
#        "methodSettings": {
#            "umap": {
#                "minimumDistance": 0.2,
#                "distanceMetric": "euclidean"
#            },
#            "tsne": {
#                "perplexity": 30,
#                "learningRate": 200
#            }
#        }
#    },
#    "clusteringSettings": {
#        "method": "louvain",
#        "methodSettings": {
#            "louvain": {
#                "resolution": 0.5
#            }
#        }
#    }
#},

computeEmbedding <- function(scdata, config){

    # Check wheter the filter is set to true or false
    if (as.logical(toupper(config$enabled)))
        scdata.embedding <- run_computeEmbedding(scdata, config)
    else
        scdata.embedding <- scdata.embedding

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata.embedding,
        config = config,
        plotData = list(
            ###################################
            ####
            #### TODO: return plot information
            ####
            ###################################
        )
    )

    return(result)

}

# This function covers
#   - Compute embedding: t-SNE or UMAP
#   - Clustering
run_computeEmbedding <- function(scdata, config){
    
    #################
    # Embedding part
    #################

    # HARDCODE. Common hardcode, need to be reviewd and try to persist from the previous step. 
    pca_nPCs <- 30


    if (config$embeddingSettings$method=="umap"){
        
        minimumDistance <- config$embeddingSettings$methodSettings$umap$minimumDistance
        distanceMetric <- config$embeddingSettings$methodSettings$umap$distanceMetric

        message("Running embedding --> umap")

        scdata <- RunUMAP(scdata, reduction='pca', dims = 1:pca_nPCs, verbose = F, umap.method = "uwot-learn", min.dist = minimumDistance, metric = distanceMetric)


    }

    if (config$embeddingSettings$method=="tsne"){
        
        perplexity <- config$embeddingSettings$methodSettings$tsne$perplexity
        learningRate <- config$embeddingSettings$methodSettings$tsne$learningRate
        
        if (as.logical(toupper(config$auto))){
            perplexity <- min(30, ncol(scdata)/100)
            learningRate <- max(200, ncol(scdata)/12)
        }


        message("Running embedding --> tsne")
        
        scdata <- RunTSNE(scdata, reduction = 'pca', dims = 1:pca_nPCs, perplexity = perplexity, learning.rate = learningRate)

    }


    #################
    # Clustering part
    #################

    if(config$clusteringSettings$method=="louvain"){
        clustering_method <- 1 #"Louvain"
        clustering_resolution <- config$clusteringSettings$methodSettings[[clustering_method]]$resolution

        # HARDCODE
        annoy.metric = "cosine"
        scdata <- FindNeighbors(scdata, k.param = 20, annoy.metric = annoy.metric, verbose=FALSE)
        scdata <- FindClusters(scdata, resolution=clustering_resolution, verbose = FALSE, algorithm = clustering_method)

    }

    return(scdata)
}

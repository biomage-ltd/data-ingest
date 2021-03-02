################################################
## STEP 5. Doublet score filter 
#################################################
#
# doubletScores function 
#' @description Filters seurat object based on doubletScores scores
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_doubletScores)
#'          - filterSettings: slot with thresholds
#'                  - probabilityThreshold: Float. cut-off for the maximun probability scores that could have a cell. The doubletScores scores have been computed
#'                  through "scrublet" [1]
#'                  - binStep: Float. Bin size for the histogram
#' @export return a list with the filtered seurat object by doublet score, the config and the plot values

doubletScores <- function(scdata, config){

    # Check wheter the filter is set to true or false
    if(!as.logical(toupper(config$enabled)))
        return(scdata)

    # Check if the experiment has MT-content
    if(!"doubletScores_scores"%in%colnames(scdata@meta.data)){
        message("Warning! No doubletScores scores has been computed for this experiment!")
        return(scdata)
    }
    

    # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
    # to get a list with two elements {probabilityThreshold and binStep}.
    if(as.logical(toupper(config$auto)))
        config$filterSettings <- generate_default_values_doubletScores(scdata, config)

    # Information regarding number of UMIs per cells is pre-computed during the 'CreateSeuratObject' function. 
    scdata <- subset(scdata, subset = doubletScores_scores <= config$filterSettings$probabilityThreshold)
    
    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata,
        config = config,
        plotData = list(
            plot1 = c(),
            plot2 = c()
        )
    )

    return(result)

}


# [1] Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. 
# Cell Syst. 2019 Apr 24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. Epub 2019 Apr 3. PMID: 30954476; PMCID: PMC6625319.
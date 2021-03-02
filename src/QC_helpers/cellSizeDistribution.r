################################################
## STEP 1. Cell size distribution filter
#################################################
#
#' @description Filters seurat object based on cell size distribution
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_cellSizeDistribution)
#'          - filterSettings: slot with thresholds
#'                  - minCellSize: Integer. Threshold that contain the minimun number of UMIs per cell
#'                  - binStep: Integer. Bin size for the histogram
#' @export return a list with the filtered seurat object by cell size ditribution, the config and the plot values

cellSizeDistribution <- function(scdata, config){
    
    # Check wheter the filter is set to true or false
    if(!as.logical(toupper(config$enabled)))
        return(scdata)

    # Check if it is required to compute sensible values. From the function 'generate_default_values_cellSizeDistribution', it is expected
    # to get a list with two elements {minCellSize and binStep}
    if(as.logical(toupper(config$auto)))
        config$filterSettings <- generate_default_values_cellSizeDistribution(scdata, config)

    # Information regarding number of UMIs per cells is pre-computed during the 'CreateSeuratObject' function. 
    scdata <- subset(scdata, subset = nCount_RNA >= minCellSize)

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
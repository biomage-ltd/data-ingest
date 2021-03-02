################################################
## STEP 2. Mitochondrial content filter 
#################################################
#
#' @description Filters seurat object based on mitochondrialContent
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_mitochondrialContent)
#'          - filterSettings: slot with thresholds
#'                  - method: String. Method to be used {absolute_threshold} 
#'                  - methodSettings: List with the method as key and contain all the filterSettings for this specific method. 
#'                          * absolute_threshold: based on a cut-off threshold
#'                                  - maxFraction: Float. maximun pct MT-content that we considere for a alive cell
#'                                  - binStep: Float. Bin size for the histogram
#'                          * we are supposed to add more methods ....
#' @export return a list with the filtered seurat object by mitochondrial content, the config and the plot values

mitochondrialContent <- function(scdata, config){

    # Check wheter the filter is set to true or false
    if(!as.logical(toupper(config$enabled)))
        return(scdata)

    # Check if the experiment has MT-content
    if(!"percent.mt"%in%colnames(scdata@meta.data)){
        message("Warning! No MT-content has been computed for this experiment!")
        return(scdata)
    }

    # Check if it is required to compute sensible values. From the function 'generate_default_values_mitochondrialContent', it is expected
    # to get a list with two elements {method and methodSettings}. Depending on the method, the element methodSettings should contain different information. 
    # In the case of 'absolute_threshold', the schema is: 
    # list("method" = "absolute_threshold", "methodSettings" = list("absolute_threshold" = list("maxFraction" = maxFraction_value, "binStep" = binStep_value)))
    if(as.logical(toupper(config$auto)))
        config$filterSettings <- generate_default_values_mitochondrialContent(scdata, config)


    # Check methods
    if(config$filterSettings$method=="absolute_threshold"){
        # Information regarding number of UMIs per cells is pre-computed during the 'CreateSeuratObject' function. 
        scdata <- subset(scdata, subset = percent.mt <= config$filterSettings$method$absolute_threshold$maxFraction)
    }

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
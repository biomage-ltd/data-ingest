################################################
## STEP 3. Classifier filter 
#################################################
#
#' @description Filters seurat object based on classifier filter
#' @param config list containing the following information
#'          - enable: true/false. Refering to apply or not the filter.
#'          - auto: true/false. 'True' indicates that the filter setting need to be changed depending on some sensible value (it requires
#'          to call generate_default_values_classifier)
#'          - filterSettings: slot with thresholds
#'                  - minProbabiliy: 
#'                  - filterThreshold: 
#' @export return a list with the filtered seurat object by probabilities classifier, the config and the plot values
classifier <- function(scdata, config){

    # Check wheter the filter is set to true or false
    if(!as.logical(toupper(config$enabled)))
        return(scdata)
    
    # Check if it is required to compute sensible values. From the function 'generate_default_values_classifier', it is expected
    # to get a list with two elements {minProbabiliy and filterThreshold}.
    if(as.logical(toupper(config$auto)))
        config$filterSettings <- generate_default_values_classifier(scdata, config)

    ################################################
    ## TODO: Implement filtering
    #################################################


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
################################################
## STEP 2. Mitochondrial content filter 
#################################################
# This is a simplest filter that looks at a threshold value for the mitochondrial content.
# To separate cells with low MT-content from the ones that have a high MT-content what makes us think that are dead.  
# This can be a useful first guess. The settings for such a filter can also contain a simple "probabilityThreshold" setting. 

# The most uses values in MT-content are between [0.1, 0.2]. There are not too much literature about how to compute
# a threshold. For now, we will offer two methods:
# --> Absolute threshold: In order to be not too extrictive the threshold is set to 0.1
# --> Quantile 95th: get the value which leaves the 95% of the values of MT-content to the left. 
generate_default_values_mitochondrialContent <- function(scdata, config) {
   if(config$filterSettings$method == "Absolute threshold")
        # HARDCODE
        threshold <- 0.1
   
   if(threshold$filterSettings$method == "Quantile 95th")
        threshold <- quantile(scdata$percent.mt, 0.95)

    if(is.null(threshold))
        stop("Enter a valid method: {Absolute threshold, Quantile 95th}.")

  return(threshold)
}

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


    # Check if the experiment has MT-content
    if(!"percent.mt"%in%colnames(scdata@meta.data)){
        message("Warning! No MT-content has been computed for this experiment!")
        return(scdata)
    }

    probabilityThreshold <- config$filterSettings$probabilityThreshold
    
    # Check if it is required to compute sensible values. From the function 'generate_default_values_doubletScores', it is expected
    # to get a list with only one element --> probabilityThreshold.
    if (as.logical(toupper(config$auto)))
        probabilityThreshold <- generate_default_values_mitochondrialContent(scdata, config)


    # Check wheter the filter is set to true or false
    if (as.logical(toupper(config$enabled)))
        # Information regarding MT-content is pre-computed during the 'data-ingest'. 
        scdata <- subset(scdata, subset = percent.mt <= probabilityThreshold)
    
    # update config
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = scdata, # scdata filter
        config = config,
        plotData = list(
            # plot 1: histgram of MT-content
            # AAACCCAAGCGCCCAT-1 AAACCCAAGGTTCCGC-1 AAACCCACAGAGTTGG-1
            #              0.161              0.198              0.284  ...
            plot1 = scdata$percent.mt,
            # plot 2: There are two alternavitive:
            #           - Scatter plot with UMIs in the x-axis and MT-content in the y-axis
            #           --> code: plot2 = list(u=scdata$nCount_RNA.mt, "MT-content" = scdata$percent.mt)
            #           - Barplot representing in the x-axis the log10(UMIs) and in the y-axis the MT-content. This option is the one 
            #           that is shown in the mockup.
            #           --> code: plot2 = list(log_10_UMIs=log10(scdata$nCount_RNA), MT_content =mscdata$percent.mt)
            plot2 = list()
        )
    )

    return(result)
}

# There are some config parameters that depends on the data it-self. In this file we are going to create the functions
# that allow us to compute the best config parameter for Data Processing in the doubletScores step.

# To identify intelligently the treshold we are going to use the logic inside scDblFinder, which creates a classification
# (singlet our doublet) [ref: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/2_scDblFinder.html#thresholding-and-local-calibration]
# To set the auto value we are going to use as a threshold the maximun score that is given to a singlet. 

doubletScores_config <- function(scdata, config){

    # Minimun score that has a singlet 
    probabilityThreshold <-  max(scdata$doublet_scores[scdata$doublet_class=="singlet"])
    # update config
    config$filterSettings$probabilityThreshold <- probabilityThreshold

    return(config)
}



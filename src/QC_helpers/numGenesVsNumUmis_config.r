# There are some config parameters that depends on the data it-self. In this file we are going to create the functions
# that allow us to compute the best config parameter for Data Processing in the numGenesVsNumUmis step.

numGenesVsNumUmis_config <- function(scdata, config){

    # Sensible values are based on the funciton "gene.vs.molecule.cell.filter" from the pagoda2 package
    p.level <-  min(0.001, 1/ncol(scdata))
    # update config
    config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

    return(config)
}



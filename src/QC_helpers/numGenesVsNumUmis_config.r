numGenesVsNumUmis_config <- function(scdata, config){

    # Sensible values are based on the funciton "gene.vs.molecule.cell.filter" from the pagoda2 package
    p.level <-  min(0.001, 1/ncol(scdata))
    # update config
    config$filterSettings$regressionTypeSettings[[config$filterSettings$regressionType]]$p.level <- p.level

    return(config)
}



# There are some config parameters that depends on the data it-self. In this file we are going to create the functions
# that allow us to compute the best config parameter for Data Processing in the cellSizeDistribution step.


# CalculateBarcodeInflections calculates an adaptive inflection point ("knee")
# of the barcode distribution for each sample group. This is
# useful for determining a threshold for removing low-quality
# samples.
generate_default_values_cellSizeDistribution <- function(seurat_obj, config) {
  # Q: should we check for precalculated values? e.g.:
  # is.null(Tool(seurat_obj, slot = "CalculateBarcodeInflections"))
  # Returns Seurat object with a new list in the `tools` slot,
  # `CalculateBarcodeInflections` including inflection point calculatio

  seurat_obj_tmp <- CalculateBarcodeInflections(
                                        object = seurat_obj,
                                        barcode.column = "nCount_RNA",
                                        group.column = "orig.ident",
                                        # [HARDCODED]
                                        threshold.low = 1e2,
                                        threshold.high = NULL
  )
  # extracting the inflection point value which can serve as minCellSize threshold
  # all other side effects to the scdate object will be discarded
  # [TODO unittest]: this calculation assumes a single sample, i.e. only one group in orig.ident!
  # otherwise it will return multiple values for each group!
  # returned is both the rank(s) as well as inflection point
  # orig.ident          nCount_RNA  rank
  # SeuratProject       1106        10722
  tmp <- Tool(seurat_obj_tmp, slot = "CalculateBarcodeInflections")$inflection_points
  # extracting only inflection point(s)
  return(tmp$nCount_RNA)
}

cellSizeDistribution_config <- function(seurat_obj, config) {
        
    minCellSize <- generate_default_values_cellSizeDistribution(seurat_obj, config)
    # update config
    config$filterSettings$minCellSize <- minCellSize

    return(config)
}

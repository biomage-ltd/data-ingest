################################################
## STEP 1. Cell size distribution filter
#################################################
# cell_size_distribution_filter function 
# This is a simplest filter that looks at the shape of the cell size (# of UMIs per cell) distribution and looks for some local minima, minimum of second derivative, etc. 
# To separate the large cell barcodes that correspond to real cells from the tail containing 'empty droplets'. 
# This can be a useful first guess. The settings for such a filter can also contain a simple "min cell size" setting. 
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

cellSizeDistribution <- function(seurat_obj, config) {
    
    minCellSize <- as.numeric(config$filterSettings$minCellSize)
    # extract plotting data of original data to return to plot slot later
    plot1_data  <- seurat_obj$nCount_RNA
    # plot 1: histgram of UMIs, hence input is all UMI values, e.g.
    # AAACCCAAGCGCCCAT-1 AAACCCAAGGTTCCGC-1 AAACCCACAGAGTTGG-1
    #               2204              20090               5884  ...
    names(plot1_data) <- rep("u", length(plot1_data))
    #    u    u    u    u    u    u
    # 3483 6019 3892 3729 4734 3244
    plot2_data <- order(seurat_obj$nCount_RNA)
    names(plot2_data) <- rep("rank", length(plot2_data))

    # Check if it is required to compute sensible values. From the function 'generate_default_values_cellSizeDistribution', it is expected
    # to get a list with two elements {minCellSize and binStep}
    if (as.logical(toupper(config$auto)))
        # config not really needed for this one (maybe later for threshold.low/high):
        minCellSize <- generate_default_values_cellSizeDistribution(seurat_obj, config)

    # Check wheter the filter is set to true or false
    # if true overwrite seurat object (side effect!)
    if (as.logical(toupper(config$enabled))) {
        # Information regarding number of UMIs per cells is pre-computed during the 'CreateSeuratObject' function. 
        seurat_obj <- subset(seurat_obj, subset = nCount_RNA >= minCellSize)
    }

    # update config
    config$filterSettings$minCellSize <- minCellSize
    

    # the result object will have to conform to this format: {data, config, plotData : {plot1, plot2}}
    result <- list(
        data = seurat_obj,
        config = config,
        plotData = list(
            cellSizeDistributionHistogram = plot1_data,
            # Q: are both plots updated for this filter?
            # Q: what is the format of plot2?
            # knee-plot: this is on a log-log scale, are logs calucated here or on the UI?
            # cells are ordered on the x-axis according to the number of distinct UMIs observed. 
            # The y-axis displays the number of distinct UMIs for each barcode (here barcodes are proxies for cells).
            # cellRank_sorted.json: [{"u": 0, "rank": 17852}, {"u": 1, "rank": 17412},...]
            # plot2 = list(u = seurat_obj$nCount_RNA, rank = order(seurat_obj$nCount_RNA))
            cellSizeDistributionKneePlot = list(plot1_data, plot2_data)
        )
    )

    return(result)
}

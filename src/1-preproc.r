library(Seurat)
library(magrittr)
library(parallel)
library(Matrix)
require(data.table)
library(RJSONIO)
library(ggplot2)
library(MASS)
library(DropletUtils)


# Read10X_data function from Seurat V3 [1]
#' @description Flexible function to load 10x data. It allows different format:
#'      - matrix.mtx, barcodes.tsv, genes.tsv  
#'      - matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz  
#'      - matrix.mtx, barcodes.tsv, features.tsv  
#'      - ... any possible combination that have at least  matrix, barcodes and genes/features data.
#'      It also works with multiplesamples
#' @param data.dir file to look for
#' @param gene.column gene column number
#' @param unique.features not duplicated features
#' @param strip.suffix boolean to remove -1 at the end of barcodes (default FALSE)
#' @param sample in case of multiplesample, it refers to the name of the samples
#' 
#' @return dgCMatrix with genes as rows and cells as column. In the case of multisample,
#'  the cell names have a prefix _ corresponding to the type metadata from the json.

Read10X_data <- function (data.dir = NULL, gene.column = 2, unique.features = TRUE, 
                          strip.suffix = FALSE, samples=NULL) {
  full.data <- list()
  # Create a matrix per sample
  for (i in 1:(ifelse(is.null(samples), 1, length(samples)))) {
    run <- data.dir
    
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    
    # index of json files to not take into account
    idt_json = grep("json", list.files(run))

    # Checking the three files for the sample i
    files_10x <- ifelse(rep(is.null(samples[i]), 3), list.files(run)[-idt_json], list.files(run, recursive = T, pattern = samples[i]))
    # Add the name or leave it empty if there is only one sample
    sample_name = ifelse(is.null(samples[i]), "", paste(samples[i], "_", sep = ""))
    
    if(any(is.na(files_10x)))
      stop("Multisample loading was selected but there is no files with sample ", samples[i], 
           ". Please, be sure that the directory has three files with the sample name in the filename.")
    
    # Lookin for the trhee files: barcode, matrix and features/genes
    barcode.idt <- grep("barcodes.tsv", files_10x)
    matrix.idt <- grep("matrix.mtx", files_10x)
    features.idt <- grep("features.tsv", files_10x)
    
    if(length(features.idt)==0) features.idt <- grep("genes.tsv", files_10x)
    
    if(length(barcode.idt)==0){
      stop("Barcode file missing. Expecting ", 
           basename(path = run), "/", sample_name, "barcodes.tsv")
    }
    
    if(length(matrix.idt)==0){
      stop("Expression matrix file missing. Expecting ", 
           basename(path = run), "/", sample_name, "matrix.mtx")
    }
    
    if(length(features.idt)==0){
      stop("Gene name or features file missing. Expecting ", 
           basename(path = run), "/", sample_name, "genes.tsv or features.tsv")
    }
    
    # Get the full path
    barcode.loc <- file.path(run, files_10x[barcode.idt])
    matrix.loc <- file.path(run, files_10x[matrix.idt])
    features.loc <- file.path(run, files_10x[features.idt])
    
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    
    # Adding the prefix
    colnames(x = data) <- paste0(sample_name, "", cell.names)
    feature.names <- read.delim(file = features.loc, header = FALSE, 
                                stringsAsFactors = FALSE)
    
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}

# create_dataframe function
#' @description read matrix based on config. Possibles input
#'      - 10x scdata
#'      - table
#' @param config experiment settings.
#'
#' @return list with an element named raw with dgCMatrix with genes as rows and cells as column

create_dataframe <- function(config){
  data_type <- config$input["type"]
  scdata <- list()
  
  #Check config format
  if(!"samples"%in%names(config))
    stop("The format of the config is wrong. There should be a category for the sample.
     Current content of config:", names(config))

  if(any(!c("samples_info", "multisample")%in%names(config$samples)))
    stop("The format of the config is wrong. There should be a category inside sample with the multisample information.
     Current content of config:", names(config$samples))

  if(!"type"%in%names(config$samples$samples_info))
    stop("The format of the config is wrong. There should be a category for the samples_info inside sample specifing the sample names in the case of 
    multisample or empty otherwise.") 

  if (data_type == "10x"){
    message("Loading 10x data set from input folder.")
    samples <- NULL
    if(as.logical(config$samples[["multisample"]]))
      samples <- config$samples$samples_info$type
    scdata$raw <- Read10X_data("/input", unique.features=TRUE, gene.column = 1, samples=samples)
  }
  
  if (data_type == "table") {
    path <- config$input["path"]
    message(paste("Loading table-type data set from", path))
    scdata$raw <- as.matrix(read.table(paste("/input", path, sep = "")))
  }
  
  message(
    paste(
      "Found", nrow(scdata$raw), "genes and", ncol(scdata$raw), "cells."
    )
  )
  
  return(scdata)
}

# prepare_scrublet_table function 
#' @description Save raw values before doublet filtering
#' @param scdata list with an elment named `filtered` with raw values before doublet filtering
#' 
#' @export save matrix as pre-doublet-matrix.csv in output directory

prepare_scrublet_table <- function(scdata) {
  table <- data.table(
    as.matrix(
      t(
        scdata$filtered
      )
    )
  , keep.rownames=T)
  
  path <- "/output/pre-doublet-matrix.csv"
  
  file.create(path)
  data.table::fwrite(table, file = path, row.names = F)
}

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")

message("Creating raw dataframe...")
# Note, scdata is not a seurat-object:
# List of 1
#  $ raw     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
scdata <- create_dataframe(config)


message("calculating probability of barcodes being background noise...")

# [HARDCODE]
# threshold_emptydrops <- 100
threshold_emptydrops <- 2000 # use this to simulate unfiltered data
# automatically checking if data contains enough (50) barcodes that 
# are confidently considered to be empty to serve as training set
# not necessary as of now: sce <- as.SingleCellExperiment(scdata)
nr_barcodes_ambient <- sum(colSums(scdata$raw) < threshold_emptydrops) 
range(colSums(scdata$raw) ) 
if (nr_barcodes_ambient < 50 ) {
  flag_filtered <- TRUE
  message("Found not enough [", nr_barcodes_ambient, "] ambient barcodes,
          so data is considered to be pre-filtered")
  message("Smallest nr of UMI for any barcode: ", min(colSums(scdata$raw)))
  message("EmptyDrops filter will not be applied")
} else {
  flag_filtered <- FALSE
  # if non-filtered, run emptyDrops and save results in file
  message("Found enough [", nr_barcodes_ambient, "] ambient barcodes, 
          so data is considered to be not filtered")
  emptydrops_out <- emptyDrops(scdata$raw, lower = threshold_emptydrops) 
  # TODO: save flag_filtered in DynamoDB or config (now: checking if file below does exists)
  saveRDS(emptydrops_out, file = "/output/pre-emptydrops-data.rds", compress = FALSE)
  # scdata@tools$CalculateEmptyDrops <- emptydrops_out
}

message("Filtering cells by size")
# Default parameter by Seurat "min.features" of CreateSeuratObject fn (url: https://satijalab.org/seurat/archive/v3.2/pbmc3k_tutorial.html)
# [HARDCODE]
min.features <- 200
scdata$filtered <- scdata$raw[, Matrix::colSums(scdata$raw>0)>= min.features]

message("Filtering cells by molecules/gene...")
# Default parameter by Seurat "min.cells"  of CreateSeuratObject fn (url: https://satijalab.org/seurat/archive/v3.2/pbmc3k_tutorial.html)
# [HARDCODE]
min.cells <- 3
scdata$filtered <- scdata$filtered[Matrix::rowSums(scdata$filtered>0)> min.cells,]

message("Exporting pre-scrublet scdata...")
prepare_scrublet_table(scdata)
saveRDS(scdata, file = "/output/pre-doublet-data.rds", compress = FALSE)

####### CITATIONS #######
##
## [1] Stuart and Butler et al. Comprehensive Integration of Single-Cell Data. Cell (2019) [Seurat V3]
##

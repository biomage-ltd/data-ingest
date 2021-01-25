library(Seurat)
library(magrittr)
library(parallel)
library(Matrix)
require(data.table)
library(RJSONIO)
library(ggplot2)
library(MASS)


# Read10X_data function
#' @description Flexible function to load 10x data. It allows different format:
#'      - matrix.mtx, barcodes.tsv, genes.tsv  
#'      - matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz  
#'      - matrix.mtx, barcodes.tsv, features.tsv  
#'      - ... any possible combination that have at least  matrix, barcodes and genes/features data.
#' @param data.dir file to look for
#' @param gene.column gene column number
#' @param unique.features not duplicated features
#' @param strip.suffix boolean to remove -1 at the end of barcodes (default FALSE)
#' 
#' @return dgCMatrix with genes as rows and cells as column

Read10X_data <- function (data.dir = NULL, gene.column = 2, unique.features = TRUE, 
          strip.suffix = FALSE) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    
    files_10x <- list.files(run)
    
    barcode.idt <- grep("barcodes.tsv", files_10x)
    matrix.idt <- grep("matrix.mtx", files_10x)
    features.idt <- grep("features.tsv", files_10x)
    
    if(length(features.idt)==0) features.idt <- grep("genes.tsv", files_10x)
    
    
    if(length(barcode.idt)==0){
      stop("Barcode file missing. Expecting ", 
           basename(path = run), "/barcodes.tsv")
    }
    
    if(length(matrix.idt)==0){
      stop("Expression matrix file missing. Expecting ", 
           basename(path = run), "/matrix.mtx")
    }
    
    if(length(features.idt)==0){
      stop("Gene name or features file missing. Expecting ", 
           basename(path = run), "/genes.tsv or features.tsv")
    }
    
    barcode.loc <- file.path(run, files_10x[barcode.idt])
    matrix.loc <- file.path(run, files_10x[matrix.idt])
    features.loc <- file.path(run, files_10x[features.idt])
    
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
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
#'      - 10x data
#'      - table
#' @param config experiment settings.
#'
#' @return list with an element named raw with dgCMatrix with genes as rows and cells as column

create_dataframe <- function(config){
    data_type <- config$input["type"]
    data <- list()

    if (data_type == "10x"){
        message("Loading 10x data set from input folder.")

        data$raw <- Read10X_data("/input", unique.features=TRUE)
    }

    if (data_type == "table") {
        path <- config$input["path"]

        message(paste("Loading table-type data set from", path))

        data$raw <- as.matrix(read.table(paste("/input", path, sep = "")))
    }

    message(
        paste(
            "Found", nrow(data$raw), "genes and", ncol(data$raw), "cells."
        )
    )
  
    return(data)
}

# prepare_scrublet_table function 
#' @description Save raw values before doublet filtering
#' @param data list with an elment named `filtered` with raw values before doublet filtering
#' 
#' @export save matrix as pre-doublet-matrix.csv in output directory

prepare_scrublet_table <- function(data) {
    table <- data.table(
        as.matrix(
            t(
                data$filtered
            )
        )
    )

    path <- "/output/pre-doublet-matrix.csv"

    file.create(path)
    data.table::fwrite(table, file = path)
}

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")

message("Creating raw dataframe...")
data <- create_dataframe(config)

message("Filtering cells by size")
data$filtered <- data$raw[, Matrix::colSums(data$raw)>=1e3]

message("Filtering cells by molecules/gene...")
data$filtered <- data$filtered[Matrix::rowSums(data$filtered>0)>5,]

message("Exporting pre-scrublet data...")
prepare_scrublet_table(data)
saveRDS(data, file = "/output/pre-doublet-data.rds", compress = FALSE)

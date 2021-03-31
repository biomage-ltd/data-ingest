################################################
## 1_Preproc.r
##  - Read input folder {10x data}
##  - Prepare rds with a list of raw counts matrix per sample to compute emptyDrops 
##  - Prepare a csv per sample to run scrublet
################################################



suppressWarnings(library(Seurat))
suppressWarnings(library(magrittr))
suppressWarnings(library(parallel))
suppressWarnings(library(Matrix))
suppressWarnings(require(data.table))
suppressWarnings(library(RJSONIO))
suppressWarnings(library(MASS))


# checking_10x_structure
#' @description The design of the input data needs to be in a particular way. With this function we are going to check if
#' the input folder is designed correctly:
#' intput/
#' ------ sample_name/
#' ------------------- features.tsv.gz
#' ------------------- barcodes.tsv.gz
#' ------------------- matrix.mtx.gz
#' @param samples samples to check
#'
#' @return TRUE if the design is correct FALSE otherwise
check_10x_input <- function(samples){
  
  cell_ranger_v2 <- c("genes.tsv", "barcodes.tsv", "matrix.mtx")
  cell_ranger_v3 <- c("features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz")

  for(sample in samples){
    check_v2 <- all(cell_ranger_v2%in%list.files(paste("/input", sample, sep="/"), full.names = F))
    check_v3 <- all(cell_ranger_v3%in%list.files(paste("/input", sample, sep="/"), full.names = F))

    if(!check_v2 & !check_v3){
      return(FALSE)
    }
  }

  return(TRUE)

}


# create_dataframe function
#' @description read matrix based on config. Possibles input
#'      - 10x scdata
#'      - table
#' @param config experiment settings.
#'
#' @return list with an element per sample with dgCMatrix with genes as rows and cells as column

create_dataframe <- function(config){
  data_type <- config$input["type"]
  scdata <- list()

  # Check config format. It requires to have a 'samples' key. For the multisample experiment, it needs to be the same
  # as the name of the folders that are inside the input folder. In the case of unisample it should be [].
  if(!"samples"%in%names(config))
    stop("The format of the config is wrong. There should be a category for the sample.
     Current content of config:", names(config))
  
  # If no samples have been added, we will take all the samples of the project. 
  if(length(config$samples)==0){
    warning("All samples in input folder will be included in the analysis!")
    samples <- list.dirs("/input", full.names = F)[-1]
  }else{
    samples <- config$samples
    if (!all(samples%in%list.dirs("/input", full.names = F)))
      stop("Check samples to be used in the analysis, 
      since there are some of them that hasn't got a folder with the files: ",
           samples[!samples%in%list.dirs("/input", full.names = F)])
  }
  
  message("Samples to include in the analysis: ", paste(samples, collapse = " - "))
  
  if (data_type == "10x"){
    message("Loading 10x data set from input folder.")

    # Checking design
    if(!check_10x_input(samples)){
      stop("Please! Check in files inside the input folder.", 
      "There should be the files {genes.tsv, matrix.mtx, barcodes.tsv} or {features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz}")
    }

    for(sample in samples){
      scdata[[sample]] <- Read10X(paste("/input", sample, sep = "/"))  
      message(
        paste(
          "Found", nrow(scdata[[sample]]), "genes and", ncol(scdata[[sample]]), "cells in sample", sample, "."
        )
      )
    }
  }
  
  # So far, we haven't tested input format in table type. 
  if (data_type == "table") {
    message(paste("Loading table-type data set from", path))
    for(sample in samples){
      scdata[[sample]] <- as.matrix(read.table(paste("/input", sample, sep = "/")))
      message(
        paste(
          "Found", nrow(scdata$raw[[sample]]), "genes and", ncol(scdata$raw[[sample]]), "cells in sample", sample, "."
        )
      )
    }
  }
  
  return(scdata)
}

# prepare_scrublet_table function 
#' @description Since scrublet cannot run with the original raw matrix (since there are a large amount of cells with empty reads), we need to prefiltered with a minimun
#' threshold. For that, we have used the one propose in Seurat tutorial (https://satijalab.org/seurat/archive/v3.2/pbmc3k_tutorial.html), which is related with cell size 
#' and gene size.
#' 
#' @param scdata Sparse matrix with the counts for one sample.
#' @param sample_name Name of the sample that we are preparing.
#' 
#' @export

prepare_scrublet_table <- function(scdata, sample_name) {

  message(
    "Saving ", sample_name, "..."
  )
  # Default parameter by Seurat "min.cells" and "min.features" of CreateSeuratObject fn (url: https://satijalab.org/seurat/archive/v3.2/pbmc3k_tutorial.html)
  # [HARDCODE]
  min.cells <- 3
  min.features <- 200
  scdata.filtered <- scdata[Matrix::rowSums(scdata>0)>min.cells, Matrix::colSums(scdata>0)>=min.features]
  table <- data.table(
    as.matrix(
      t(
        scdata.filtered
      )
    )
    , keep.rownames=T)
  
  path <- paste("/output/pre-doublet-matrix-", sample_name, ".csv", sep="")
  
  file.create(path)
  data.table::fwrite(table, file = path, row.names = F)
}

message("Loading configuration...")
config <- RJSONIO::fromJSON("/input/meta.json")

# We include in variable scdata_list all the sparse matrix per sample
message("Creating raw dataframe...")
scdata_list <- create_dataframe(config)

# We store the pre-filtered scdata for scrublets per sample
message("Exporting pre-filtered scdata for scrublets...")
sapply(names(scdata_list), function(sample_name) prepare_scrublet_table(scdata_list[[sample_name]], sample_name))

# We store the raw scdata_list since for the emptyDrops since to compute the background we cannot remove any cells. 
message("Exporting raw scdata for emptyDrops...")
saveRDS(scdata_list, file = "/output/pre-doublet-scdata_list.rds", compress = FALSE)

message("Step 1 completed.")
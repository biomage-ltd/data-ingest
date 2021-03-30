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

  #Check config format
  if(!"samples"%in%names(config))
    stop("The format of the config is wrong. There should be a category for the sample.
     Current content of config:", names(config))
  
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
  
  message("Samples to include in the analysis: ", samples)
  
  if (data_type == "10x"){
    message("Loading 10x data set from input folder.")
    for(sample in samples){
      scdata[[sample]] <- Read10X(paste("/input", sample, sep = "/"))  
      message(
        paste(
          "Found", nrow(scdata[[sample]]), "genes and", ncol(scdata[[sample]]), "cells in sample", sample, "."
        )
      )
    }
  }
  
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
#' @description Save raw values before doublet filtering
#' @param scdata list with an elment per sample with raw values before doublet filtering
#' 
#' @export save matrix as pre-doublet-matrix-sample.csv in output directory

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

message("Creating raw dataframe...")
scdata_list <- create_dataframe(config)

message("Exporting pre-filtered scdata for scrublets...")
sapply(names(scdata_list), function(sample_name) prepare_scrublet_table(scdata_list[[sample_name]], sample_name))

message("Exporting raw scdata for emptyDrops...")
saveRDS(scdata_list, file = "/output/pre-doublet-scdata_list.rds", compress = FALSE)

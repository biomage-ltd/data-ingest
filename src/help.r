
# get_doublet_score function 
#' @description Get the cells with its doublet scores computed previously through scrublet
#' @param scdata matrix with barcodes as columns
#' 
#' @export save barcodes and double scores

get_doublet_score <- function(sample) {
    scores <-
        data.table::fread(
            paste("/output/doublet-scores-", sample, ".csv", sep = ""),
        )

    colnames(scores) <- c("barcodes", "doublet_scores")
    rownames(scores) <- scores$barcodes    
    return(as.data.frame(scores))
}



# check_config function 
#' @description Create metadata dataframe from config files
#' @param scdata matrix with barcodes as columns to assign metadata information
#' @param config config list from meta.json
#' 
#' @export save barcodes to keep

check_config <- function(scdata, sample, config){
    metadata <- NULL
    metadata <- data.frame(row.names = colnames(scdata), samples=rep(sample, ncol(scdata)))

    # Check if "type" exists on config file inside samples_info. If it is TRUE, 

    if("metadata" %in% names(config)){
        rest_metadata <- as.data.frame(config$metadata)
        rest_metadata$sample <- ifelse(length(config$samples)>1, config$samples, sample)
        for(var in rest_metadata){
            metadata[, var] <- rest_metadata[, var][match(metadata$sample, rest_metadata$sample)]
        }
    }
    
    return(metadata)
}
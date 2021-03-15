
# get_doublet_score function 
#' @description Get the cells with its doublet scores computed previously through scrublet
#' @param scdata matrix with barcodes as columns
#' 
#' @export save barcodes and double scores

get_doublet_score <- function(scdata) {
    scores <-
        data.table::fread(
            "/output/doublet-scores.csv",
            col.names = c("score")
        )

    scores <- as.data.frame(scores[, "barcodes" := colnames(scdata)])
    rownames(scores) <- scores$barcodes    
    return(scores)
}

# check_config function 
#' @description Create metadata dataframe from config files
#' @param scdata matrix with barcodes as columns to assign metadata information
#' @param config config list from meta.json
#' 
#' @export save barcodes to keep

check_config <- function(scdata, config){
    metadata <- NULL
    
    # Check if "type" exists on config file inside samples_info. If it is TRUE, 
    # we are in multisample experiments and we create metadata with samples names
    # and other attributes that are inside samples_info 
    if("type" %in% names(config$samples$samples_info)){
        metadata <- data.frame(row.names = colnames(scdata$filtered))
        metadata[colnames(scdata$filtered), "type"] <- unlist(lapply(strsplit(colnames(scdata$filtered), "_"), `[`, 1))
        
        rest_metadata <- as.data.frame(config$samples$samples_info)
        for(var in colnames(rest_metadata)[-which(colnames(rest_metadata)%in%"type")]){
            metadata[, var] <- rest_metadata[, var][match(metadata$type, rest_metadata$type)]
        }
    }
    
    return(metadata)
}
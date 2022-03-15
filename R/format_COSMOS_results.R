#' format_COSMOS_res
#'
#' formats the network with readable names
#'
#' @param cosmos_res  results of COSMOS run
#' @param metab_mapping a named vector with HMDB Ids as names and desired metabolite names as values.
#' @return list with network and attribute tables.
#' @importFrom stringr str_extract
#' @export
format_COSMOS_res <- function(cosmos_res,
                              metab_mapping = NULL)
{
  
  if(is.null(metab_mapping))
  {
    data("HMDB_mapper_vec",package = "cosmosR",envir = environment())
  }
  SIF <- as.data.frame(cosmos_res$weightedSIF)
  ATT <- as.data.frame(cosmos_res$nodesAttributes)
  colnames(ATT)[1] <- "Nodes"
  ATT$measured <- ifelse(ATT$NodeType %in% c("M","T","S","P"),1,0)
  ATT$Activity <- ATT$AvgAct
  
  for(i in c(1,3))
  {
    SIF[,i] <- sapply(SIF[,i], function(x, HMDB_mapper_vec){
      x <- gsub("Metab__","",x)
      x <- gsub("^Gene","Enzyme",x)
      suffixe <- stringr::str_extract(x,"_[a-z]$")
      x <- gsub("_[a-z]$","",x)
      if(x %in% names(HMDB_mapper_vec))
      {
        x <- HMDB_mapper_vec[x]
        x <- paste("Metab__",paste(x,suffixe,sep = ""),sep = "")
      }
      return(x)
    },HMDB_mapper_vec = HMDB_mapper_vec)
  }
  
  ATT[,1] <- sapply(ATT[,1], function(x, HMDB_mapper_vec){
    x <- gsub("Metab__","",x)
    x <- gsub("^Gene","Enzyme",x)
    suffixe <- stringr::str_extract(x,"_[a-z]$")
    x <- gsub("_[a-z]$","",x)
    if(x %in% names(HMDB_mapper_vec))
    {
      x <- HMDB_mapper_vec[x]
      x <- paste("Metab__",x,sep = "")
    }
    if(!is.na(suffixe))
    {
      x <- paste(x,suffixe,sep = "")
    }
    return(x)
  },HMDB_mapper_vec = HMDB_mapper_vec)
  
  return(list(SIF,ATT))
}
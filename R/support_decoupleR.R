#' Format Ligand-Receptor Resource
#'
#' This function formats a ligand-receptor resource by creating a gene set
#' with source-target pairs, converting it to a long format, and adding
#' default values for 'mor' and 'likelihood'.
#'
#' @param ligrec_ressource A data frame representing the ligand-receptor resource with columns for source and target gene symbols.
#'
#' @return A data frame containing the formatted ligand-receptor gene set with columns:
#'   \item{gene}{The gene symbol from the ligand-receptor pairs.}
#'   \item{set}{The set identifier combining source and target gene symbols.}
#'   \item{mor}{Default value set to 1 for all entries.}
#'   \item{likelihood}{Default value set to 1 for all entries.}
#'
#' @examples
#' # Create a sample ligand-receptor resource
#' ligrec_ressource <- data.frame(source_genesymbol = c("L1", "L2"),
#'                                target_genesymbol = c("R1", "R2"))
#'
#' # Format the ligand-receptor resource
#' formatted_geneset <- format_LR_ressource(ligrec_ressource)
#'
#' @export
format_LR_ressource <- function(ligrec_ressource)
{
  ligrec_geneset <- ligrec_ressource[,c("source_genesymbol","target_genesymbol")]
  ligrec_geneset$set <- paste(ligrec_geneset$source_genesymbol, ligrec_geneset$target_genesymbol, sep = "___")
  ligrec_geneset <- reshape2::melt(ligrec_geneset, id.vars = "set")[,c(3,1)]
  names(ligrec_geneset)[1] <- "gene"
  ligrec_geneset$mor <- 1
  ligrec_geneset$likelihood <- 1
  ligrec_geneset <- distinct(ligrec_geneset)
  
  return(ligrec_geneset)
}

#' Convert ULM Results to Wide Format
#'
#' This function converts the results from a ULM analysis to a wide format data frame.
#' The input is a data frame with columns for source, condition, and score. The output 
#' is a data frame where each row represents a source and each column represents a condition,
#' with the corresponding scores as values.
#'
#' @param ulm_result A data frame representing the ULM results with columns: source, condition, and score.
#'
#' @return A data frame in wide format where each row is a source and each column is a condition.
#'
#' @examples
#' # Create a sample ULM result
#' ulm_result <- data.frame(source = c("A", "A", "B", "B"),
#'                          condition = c("cond1", "cond2", "cond1", "cond2"),
#'                          score = c(0.5, 0.8, 0.3, 0.7))
#'
#' # Convert to wide format
#' wide_ulm_result <- wide_ulm_res(ulm_result)
#'
#' @export
wide_ulm_res <- function(ulm_result)
{
  ulm_result_df <- reshape2::dcast(ulm_result, formula = source~condition, value.var = "score")
  row.names(ulm_result_df) <- ulm_result_df$source
  ulm_result_df <- ulm_result_df[,-1]
  
  return(ulm_result_df)
}

#' Translate Column Using HMDB Mapper
#'
#' This function translates the values in a column using a provided Human Metabolome Database (HMDB) mapper vector.
#' It modifies the input values by replacing certain prefixes and suffixes according to specific rules.
#'
#' @param my_column A vector of values to be translated.
#' @param HMDB_mapper_vec A named vector where the names are the original identifiers and the values are the corresponding HMDB identifiers.
#'
#' @return A vector with the translated values.
#'
#' @examples
#' # Create a sample column and HMDB mapper vector
#' my_column <- c("Metab__1234_a", "Gene5678_b", "Metab__91011_c")
#' HMDB_mapper_vec <- c("1234" = "HMDB00001", "5678" = "HMDB00002", "91011" = "HMDB00003")
#'
#' # Translate the column
#' translated_column <- translate_column_HMDB(my_column, HMDB_mapper_vec)
#'
#' @export
translate_column_HMDB <- function(my_column, HMDB_mapper_vec)
{
  return(sapply(my_column, function(x, HMDB_mapper_vec) {
    x <- gsub("Metab__", "", x)
    x <- gsub("^Gene", "Enzyme", x)
    suffixe <- stringr::str_extract(x, "_[a-z]$")
    x <- gsub("_[a-z]$", "", x)
    if (x %in% names(HMDB_mapper_vec)) {
      x <- HMDB_mapper_vec[x]
      x <- paste("Metab__", x, sep = "")
    }
    if (!is.na(suffixe)) {
      x <- paste(x, suffixe, sep = "")
    }
    return(x)
  }, HMDB_mapper_vec = HMDB_mapper_vec))
}
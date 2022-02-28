#' add metabolic compartment and metab__ prefix to metabolite IDs
#' 
#' This function adds metabolic compartments to the metabolic identifiers provided by the user.
#' 
#' @param metab_input a named  vector with matebolic statistics as inputs and metabolite identifiers as names
#' @param compartment_codes a character vector, the desired compartment codes to be added. Possible values are "r", "c", "e", "x", "m", "l", "n" and "g"
#' @export
#' @return a named vector with the compartment code and prefixed added to the names
prepare_metab_inputs <- function(metab_input, compartment_codes)
{
  comps <- c("r", "c", "e", "x", "m", "l", "n", "g")
  
  if(length(compartment_codes[which(!(compartment_codes %in% comps))]) > 0)
  {
    ignored <- compartment_codes[which(!(compartment_codes %in% comps))]
    print("The following compartment codes are not found in the PKN and will be ignored")
    print(ignored)
  } 
  
  compartment_codes <- compartment_codes[which(compartment_codes %in% comps)]
  
  if(length(compartment_codes) == 0)
  {
    print("There are no valid compartment left. No compartiment codes will be added.")
    names(metab_input) <- paste("Metab__",names(metab_input),sep = "")
    return(metab_input)
  } else
  {
    print("Adding compartment codes.")
    metab_input_list <- lapply(compartment_codes, function(compartment_code, curr_metab_input){
      names(curr_metab_input) <- paste(names(curr_metab_input),compartment_code, sep = "_")
      names(curr_metab_input) <- paste("Metab__",names(curr_metab_input),sep = "")
      return(curr_metab_input)
    }, curr_metab_input = metab_input)
    metab_input <- unlist(metab_input_list)
    return(metab_input)
  }
}


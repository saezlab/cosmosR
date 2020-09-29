
#' limit_string_vec
#' 
#' writes the first 6 elements and the number of other elements in a string.
#' 
limit_string_vec <- function(str_vec){
    paste(paste(head(str_vec),collapse = ", "), "and",
    length(str_vec) - length(head(str_vec)), "more.")
}

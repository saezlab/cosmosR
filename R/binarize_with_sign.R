#' binarize_with_sign
#'
#' binarizes the data based on the threshold using the absolute value of the data.
#' Gives a sign based on the sign of the input data.
#' @param data numerical vector
#' @param threshold threshold used to binarise the data
binarize_with_sign <- function(data,threshold){
    
    out = sign(data) * as.numeric(abs(data) > threshold)
    return(out)
}

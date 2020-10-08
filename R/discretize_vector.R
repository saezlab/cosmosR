#' discretize input vector
#' 
#' convert a vector of continuous values to discrete {-1, 0, 1}
#' 
#' @param input_vec numerical value


discretize_input <- function(input_vec){
    
    
    if(!all(input_vec %in% c(-1,0,1))){
        message("Input nodes should have values from {-1, 0, 1}. We discretize your input with sign().")
        res = sign(input_vec)
    }else{
        res = input_vec
    }
    
   return(res)
}
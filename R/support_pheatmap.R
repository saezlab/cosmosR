#' Create Heatmap Color Palette
#'
#' This function generates a color palette suitable for heatmaps based on the values in a matrix. It uses the `createLinearColors` function to generate separate color gradients for positive and negative values.
#'
#' @param my_matrix A numeric matrix for which the heatmap color palette is to be generated.
#'
#' @return A character vector of colors representing the heatmap color palette based on the input matrix values.
#'
#' @examples
#' # Create a sample matrix
#' my_matrix <- matrix(c(-3, -1, 0, 1, 3), nrow = 1)
#'
#' # Generate heatmap color palette
#' heatmap_palette <- make_heatmap_color_palette(my_matrix)
#'
#' @export
make_heatmap_color_palette <- function(my_matrix)
{
  t <- as.vector(t(my_matrix))
  palette1 <- createLinearColors(t[t < 0],withZero = F , maximum = abs(min(t,na.rm = T)) * 10)
  palette2 <- createLinearColors(t[t > 0],withZero = F , maximum = abs(max(t,na.rm = T)) * 10)
  palette <- c(palette1, palette2)
}

#' Create Linear Colors Based on Numeric Input
#'
#' This function generates a gradient of colors based on the provided numeric values. The colors can be adjusted to include zero and are configurable with a specified maximum and custom color palette.
#'
#' @param numbers A numeric vector for which the color gradient is to be generated.
#' @param withZero A logical value indicating whether zero should be included in the color gradient. Default is TRUE.
#' @param maximum An integer specifying the maximum number of colors to be generated in the gradient. Default is 100.
#' @param my_colors A character vector of length three specifying the colors to be used in the gradient. Default is c("royalblue3", "white", "red").
#'
#' @return A character vector of colors representing the gradient based on the input numeric values.
#'
#' @examples
#' # Generate colors for a set of numbers including zero
#' numbers <- c(-50, -20, 0, 20, 50)
#' colors <- createLinearColors(numbers, withZero = TRUE, maximum = 100)
#'
#' @export
createLinearColors <- function(numbers, withZero = T, maximum = 100, my_colors = c("royalblue3","white","red"))
{
  if (min(numbers, na.rm = T) > 0)
  {
    if(withZero)
    {
      numbers <- c(0,numbers)
    }
    myPalette <- colorRampPalette(my_colors[c(2,3)])
    myColors <- myPalette(maximum)
  }
  else
  {
    if (max(numbers, na.rm = T) < 0)
    {
      if(withZero)
      {
        numbers <- c(0,numbers)
      }
      myPalette <- colorRampPalette(my_colors[c(1,2)])
      myColors <- myPalette(maximum)
    }
    else
    {
      myPalette_pos <- colorRampPalette(c("white","red"))
      myPalette_neg <- colorRampPalette(c("royalblue3","white"))
      npos <- length(numbers[numbers >= 0]) + 1
      nneg <- length(numbers[numbers <= 0]) + 1
      
      myColors_pos <- myPalette_pos(npos)
      myColors_neg <- myPalette_neg(nneg)
      
      #print(myColors_neg)
      #print(myColors_pos)
      
      myColors <- c(myColors_neg[-(nneg)], myColors_pos[-1])
    }
  }
  return(myColors)
}

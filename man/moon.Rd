% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decoupleRnival.R
\name{moon}
\alias{moon}
\title{moon}
\usage{
moon(
  upstream_input = NULL,
  downstream_input,
  meta_network,
  n_layers,
  n_perm = 1000,
  downstream_cutoff = 0,
  statistic = "ulm",
  return_levels = F
)
}
\arguments{
\item{upstream_input}{A named vector with up_stream nodes and their corresponding activity.}

\item{downstream_input}{A named vector with down_stream nodes and their corresponding activity.}

\item{meta_network}{A network data frame containing signed directed prior knowledge of molecular interactions.}

\item{n_layers}{The number of layers that will be propagated upstream.}

\item{n_perm}{The number of permutations to use in decoupleR's algorithm.}

\item{downstream_cutoff}{If downstream measurments should be included above a given threshold}

\item{statistic}{the decoupleR stat to consider: "wmean", "norm_wmean", or "ulm"}

\item{return_levels}{true or false, if true the layers that the protein belongs to will be returned alongside the scores}
}
\value{
A data frame containing the score of the nodes upstream of the 
downstream input based on the iterative propagation
}
\description{
Iteratively propagate downstream input activity through a signed directed network
using the weighted mean enrichment score from decoupleR package
}
\examples{
# Example input data
upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
downstream_input <- c("D" = 2, "E" = -1.5)
meta_network <- data.frame(
  source = c("A", "A", "B", "C", "C", "D", "E"),
  target = c("B", "C", "D", "E", "D", "B", "A"),
  sign = c(1, -1, -1, 1, -1, -1, 1)
)

# Run the function with the example input data
result <- moon(upstream_input, downstream_input, meta_network, n_layers = 2, statistic = "wmean")

# View the results
print(result)
}

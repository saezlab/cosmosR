% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support_decoupleR.R
\name{wide_ulm_res}
\alias{wide_ulm_res}
\title{Convert ULM Results to Wide Format}
\usage{
wide_ulm_res(ulm_result)
}
\arguments{
\item{ulm_result}{A data frame representing the ULM results with columns: source, condition, and score.}
}
\value{
A data frame in wide format where each row is a source and each column is a condition.
}
\description{
This function converts the results from a ULM analysis to a wide format data frame.
The input is a data frame with columns for source, condition, and score. The output 
is a data frame where each row represents a source and each column represents a condition,
with the corresponding scores as values.
}
\examples{
# Create a sample ULM result
ulm_result <- data.frame(source = c("A", "A", "B", "B"),
                         condition = c("cond1", "cond2", "cond1", "cond2"),
                         score = c(0.5, 0.8, 0.3, 0.7))

# Convert to wide format
wide_ulm_result <- wide_ulm_res(ulm_result)

}

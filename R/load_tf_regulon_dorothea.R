#' Retired DoRothEA regulon loader
#'
#' `load_tf_regulon_dorothea()` is retained only to give a clear migration
#' message. DoRothEA is no longer a cosmosR dependency. For the legacy
#' CARNIVAL workflow, supply a prepared CollecTRI regulon explicitly through
#' the `tf_regulon` argument of `preprocess_COSMOS_*()`. It must contain the
#' columns `tf`, `sign`, and `target`. A cached `decoupleR` CollecTRI table
#' has `source`, `target`, and `mor` columns; convert it with
#' `data.frame(tf = collectri_regulon$source, sign = collectri_regulon$mor, target = collectri_regulon$target)`.
#'
#' @param confidence Deprecated and ignored.
#' @return No value. This defunct compatibility function always errors.
#' @export
load_tf_regulon_dorothea <- function(confidence = c("A","B","C")){
    stop(
        "`load_tf_regulon_dorothea()` was retired in cosmosR 1.21.2. ",
        "Supply a prepared CollecTRI regulon to the legacy CARNIVAL ",
        "preprocessing functions using columns `tf`, `sign`, and `target`.",
        call. = FALSE
    )
}

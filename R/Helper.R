#' @title Optimal Accuracy Computing
#'
#' @description
#' Computing AUC
#'
#' @param pred Predicted labels vector, as returned by a classifier
#' @param label Ground truth (correct) 0-1 labels vector
#' @return AUC
#' @keywords internal
#' @export
Optimal_Accuracy <- function(pred, label) {
  # Check
  stopifnot(length(pred) == length(label))
  stopifnot(all(is.finite(pred)))
  stopifnot(all(is.finite(label)))
  # Computing
  ModelMetrics::auc(label, pred)
}

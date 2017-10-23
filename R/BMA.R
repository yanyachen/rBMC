#' @title Bayesian Model Averaging
#'
#' @description
#' Computing Ensemble Weights using Bayesian Model Averaging.
#'
#' @param Model_List A list of numeric probabilities vector
#' @param Target Ground truth (correct) 0-1 labels vector
#' @return Ensemble Weights
#' @export
Bayesian_Model_Averaging <- function(Model_List, Target) {
  # Preparation
  M <- length(Model_List)
  N <- length(Target)
  # Uniform Prior
  Prior <- rep(1 / M, M)
  # Pre - allocation
  Log_Likelihood <- numeric(M)
  Z <- -Inf
  # Iteration
  for (i in 1:M) {
    A <- Optimal_Accuracy(Model_List[[i]], Target)
    Log_Likelihood[i] <- N * (A * log(A) + (1 - A) * log(1 - A))
    Z <- max(Z, Log_Likelihood[i])
  }
  # Computing Weight
  Weight <- Prior * exp(Log_Likelihood - Z)
  Weight0 <- Weight / sum(Weight)
  # Return
  return(Weight0)
}

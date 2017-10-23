#' @title Bayesian Model Combination Naive
#'
#' @description
#' Computing Ensemble Weights using Bayesian Model Combination Naive Version.
#'
#' @param Model_List A list of numeric probabilities vector
#' @param Target Ground truth (correct) 0-1 labels vector
#' @param Sample_Size Number of weights to be random sampled
#' @return Ensemble Weights
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @export
Bayesian_Model_Combination_Naive <- function(Model_List, Target, Sample_Size) {
  # Preparation
  M <- length(Model_List)
  N <- length(Target)
  # Uniform Prior
  Prior <- rep(1 / Sample_Size, Sample_Size)
  # Pre - allocation
  Log_Likelihood <- numeric(Sample_Size)
  Z <- -Inf
  # Sampling and Computing Prediction
  Model_Mat <- Model_List %>%
    data.table::as.data.table(.) %>%
    as.matrix(.)
  V_Mat <- foreach::foreach(i = 1:Sample_Size) %do% {
    V <- -log(stats::runif(M, min = 0, max = 1)) %>%
      magrittr::divide_by(., sum(.))
  } %>%
    data.table::as.data.table(.) %>%
    as.matrix(.)
  Pred_Mat <- Model_Mat %*% V_Mat
  # Iteration
  for (i in 1:Sample_Size) {
    A <- Optimal_Accuracy(Pred_Mat[,i], Target)
    Log_Likelihood[i] <- N * (A * log(A) + (1 - A) * log(1 - A))
    Z <- max(Z, Log_Likelihood[i])
  }
  # Computing Weight
  Weight <- Prior * exp(Log_Likelihood - Z) %>%
    magrittr::divide_by(., sum(.))
  Weight0 <- (V_Mat %*% Weight) %>%
    magrittr::divide_by(., sum(.)) %>%
    as.numeric(.)
  # Return
  return(Weight0)
}



#' @title Bayesian Model Combination Helper Function
#'
#' @description
#' Computing Ensemble Weights in one round for Bayesian Model Combination.
#'
#' @param Model_List A list of numeric probabilities vector
#' @param Target Ground truth (correct) 0-1 labels vector
#' @param Prior Prior for the ensemble weights
#' @param Sample_Size Number of weights to be random sampled
#' @return Ensemble Weights
#' @keywords internal
#' @importFrom magrittr %>%
#' @export
Bayesian_Model_Combination_Helper <- function(Model_List, Target, Prior = NULL, Sample_Size) {
  # Preparation
  M <- length(Model_List)
  N <- length(Target)
  # Prior
  Weight <- Prior
  Weight_Sum <- sum(Weight)
  # Pre - allocation
  Log_Likelihood <- numeric(Sample_Size)
  Z <- -Inf
  # Sampling and Computing Prediction
  Model_Mat <- Model_List %>%
    data.table::as.data.table(.) %>%
    as.matrix(.)
  V_Mat <- gtools::rdirichlet(Sample_Size, Weight) %>% t(.)
  Pred_Mat <- Model_Mat %*% V_Mat
  # Iteration
  for (i in 1:Sample_Size) {
    A <- Optimal_Accuracy(Pred_Mat[,i], Target)
    if (A == 1) {
      return(V_Mat[,i])
    } else if (A == 0) {
      next
    } else {
      Log_Likelihood[i] <- N * (A * log(A) + (1 - A) * log(1 - A))
    }
    if (Log_Likelihood[i] > Z) {
      if (i != 1L) {
        Weight <- Weight * exp(Z - Log_Likelihood[i])
      }
      Z <- Log_Likelihood[i]
    }
    W <- exp(Log_Likelihood[i] - Z)
    Weight <- Weight * Weight_Sum / (Weight_Sum + W) + W * V_Mat[, i]
    Weight_Sum <- Weight_Sum + W
  }
  # Computing Weight
  Weight0 <- Weight / sum(Weight)
  # Return
  return(Weight0)
}



#' @title Bayesian Model Combination Function
#'
#' @description
#' Computing Ensemble Weights using Bayesian Model Combination iteratively.
#'
#' @param Model_List A list of numeric probabilities vector
#' @param Target Ground truth (correct) 0-1 labels vector
#' @param Sample_Size Number of weights to be random sampled
#' @param Sample_Time Number of times for weights to be random sampled
#' @param Prior Prior for the ensemble weights
#' @param Verbose Whether or not to print progress.
#' @return Ensemble Weights
#' @export
Bayesian_Model_Combination <- function(Model_List, Target, Sample_Size, Sample_Time, Prior = NULL, Verbose = FALSE) {
  # Prior
  if (is.null(Prior)) {
    M <- length(Model_List)
    Weight <- rep(1 / M, M)
  } else {
    Weight <- Prior
  }
  # Iteration
  for (i in seq_len(Sample_Time)) {
    if (Verbose & i == 1) {
      cat("Weight Initialization:", "\n")
      cat(round(Weight, 3), "\n")
    }
    Weight <- Bayesian_Model_Combination_Helper(Model_List, Target, Weight, Sample_Size)
    if (Verbose) {
      cat("Round:", i, "\n")
      cat(round(Weight, 3), "\n")
      cat(Optimal_Accuracy(as.numeric(as.matrix(data.table::as.data.table(Model_List)) %*% Weight),
                           Target),
          "\n")
    }
  }
  return(Weight)
}

utils::globalVariables(c("."))

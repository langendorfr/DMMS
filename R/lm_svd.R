#' @title SVD solution to the Linear Least Squares Problem
#'
#' @description Solution to the linear least squares problem using Singular Value Decomposition (SVD). 
#'
#' @param y Target time series at time t+1.
#' 
#' @param x Time series at time t.
#' 
#' @param weights Weights for the exponentially-weighted linear model.
#'
#' @return Values in the Jacobian.
#' 
#' @references Deyle, E. R., May, R. M., Munch, S. B., & Sugihara, G. (2016, January). Tracking and forecasting ecosystem interactions in real time. In Proc. R. Soc. B (Vol. 283, No. 1822, p. 20152258). The Royal Society.
#'
#' @examples 
#' x <- matrix(sample(1:10, size = 30, replace = TRUE), nrow = 10, ncol = 3)
#' y <- 2*x[,1] + 0.25*x[,2] + rnorm(n = 10, mean = 0, sd = 1)
#' weights <- apply(x, MARGIN = 1, FUN = var)
#' 
#' lm_svd(y, x, weights)
#' 
#' @export

lm_svd <- function(y, x, weights, intercept = TRUE) 
{
  # Assign all ones for the weights if none are supplied
  if (missing(weights)) {
    weights <- rep(1, nrow(x))
  }
  
  # Process weights
  weights_processed <- sqrt(weights/sum(weights)) #/sum(weights)
  
  # Add column of ones for the constant term in the linear model unless intercept = FALSE
  if (intercept == TRUE) {
    A <- cbind(1, x) * weights_processed
  } else {
    A <- x * weights_processed
  }
  
  # Singular Value Decomposition
  A_svd <- svd(A)
  
  # Calculate Sigma^-1
  s <- A_svd$d
  s_inv <- matrix(0, nrow = ncol(A), ncol = ncol(A))
  for(i in seq_along(s)){
    if(s[i] >= max(s) * 1e-5) { # Remove small singular values
      s_inv[i, i] <- 1/s[i]
    }
  }
  
  # Exponentially-weighted regression using singular value decomposition
  coeff <- t(A_svd$v %*% s_inv %*% t(A_svd$u) %*% (weights_processed * y))
  
  
  if (intercept == TRUE) {
    if (is.null(colnames(x)) == TRUE) {
      colnames(coeff) <- c("Constant", colnames(data.frame(x)))  
    } else {
      colnames(coeff) <- c("Constant", colnames(x))
    }
  } else {
    if (is.null(colnames(x)) == TRUE) {
      colnames(coeff) <- colnames(data.frame(x))  
    } else {
      colnames(coeff) <- colnames(x)
    }
  }
  
  return(coeff)
}

#' @title SVD Solution to the Linear Least Squares Problem
#'
#' @description Solution of the linear least squares problem using Singular Value Decomposition (SVD). 
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
#' lm_svdsolve(block[libs, 1], block[libs, 2:dim(block)[2]], weights)
#' 
#' @export

lm_svdsolve <- function(y, x, weights) 
{
  
  # Add column of ones for the constant term in the linear model
  A <- cbind(1, x) * weights
  A_svd <- svd(A)
  
  s <- A_svd$d
  s_inv <- matrix(0, nrow = (ncol(x)+1), ncol = (ncol(x)+1))
  for(i in seq_along(s)){
    # if(s[i] >= max(s) * 1e-5) # Remove small singular values
      s_inv[i, i] <- 1/s[i]
  }
  
  # Exponentially-weighted regression using singular value decomposition
  coeff <- t(A_svd$v %*% s_inv %*% t(A_svd$u) %*% (weights * y))
  colnames(coeff) <-c("Constant", colnames(x))
  
  return(coeff)
}

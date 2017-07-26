#' @title Tracking and Forecasting Interactions in Real Time
#'
#' @description Calculates sequential Jacobians on time series data using the S-map method developed and detailed in the 2016 Proceedings of the Royal Society B paper titled Tracking and Forecasting Ecological Interactions in Real Time by Deyle, May, Munch, and Sugihara.
#'
#' @param time_series A data frame where each column is a time series. 
#' 
#' @param target Which time series (column) to calculate the sequential interactions of. 
#' 
#' @param date Optional vector of dates. 
#' 
#' @param theta Parameter to tune the relationship between distance and weight in the linear model. Defaults to 8 which was used by Deyle, May, Munch, and Sugihara.
#'
#' @return A data frame where each column is a row of the sequential Jacobian.
#' 
#' @references Deyle, E. R., May, R. M., Munch, S. B., & Sugihara, G. (2016, January). Tracking and forecasting ecosystem interactions in real time. In Proc. R. Soc. B (Vol. 283, No. 1822, p. 20152258). The Royal Society.
#'
#' @examples
#' data <- matrix(runif(300, 0, 1), nrow = 100, ncol=3))
#' 
#' track(data, target = 1)
#' 
#' @export

track <- function(time_series, target, date, theta = 8)
{
  
  # Names for the output sequential Jacobian
  coeff_names <-sapply(colnames(time_series), function(x) paste("d", colnames(time_series)[target], "/d", x, sep =""))
  
  # Combine target time series at time t+1 with all time series at time t
  block_raw <- cbind(time_series[2:nrow(time_series), target], time_series[1:(nrow(time_series)-1), ])
  block <- as.data.frame(apply(block_raw, 2, function(x) (x-mean(x)) / (sd(x)+1e-5) )) # 1e-5 is added to prevent divide by zero errors
  
  # Full library of time steps
  lib <- 1:nrow(block)  
  
  # Full set of prediction time steps
  pred <- 1:nrow(block) 
  
  # Output sequential Jacobian
  coeff <- as.data.frame(matrix(0, nrow = length(pred), ncol = ncol(time_series)))
  colnames(coeff) <- coeff_names

  # Exponentially-weighted linear model for each time point from t to t_max-1  
  for (ipred in 1:length(pred)) {
    
    # Target time point is excluded from the fitting procedure
    libs = lib[-pred[ipred]]
    
    # Calculate weights
    q <- matrix(as.numeric(block[pred[ipred], 2:dim(block)[2]]), ncol = ncol(time_series), nrow = length(libs), byrow = TRUE)
    distances <- sqrt(rowSums((block[libs, 2:dim(block)[2]] - q)^2))
    dbar <- mean(distances)
    weights <- exp(-theta * distances / dbar)
    
    # Regression using singular value decomposition and the calculated weights
    svd_fit <- DMMS::lm_svdsolve(block[libs, 1], block[libs, 2:dim(block)[2]], weights)
    
    # Add sequential Jacobians to the growing output, ignoring the constant term in the linear model
    coeff[ipred, ] <- svd_fit[-1]
  }
  
  # Combine output with the observation dates if supplied, or a simple counter ID if not
  if (missing(date)) {
    coeff <- cbind(pred, coeff)
    colnames(coeff)[1] <- "ID"    
  } else {
    coeff <- cbind(date[1:(length(date)-1)], coeff)
    colnames(coeff)[1] <- "Date" 
  }
    
  return(coeff)
}

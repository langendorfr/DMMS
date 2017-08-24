#' @title Tracking and Forecasting Interactions in Real Time
#'
#' @description Calculates sequential Jacobians on time series data using the S-map method developed and detailed in the 2016 Proceedings of the Royal Society B paper titled Tracking and Forecasting Ecological Interactions in Real Time by Deyle, May, Munch, and Sugihara.
#'
#' @param time_series A data frame where each column is a time series. 
#' 
#' @param target Which time series (column) to calculate the sequential interactions of. 
#' 
#' @param dates Optional vector of dates. 
#' 
#' @param theta Parameter to tune the relationship between distance and weight in the linear model. Defaults to 8 which was used by Deyle, May, Munch, and Sugihara.
#'
#' @param manual_block Defaults to FALSE. If your data is stacked spatial replicates, where you have already made the first column for time 1+1 and all others for time t.
#'
#' @param causal_probabilities Causal probabilities, or 1 - p-values from a test for causality, for each interactor (column in the data). 
#' 
#' @param causal_iterations Defaults to 100. Determines how many times the causal filter is applied. For large datasets this can cause slowdowns which is why the default is only 100. 
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

track <- function(time_series, target, dates, theta = 8, manual_block = FALSE, causal_probabilities, causal_iterations = 1e2, causal_ID)
{
  
  # Names for the output sequential Jacobian
  if (manual_block == TRUE) {
    coeff_names <- sapply(colnames(time_series)[2:ncol(time_series)], function(x) paste("d", colnames(time_series)[1], "/d", x, sep ="")) 
  } else {
    coeff_names <- sapply(colnames(time_series), function(x) paste("d", colnames(time_series)[target], "/d", x, sep =""))
  }

  # Combine target time series at time t+1 with all time series at time t, unless block == TRUE
  if (manual_block == TRUE) {
    block_raw <- time_series
  } else {
    block_raw <- cbind(time_series[-1, target], time_series[-nrow(time_series), ])
  }  
  
  # Causally filter the data. Note this defaults to FALSE.
  if (missing(causal_probabilities)) {
    block_processed <- block_raw  
  } else {
    dat <- causal_filter(as.matrix(block_raw[, -1]), probabilities = causal_probabilities, iterations = causal_iterations, ID = causal_ID)
    block_processed <- cbind(rep(as.matrix(block_raw[, 1]), causal_iterations), dat$data)  
  }
  
  
  
  
  
  
  # "All time series were normalized to have a mean of 0 and standard deviation of 1." - Deyle, May, Munch, and Sugihara
  block <- as.data.frame(apply(block_processed, 2, function(x) (x-mean(x)) / (sd(x)+1e-5) )) # 1e-5 is added to prevent divide by zero errors
  
  # Full library of time steps
  lib <- 1:nrow(block)  
  
  # Full set of prediction time steps
  pred <- 1:nrow(block) 
  
  # Output sequential Jacobian
  if (manual_block == TRUE) {
    coeff <- as.data.frame(matrix(0, nrow = length(pred), ncol = ncol(time_series) - 1))
  } else {
    coeff <- as.data.frame(matrix(0, nrow = length(pred), ncol = ncol(time_series)))
  }
  
  colnames(coeff) <- coeff_names

  # Exponentially-weighted linear model for each time point from t to t_max-1  
  for (ipred in 1:length(pred)) {
    
    # Target time point is excluded from the fitting procedure
    libs = lib[-pred[ipred]]
    
    # If a causal filter was supplied determine which data points will be included in calculating the weights for the linear models
    if (missing(causal_probabilities)) {
      weights_filter <- matrix(1, nrow = nrow(block) - 1, ncol = ncol(block) - 1)
    } else {
      weights_filter <- dat$filter[-pred[ipred], ]
    }
    
    # Calculate weights
    q <- matrix(as.numeric(block[pred[ipred], 2:ncol(block)]), ncol = ncol(coeff), nrow = length(libs), byrow = TRUE)
    
    block_filter <- block[libs, 2:ncol(block)] * weights_filter
    q_filter <- q * weights_filter
    
    distances <- sqrt(rowSums((block_filter - q_filter)^2))
    dbar <- mean(distances)
    weights <- exp(-theta * distances / dbar)
    
    # Regression using singular value decomposition and the calculated weights
    y <- block[libs, 1]
    x <- block[libs, 2:ncol(block)]
    # svd_fit <- DMMS::lm_svd(y, x, weights)
    svd_fit <- lm_svd(y, x, weights)
    
    # Add sequential Jacobians to the growing output, ignoring the constant term in the linear model
    coeff[ipred, ] <- svd_fit[-1]
  }
  
  
  # Needed to find the summarized Jacobians for each true data point if causal filtering was used
  if (missing(causal_probabilities)) {
    data_size <- nrow(time_series)
  } else {
    data_size <- nrow(time_series)/causal_iterations
  }
  
  # Return the average for each data point if causal filtering was used
  if (missing(causal_probabilities)) {
  } else {
    coeff_full <- coeff
    coeff <- data.frame(matrix(NA, nrow = data_size, ncol = ncol(coeff)))
    colnames(coeff) <- colnames(coeff_full)
    
    for (r in 1:data_size) {
      subset <- seq(from = r, to = nrow(coeff_full), by = data_size)
      coeff[r, ] <- apply(coeff_full[subset,], MARGIN = 2, FUN = mean)
    }
    
  }
  
  
  
  
  # Combine output with the observation dates if supplied, or a simple counter ID if not
  if (missing(dates)) {
    coeff <- cbind(1:data_size, coeff)
    colnames(coeff)[1] <- "ID"    
  } else {
    coeff <- cbind(dates[1:(length(dates)-1)], coeff)
    colnames(coeff)[1] <- "Date" 
  }
  

    
  return(coeff)
}

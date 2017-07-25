#' @title Tracking and Forecasting Interactions in Real Time
#'
#' @description Calculates sequential Jacobians on time series data using the S-map method developed and detailed in the 2016 Proceedings of the Royal Society B paper titled Tracking and Forecasting Ecological Interactions in Real Time by Deyle, May, Munch, and Sugihara.
#'
#' @param time_series A data frame where each column is a time series. 
#' 
#' @param target Which time series (column) to calculate the sequential interactions of.
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

track <- function(time_series, target)
{

  

  return()
}

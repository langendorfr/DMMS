#' @title Causal Filtering
#'
#' @description Filters data with supplied probabilities of causality. This can serve a priori in place of post-hoc model selection. 
#' 
#' @param data Data.
#' 
#' @param probabilities Causal probabilities, or 1 - p-values from a test for causality, for each interactor (column in the data). 
#' 
#' @param iteractions Number of times to apply the causal filter.
#'
#' @return Values in the Jacobian.
#'
#' @examples 
#' 
#' @export

causal_filter <- function(data, probabilities, iterations, ID)
{
  # IDs allow for the filter to be applied to stacked spatial replicates
  if (missing(ID)) {
    # ID = rep(1, nrow(data)*iterations)
    ID = rep(1, nrow(data))
  }
  
  # Initialize filtered data and record of which data was filtered
  new_data <- data.frame(matrix(NA, nrow = nrow(data)*iterations, ncol = ncol(data)))
  causal_filter <- data.frame(matrix(NA, nrow = nrow(data)*iterations, ncol = ncol(data)))
  
  for (i in 1:iterations) {
    print(c(i, iterations, i/iterations))
    # Every data point should be included the same number of times, so iterate in units the size of the original data
    rows <- (i*nrow(data) - nrow(data) + 1):(i * nrow(data))
    
    # Keep track of the filter so weights for the linear model can be calculated only using true data points
    filter_check <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))   
    
    # Apply the filter using the probabilities in probabilities  
    data_filtered <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data))) 
    for (r in 1:nrow(data)) {
      for (c in 1:ncol(data)) {
        if (runif(1) <= probabilities[c]) {
          data_filtered[r, c] <- data[r, c]
          filter_check[r, c] <- 1
        } else {
          ID_rows <- which(ID == ID[r]) # Only pick randomly from the same spatial replicate (automatically includes everything if there are no such replicates)
          data_filtered[r, c] <- sample(as.matrix(data[ID_rows , c]), size = 1) # Note that the original data point can be randomly selected
          filter_check[r, c] <- 0
        }
      }
    }
    
    # Add the causally-filtered data to the growing dataset
    new_data[rows, ] <-  data_filtered
    causal_filter[rows, ] <- filter_check
  }
  
  return(list(data = new_data, filter = causal_filter))
  
}

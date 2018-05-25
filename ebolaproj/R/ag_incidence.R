#' ag_incidence
#'
#' agregate an incidence for periods of delta days, it cuts the most recent days if incidence has not the exact number of days require.
#'
#' @param I matrix of incidence per unit time (nrow: locations, ncol: nb of time points).
#' @param delta integer, period over which to aggregate, e.g. if I is daily incidence, and delta is 7, the function return weekly incidence.
#'
#' 
#' @export
#' 
# agregate an incidence for periods of delta days, it cuts the most recent days if incidence has not the exact number of days require
ag_incidence <- function(I,delta){
  n <- nrow(I)
  n2 <- floor(n/delta)
  I_delta <- matrix(NA,n2,ncol(I))
  for (i in 2:ncol(I)){
    temp <- matrix(I[1:(n2*delta),i],delta,n2,byrow = FALSE)
    I_delta[,i] <- colSums(temp)
  }
  return(I_delta)
}

#' You can also document internal package function with roxygen
#'
#' Just make sure you add 'noRd' so no latex help files are being created.
#' @noRd
some_internal_function <- function() {
  NA
}


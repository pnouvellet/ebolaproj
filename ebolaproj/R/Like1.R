#' get likelihood  
#'
#' get likelihood  
#' internal to the MCMC
#' 
#' @param lambda: 'force of infection' matrix (incidence weighted by serial interval),
#'                  column number of days, row: number of locations   
#'                  
#' @param I matrix of observed incidence, same dimension as lambda
#' 
#' @param R0 vector of reproduction numbers per locations
#'
#' @export L log likelihood
# 

Like1<-function(lambda,I,R0){
  R <- R0%*%matrix(1,1,ncol(I))
  L <- sum(rowSums(-R*lambda+I*log(R*lambda),na.rm=TRUE),na.rm=TRUE) # poisson likelihood (or bits we are interested in!)
  return(L)
}

#' You can also document internal package function with roxygen
#'
#' Just make sure you add 'noRd' so no latex help files are being created.
#' @noRd
some_internal_function <- function() {
  NA
}


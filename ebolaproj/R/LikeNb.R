#' get likelihood  
#'
#' get NegBin likelihood  
#' internal to the MCMC
#' 
#' @param lambda: 'force of infection' matrix (incidence weighted by serial interval),
#'                  column number of days, row: number of locations   
#'                  
#' @param I matrix of observed incidence, same dimension as lambda
#' 
#' @param R0 vector of reproduction numbers per locations
#'
#' @details  L log likelihood
#' @export
#' 

LikeNb<-function(lambda,I,R0,over_disp){
  R <- R0%*%matrix(1,1,ncol(I))
  # L <- sum(rowSums(-R*lambda+I*log(R*lambda),na.rm=TRUE),na.rm=TRUE) # poisson likelihood (or bits we are interested in!)
  L <- sum(dnbinom(x = I, size = over_disp, mu = R*lambda, log = TRUE) , na.rm = TRUE)
  return(L)
}

#' You can also document internal package function with roxygen
#'
#' Just make sure you add 'noRd' so no latex help files are being created.
#' @noRd
some_internal_function <- function() {
  NA
}


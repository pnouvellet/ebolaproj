#' SerialInterval
#'
#' This is an example of how to create and document exported functions.
#'
#' @param input you should always document the paramters.
#'              Including the expected data type.
#'
#' @export
#' 
# SI
SI_gamma_dist_EpiEstim <- function(mu,cv,SItrunc){
  SI_Distr <- sapply(0:SItrunc, function(e) DiscrSI(e,mean_SI,mean_SI*CV_SI) )
  SI_Distr <- SI_Distr / sum(SI_Distr)
  return(list(dist = SI_Distr, SItrunc = SItrunc))
}

#' You can also document internal package function with roxygen
#'
#' Just make sure you add 'noRd' so no latex help files are being created.
#' @noRd
some_internal_function <- function() {
  NA
}


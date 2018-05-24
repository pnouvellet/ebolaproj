#' Project forward
#'
#' Project forward according to model
#' 
#' 
#' @param Results list conatining the results of the MCMC (e.g. from MCMC_RtProj function), contains
#'                  'theta' a matrix for the posterior samples of Rt and initial conditions (ncol: twice number of locations)
#'                  with Rts in first column followed by initial conditions.   
#'                  
#' @param Nsim integer, number of simulation - has to be smaler than number of posterior samples
#' 
#' @param week_forward integer, number of week simulated forward - this include the time window used for inference, so
#'                  if 2 weeks were used to infer Rt, and we want projection 2 weeks after that, week_forward=4
#'                  
#' @param N_geo integer, number of locations
#' 
#' @param SI serial interval, output from SI_gamma_dist_EpiEstim function
#'
#' @export I_predict an array of size [Nsim,N_geo,7*weekforward] of simulated incidences
# 
# agregate an incidence for periods of delta days, it cuts the most recent days if incidence has not the exact number of days require
# 
Proj_Pois <- function(Results,Nsim,week_forward,N_geo,SI){
  
  if (Nsim<=nrow(Results$theta)) warning('Nsim must be smaller than size of posterior samples')

  # allocate output
  I_predict <- array(data = 0, dim = c(Nsim,N_geo,7*week_forward+100))
  
  fR <- sample(x = 1:nrow(Results$theta), size = Nsim, replace = FALSE)   # samples for the posterior
  R0 <- Results$theta[fR,1:N_geo]                                         # R samples
  Ini <- Results$theta[fR,(N_geo+1):(2*N_geo)]                            # initial conditions samples
  if (N_geo==1){
    R0 <- matrix(R0,Nsim)
    Ini <- matrix(Ini,Nsim)
  }
  ws <- rev(SI$dist)                                                      # reversed serial interval
  
  # reconstruct initial conditions prior to time window 
  for (k in 1:Nsim){                                                      # for each simulations
    I0 <- matrix(0,N_geo,100)             # allocate
    I0[,1] <- Ini[k,]                     # first value = sample of initial conditions posterior
    for (i in 2:100){
      f <- max(c(1,(i-SI$SItrunc)))       # first few days, serial interval is left censored
      I0[,i] <- R0[k,]*(I0[,f:i]%*%ws[((SI$SItrunc+1)-(i-f)):(SI$SItrunc+1)])  # mean expected number of cases
    }        
    
    I=cbind(I0,matrix(0,N_geo,7*week_forward))      
    # for the time window and beyond get the mean expect number of cases day after day and draw from Poisson
    for (i in (100+1):ncol(I)){
      lambda=I[,(i-SI$SItrunc):i]%*%ws
      I[,i]=rpois(N_geo,R0[k,]*lambda)
    }
    I_predict[k,,] <- I
  }
  return(I_predict)

  
}

#' You can also document internal package function with roxygen
#'
#' Just make sure you add 'noRd' so no latex help files are being created.
#' @noRd
some_internal_function <- function() {
  NA
}


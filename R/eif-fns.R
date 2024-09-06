# program: eif-fns.R
# purpose: calculates efficient influence function for mediated pcs
# author: max rubinstein
# date modified: 8/26/24

#' CalculateEIF: calculates efficient influence function for mediated probabilities
#' of causation
#'
#' @param data a dataframe containing treatment variable A, mediator M, outcome Y,
#' and mediators; outcome regression functions labeled mu00, mu01, mu10, mu11 
#' (where the format is muam for values a, m); mediator probabilities gamma1
#' and gamma0 (where the format is gammaa for values a); and propensity scores pi
#' 
#' @note the input nuisance functions can either be estimates from some other model
#' (e.g. parametric or nonparametric models); or be the true functions from some
#' simulation. this function is agnostic to this choice 
#'
#' @return dataframe with the efficient influence functions for the mediated
#' probabilities of causation (total, indirect, direct) along with various
#' subcomponents of these functions
CalculateEIF <- function(data) {

  data$eif.mu00x = with(data, ((1 - A) * (1 - M) / ((1 - gamma0) * (1 - pi))) * (Y - mu00))
  data$eif.mu01x = with(data, ((1 - A) * M / (gamma1 * (1 - pi))) * (Y - mu01))
  data$eif.mu10x = with(data, (A * (1 - M) / ((1 - gamma1) * pi)) * (Y - mu10))
  data$eif.mu11x = with(data, (A * M / (gamma1 * pi)) * (Y - mu11))
  
  data$eif.gamma1x = with(data, (A / pi) * (M - gamma1))
  data$eif.gamma0x = with(data, ((1 - A) / (1 - pi)) * (M - gamma0))
  
  data$mu10.ratio = with(data, mu10 / mu11)
  data$mu01.ratio = with(data, mu01 / mu11)
  data$mu00.ratio = with(data, mu00 / mu11)
  data$gamma.ratio = with(data, gamma0 / gamma1)
  
  data$prod00 = with(data, mu00.ratio * gamma.ratio)
  data$prod10 = with(data, mu10.ratio * gamma.ratio)
  data$prod01 = with(data, mu01.ratio * gamma.ratio)

  data$eif.mu00.ratio = with(data, (mu11 * eif.mu00x - mu00 * eif.mu11x) / mu11^2)
  data$eif.mu10.ratio = with(data, (mu11 * eif.mu10x - mu10 * eif.mu11x) / mu11^2)
  data$eif.mu01.ratio = with(data, (mu11 * eif.mu01x - mu01 * eif.mu11x) / mu11^2)
  
  data$eif.gamma.ratio = with(data, (gamma1 * eif.gamma0x - gamma0 * eif.gamma1x) / gamma1^2)

  data$eif.prod00 = with(data, (gamma0 / gamma1) * eif.mu00.ratio + (mu00 / mu11) * eif.gamma.ratio)
  data$eif.prod10 = with(data, (gamma0 / gamma1) * eif.mu10.ratio + (mu10 / mu11) * eif.gamma.ratio)
  data$eif.prod01 = with(data, (gamma0 / gamma1) * eif.mu01.ratio + (mu01 / mu11) * eif.gamma.ratio)
  
  data$eif.pniex  = with(data, eif.prod10 - eif.mu10.ratio - eif.gamma.ratio)
  data$eif.ptex   = with(data, eif.prod00 - eif.prod01 - eif.mu00.ratio)
  data$eif.pndex  = with(data, eif.ptex - eif.pniex)  
  
  data$pnie  = with(data, 1 - mu10.ratio - gamma.ratio + prod10)
  data$pte   = with(data, 1 - mu00.ratio - prod01 + prod00)
  data$pnde  = with(data, pte - pnie)
  
  data$eif.pnie = with(data, eif.pniex + pnie)
  data$eif.pte  = with(data, eif.ptex  + pte)
  data$eif.pnde = with(data, eif.pte - eif.pnie)
  
  return(data)
}

#' Calculate PI: calculates plug-in estimates of mediated probabilities of causation
#'
#' @param data a dataframe containing outcome regression functions labeled mu00, mu01, mu10, mu11 
#' (where the format is muam for values a, m); and mediator probabilities gamma1
#' and gamma0 (where the format is gammaa for values a)
#'
#' @note the input nuisance functions can either be estimates from some other model
#' (e.g. parametric or nonparametric models); or be the true functions from some
#' simulation. this function is agnostic to this choice 
#'
#' @return dataframe with the plug in estimates  for the mediated
#' probabilities of causation (total, indirect, direct) along with various
#' subcomponents of these functions
CalculatePI <- function(data) {
  data$mu1 = with(data, mu11 * gamma1 + mu10 * (1-gamma1))
  data$mu0 = with(data, mu01 * gamma0 + mu00 * (1-gamma0))
  data$cate.gamma = with(data, gamma1 - gamma0)
  data$cate.mu = with(data, mu1 - mu0)
    
  data$pcause.y = with(data, 1 - mu0 / mu1)
  data$pcause.m = with(data, 1 - gamma0 / gamma1)
  data$pcause.nie = with(data, (1 - mu10/mu11) * (1 - gamma0/gamma1))
  data$pcause.nte = with(data, 1 - ((mu00/mu11) * (1 - gamma0/gamma1) + (mu01/mu11) * (gamma0/gamma1)))
  data$pcause.nde = with(data, pcause.nte - pcause.nie)
  return(data)
}


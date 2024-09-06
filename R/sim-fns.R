# program: sim-funs.R
# purpose: functions to simulate data and nuisance parameter estimates to verify expected performance of projection estimators
# author: max rubinstein
# date modifid: august 26, 2024


#' GenerateData: simulate data for a population of size N  
#'
#' @param N population size
#'
#' @return dataset containing exposure A, mediator M, outcomes Y, and covariate X,
#' along with propensity score pi, mediator probabilities gamma, outcome regressions
#' mu, and the corresponding potential outcomes 
GenerateData <- function(N) {
  data = tibble(
    X = runif(N, 0, 1),
    pi = plogis(-1 - 0.5*X),
    A = rbinom(N, 1, pi),
    gamma0 = 0.25 - 0.1*X - 0.1*X^2,
    gamma1 = 0.65 - 0.1*X - 0.2*X^2,
    M0 = rbinom(N, 1, gamma0),
    M1 = if_else(M0 == 1, 1, rbinom(N, 1, gamma1)),
    M = A*M1 + (1 - A)*M0,
    mu11 = 0.50 - 0.1*X - 0.2*X^2,
    mu10 = 0.30 - 0.1*X - 0.1*X^2,
    mu01 = 0.60 - 0.1*X - 0.1*X^2,
    mu00 = 0.40 - 0.1*X - 0.2*X^2,
    Y00 = rbinom(N, 1, mu00),
    Y01 = rbinom(N, 1, mu01),
    Y10 = if_else(Y00 == 1, 1, rbinom(N, 1, mu10)),
    Y11 = if_else(Y01 == 1 | Y10 == 1, 1, rbinom(N, 1, mu11)),
    Y1 = M1*Y11 + (1-M1)*Y10,
    Y0 = M0*Y01 + (1-M0)*Y00,
    Y1.M0 = M0*Y11 + (1-M0)*Y10,
    Y0.M1 = M1*Y01 + (1-M1)*Y00,
    Y = A*Y1 + (1-A)*Y0
  ) %>%
    mutate(
      gamma1 = gamma1 * (1 - gamma0) + gamma0,
      mu11 = mu11 * ((1 - mu01) * (1 - mu10) * (1 - mu00)) +
        1 - ((1 - mu01) * (1 - mu10) * (1 - mu00)),
      mu10 = mu10 * (1 - mu00) + mu00
    )
  data  
}

#' AddNoise: simulate estimation error on nuisance functions
#'
#' @param data dataset containing true nuisance functions
#' @param mu_rate rate of convergence on outcome regression function
#' @param pi_rate rate of convergence on propensity score models
#' @param gamma_rate rate of convergence of mediator probabilities
#'
#' @return dataset with simulated estimation error on nuisance components
AddNoise <- function(data, mu_rate, pi_rate, gamma_rate) {
  n = nrow(data)
  data %>%
    mutate_at(vars(matches("mu10")), ~expit(logit(.) + rnorm(n, 1/n^mu_rate, 1/n^mu_rate))) %>%
    mutate_at(vars(matches("mu01")), ~expit(logit(.) + rnorm(n, -1/n^mu_rate, 1/n^mu_rate))) %>%
    mutate_at(vars(matches("mu11")), ~expit(logit(.) + rnorm(n, 2/n^mu_rate, 1/n^mu_rate))) %>%
    mutate_at(vars(matches("mu00")), ~expit(logit(.) + rnorm(n, -1/n^mu_rate, 1/n^mu_rate))) %>%
    mutate_at(vars(matches("pi")), ~expit(logit(.) + rnorm(n, -1/n^pi_rate, 1/n^pi_rate))) %>%
    mutate_at(vars(matches("gamma1")), ~expit(logit(.) + rnorm(n, -1/n^gamma_rate, 1/n^gamma_rate))) %>%
    mutate_at(vars(matches("gamma0")), ~expit(logit(.) + rnorm(n, 2/n^gamma_rate, 1/n^gamma_rate))) #%>%
}

#' EasySim: simulate estimation of projections of probabilities of causation
#'
#' @param popdata simulated population data
#' @param sample_size sample sizes
#' @param number_simulations number of simulations to perform
#' @param mu_rate convergence rate on outcome models
#' @param pi_rate convergence rate on propensity score models
#' @param gamma_rate convergence rate on mediator probability models
#' @param oracle if TRUE, uses the true EIF; if FALSE, uses the simulated estimated
#' EIF
#'
#' @return list of dataframes with truth, estimates, and confidence intervals 
#' for each simulation
EasySim <- function(popdata, 
                    sample_size, 
                    number_simulations, 
                    mu_rate, pi_rate, gamma_rate,
                    oracle = FALSE) {
  
  data = CalculateEIF(popdata)
  
  t0 = lm(pnie ~ X, data)
  t1 = lm(pnde ~ X, data)
  t2 = lm(pte ~ X, data)
  
  truth = map_dbl(list(t0, t1, t2), ~c(1, 0.75) %*% coef(.x))
  
  res = list()
  truth0 = truth
  if (oracle) {
    truth0 = rep(truth, 2)
  }
  
  for (i in 1:number_simulations) {
    samp = sample_n(popdata, sample_size) 
      
    if (oracle) {
      samp0 = CalculateEIF(samp)
    }

    samp1 = samp %>%
      AddNoise(mu_rate, pi_rate, gamma_rate) %>%
      CalculateEIF()
    
    m0 = lm(eif.pnie ~ X, samp1)
    m1 = lm(eif.pnde ~ X, samp1)
    m2 = lm(eif.pte ~ X, samp1)
    
    iter_list = list(m0, m1, m2)
    
    estimands = c("pnie", "pnde", "pte")
    
    if (oracle) {
      m0a = lm(eif.pnie ~ X, samp0)
      m1a = lm(eif.pnde ~ X, samp0)
      m2a = lm(eif.pte ~ X, samp0)
      
      iter_list = append(iter_list, list(m0a, m1a, m2a))
      estimands = c(estimands, "pnie-oracle", "pnde-oracle", "pte-oracle")
    }
    
    ests = map_dbl(iter_list, ~c(1, 0.75) %*% coef(.x))
    vcov.est = map(iter_list, ~sandwich::vcovHC(.x, type = "HC0"))
    varest = map_dbl(vcov.est, ~c(1, 0.75) %*% .x %*% c(1, 0.75))
    res[[i]] = tibble(
      estimand = estimands,
      ests = ests,
      varest = varest,
      lci = ests - 1.96 * sqrt(varest),
      uci = ests + 1.96 * sqrt(varest),
      truth = truth0,
      covered = if_else(lci < truth & uci > truth, 1, 0)
    )    
  }
  res
}

expit <- function(x) exp(x) / (1 + exp(x))

logit <- function(x) log(x / (1-x))

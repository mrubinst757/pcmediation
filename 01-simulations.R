source("R/sim-fns.R")
source("R/eif-fns.R")
library(kableExtra)

# Convergence checks
popdata = GenerateData(1e7)

# Simulations: bias, var, RMSE, coverage for n = 1000  
test1 = EasySim(popdata = popdata, 
               sample_size = 1000, 
               number_simulations = 1000, 
               pi_rate = 0.3, 
               mu_rate = 0.3, 
               gamma_rate = 0.3)

test2 = EasySim(popdata = popdata, 
                sample_size = 1000, 
                number_simulations = 1000, 
                pi_rate = 0.1, 
                mu_rate = 0.1, 
                gamma_rate = 0.1)

table_input = map(list(test1, test2), ~.x %>%
                    invoke(rbind, .) %>%
                    group_by(Estimand = estimand) %>%
                    summarize(Truth = mean(truth),
                              Bias = mean(ests - truth),
                              RMSE = sqrt(mean((ests - truth)^2)),
                              Coverage = mean(covered))) %>%
  map2(c(0.3, 0.1), ~mutate(.x, Alpha = .y)) %>%
  invoke(rbind, .) %>%
  nest(res = c(Bias, RMSE, Coverage)) %>%
  pivot_wider(names_from = Alpha, values_from = res) %>%
  unnest() %>%
  mutate(Estimand = case_when(
    Estimand == "pnde" ~ "Direct",
    Estimand == "pnie" ~ "Indirect",
    Estimand == "pte" ~ "Total"
  ))

table_input %>%
  slice(c(2,1,3)) %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  knitr::kable(format = "latex", 
               booktabs = TRUE, 
               escape = FALSE, 
               align = "c",
               caption = "Simulation results \\label{tab:1}",
               col.names = c("Estimand", "Truth", rep(c("Bias", "RMSE", "Coverage"), 2))) %>%
  add_header_above(header = c(" " = 2, 
                              "$\\\\alpha = 0.3$" = 3, 
                              "$\\\\alpha = 0.1$" = 3), 
                   escape = FALSE)
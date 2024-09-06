# program: plot-funs.R
# purpose: generate plots for pcmediation paper
# author: max rubinstein
# date modified: 8/26/24

# generates plots of nuisance functions
Plot0 <- function(data) {
  data %>%
    sample_n(2000) %>%
    select(gamma0, gamma1, pi, mu11, mu10, mu01, mu00, X) %>%
    gather(Function, value, -X) %>%
    mutate(Group = case_when(
      grepl('mu', Function) ~ "Outcome models",
      grepl("gamma", Function) ~ "Mediator models",
      grepl("pi", Function) ~ "Propensity score model"
    )) %>%
    mutate(Function = case_when(
      Function == "mu11" ~ "P(Y(1, 1) = 1 | X)",
      Function == "mu01" ~ "P(Y(0, 1) = 1 | X)",
      Function == "mu10" ~ "P(Y(1, 0) = 1 | X)",
      Function == "mu00" ~ "P(Y(0, 0) = 1 | X)",
      Function == "gamma1" ~ "P(M(1) = 1 | X)",
      Function == "gamma0" ~ "P(M(0) = 1 | X)",
      Function == "pi" ~ "P(A = 1 | X)"
    )) %>%
    group_by(Function) %>%
    mutate(label = if_else(
      X == max(X), Function, NA_character_
    )) %>%
    mutate(Group = factor(Group, levels = c("Outcome models", "Mediator models", "Propensity score model"))) %>%
    ggplot(aes(x = X, y = value, color = Function)) +
    geom_line() +
    theme_minimal() +
    facet_wrap(~Group) +
    ggrepel::geom_label_repel(aes(label = label), size = 5) +
    theme(legend.position="none") +
    ylab("Value") +
    theme(text = element_text(size = 20))
}

# generates plots of conditional causal effects (ATE and PC)
Plot1 <- function(data) {
  data %>%
    sample_n(2000) %>%
    select(X, matches("cate"), matches("pcause")) %>%
    gather(Function, value, -X) %>%
    mutate(Group = case_when(
      grepl("cate", Function) ~ "Average treatment effect",
      grepl("pcause\\.n", Function) ~ "Mediated probabilities of causation",
      TRUE ~ "Probability of causation"
    )) %>%
    group_by(Function) %>%
    mutate(Group = factor(Group, levels = c("Average treatment effect", 
                                            "Probability of causation", 
                                            "Mediated probabilities of causation"))) %>%
    mutate(Function = case_when(
      Function == "cate.gamma" ~ "P(M(1) - M(0) = 1 | X)",
      Function == "cate.mu" ~ "P(Y(1) - Y(0) = 1 | X)",
      Function == "pcause.m" ~ "P(M(0) = 0 | M(1) = 1, X)",
      Function == "pcause.y" ~ "P(Y(0) = 0 | Y(1) = 1, X)",
      Function == "pcause.nte" ~ "Total Med. Pr. Causation",
      Function == "pcause.nde" ~ "Direct Pr. Causation",
      Function == "pcause.nie" ~ "Indirect Pr. Causation"
    )) %>%
    mutate(label = if_else(
      X == max(X), Function, NA_character_
    )) %>%
    ggplot(aes(x = X, y = value, color = Function)) +
    geom_line() +
    theme_minimal() +
    facet_wrap(~Group) +
    ggrepel::geom_label_repel(aes(label = label), size = 5) +
    theme(legend.position="none") +
    ylab("Value") +
    theme(text = element_text(size = 20))
}

# plots the total mediated prob. causation and prob. of indirect causation
Plot2 <- function(data) {
  m0 = lm(pcause.nie ~ X, data)
  m2 = lm(pcause.nte ~ X, data)
  
  colors = RColorBrewer::brewer.pal(4, "Paired")
  
  data %>%
    sample_n(2000) %>%
    select(X, matches("pcause\\.n")) %>%
    gather(Function, value, -X) %>%
    filter(Function != "pcause.nde") %>%
    group_by(Function) %>%
    mutate(Function = case_when(
      Function == "pcause.nte" ~ "Total Med. Pr. Causation",
      Function == "pcause.nde" ~ "Pr. Direct Causation",
      Function == "pcause.nie" ~ "Pr. Indirect Causation"
    )) %>%
    mutate(label = if_else(
      X == max(X), Function, NA_character_
    )) %>%
    ggplot(aes(x = X, y = value, color = Function)) +
    geom_line(lwd = 1.1) +
    theme_minimal() +
    ggrepel::geom_label_repel(aes(label = label), size = 5) +
    theme(legend.position="none") +
    ylab("Value") +
    geom_abline(slope = m0[[1]][2], intercept = m0[[1]][1], color = colors[1], lwd = 1.1) +
    geom_abline(slope = m2[[1]][2], intercept = m2[[1]][1], color = colors[3], lwd = 1.1) +
    scale_color_manual(values = colors[c(2,4)]) +
    theme(text = element_text(size = 14))
}

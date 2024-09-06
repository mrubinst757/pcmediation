source("R/sim-fns.R")
source("R/eif-fns.R")
source("R/plot-fns.R")

# plots
data = GenerateData(2e6) %>%
  CalculatePI()

# plot 0: dgp
plot0 = Plot0(data)

ggsave("plots/plot0.png", plot0, height = 7, width = 12)

# plot 1: estimands
plot1 = Plot1(data)

ggsave("plots/plot1.png",  plot1, height = 7, width = 12)

# plot 2: projections
plot2 = Plot2(data)

ggsave("plots/plot2.png",  plot2, height = 7, width = 12)

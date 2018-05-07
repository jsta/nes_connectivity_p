# setwd("scripts")
source("99_utils.R")
source("01_prepdata.R")

# ---- median_effect_size ----
library(dplyr)
library(magrittr)

dt           <- read.csv("../figures/table_1.csv")
stan_results <- readRDS("../data/nws/ll_vollenweider.rds")

median_x     <- median(
  posterior_samples(stan_results, pars = "b_x_Intercept")$b_x_Intercept)

median_tau   <- stan_results$data %>% 
  summarize(median_tau = median(retention_time_yr))
median_tau <- median_tau$median_tau

k <- apply(
  posterior_samples(stan_results, pars = c("b_k_part1", "b_k_part2")), 
  2, median)

R1 <- round((1 - (1/(1 + k[1] * (median_tau^median_x)))), 2) 
R2 <- round(1 - (1/(1 + k[2] * (median_tau^median_x))), 2)

paste0("At median water residence time, k, and tau lakes with shorter and 
      longer link lengths respectively had a P retention of ", 
       R1, " and ", R2, ".")

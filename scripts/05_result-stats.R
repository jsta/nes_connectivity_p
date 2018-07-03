# setwd("scripts")
source("99_utils.R")
source("01_prepdata.R")

# ---- median_effect_size ----
library(dplyr)
library(magrittr)

dt           <- read.csv("../figures/table_1.csv")
stan_results <- readRDS("../data/nws/ll_vollenweider.rds")

median_x     <- median(
  brms::posterior_samples(stan_results, pars = "b_x_Intercept")$b_x_Intercept)
# lower_x     <- as.numeric(quantile(
#   brms::posterior_samples(stan_results, pars = "b_x_Intercept")$b_x_Intercept, 
#   0.025))
# higher_x     <- as.numeric(quantile(
#   brms::posterior_samples(stan_results, pars = "b_x_Intercept")$b_x_Intercept, 
#   0.975))

median_tau   <- stan_results$data %>% 
  summarize(median_tau = median(retention_time_yr))
median_tau <- median_tau$median_tau

quant_tau <- as.numeric(quantile(
  stan_results$data$retention_time_yr, c(0.025, 0.5, 0.975)))

k <- apply(
  brms::posterior_samples(stan_results, pars = c("b_k_part1", "b_k_part2")), 
  2, quantile)

p_retention <- function(k, tau, x){
  round((1 - (1/(1 + k * (tau^x)))), 3) 
}

lower_k_grid <- expand.grid(k[2:4, 1], quant_tau)
upper_k_grid <- expand.grid(k[2:4, 2], quant_tau)

lower_k_grid$p_retention <- apply(lower_k_grid, 1, 
                                  function(x) p_retention(x[1], x[2], x = median_x))
upper_k_grid$p_retention <- apply(upper_k_grid, 1, 
                                  function(x) p_retention(x[1], x[2], x = median_x))

effect_sizes <- data.frame(matrix(
  lower_k_grid$p_retention - upper_k_grid$p_retention, 
       nrow = 3, ncol = 3))
row.names(effect_sizes) <- c("low", "medium", "high")
names(effect_sizes) <- c("low", "medium", "high")

knitr::kable(effect_sizes, format = "html") %>%
  kableExtra::add_header_above(c("", "", "tau", ""))

paste0()

paste0("At median water residence time, k, and tau lakes with shorter and 
      longer link lengths respectively had a P retention of ", 
       R1, " and ", R2, ".")

# ---- number_of_study_lakes ----

if(file.exists("../data")){
  pre_path <- "../data"
}else{
  pre_path <- "data"
}

dt <- read.csv(file.path(pre_path, "nes_x_lagos-ne.csv"), stringsAsFactors = FALSE)

paste("We measured connectivity in", nrow(nes_iws), 
      "out of", nrow(dt), "lakes")

# ---- nhd_version ----

library(sf)

test_gdb <- "~/.local/share/nhdR/NHDH_DC.gdb"
# st_layers(test_gdb)
View(st_read(test_gdb, "NHDSourceCitation"))
View(st_read(test_gdb, "NHDStatus"))
View(st_read(test_gdb, "NHDMetadata"))

# ---- bayesian_r2 ----

ll_results <- readRDS("../data/nws/ll_vollenweider.rds")
global_results <- readRDS("../data/global_vollenweider.rds")

lapply(list(ll_results, global_results), brms::bayes_R2)

lapply(
  dir("../data/nws/", full.names = TRUE, include.dirs = TRUE, pattern = "vollenweider.rds"), function(x) brms::bayes_R2(readRDS(x)))

lapply(
  dir("../data/iws/", full.names = TRUE, include.dirs = TRUE, pattern = "vollenweider.rds"), function(x) brms::bayes_R2(readRDS(x)))

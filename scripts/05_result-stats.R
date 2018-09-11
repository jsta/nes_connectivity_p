# setwd("scripts")
source("99_utils.R")
source("01_prepdata.R")

# ---- median_effect_size ----
library(dplyr)
library(magrittr)

# setwd("scripts")
dt           <- read.csv("../figures/table_1.csv")

p_retention_by_grp <- function(rds_path){
  stan_results <- readRDS(rds_path)
  
  median_x     <- median(
    brms::posterior_samples(stan_results, pars = "b_x_Intercept")$b_x_Intercept)
  
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
  
  lower_k_grid <- setNames(expand.grid(k[2:4, 1], quant_tau), c("k", "tau"))
  upper_k_grid <- setNames(expand.grid(k[2:4, 2], quant_tau), c("k", "tau"))
  
  lower_k_grid$p_retention <- apply(lower_k_grid, 1, 
                                    function(x) p_retention(x[1], x[2], 
                                                            x = median_x))
  upper_k_grid$p_retention <- apply(upper_k_grid, 1, 
                                    function(x) p_retention(x[1], x[2], 
                                                            x = median_x))
  
  effect_sizes <- data.frame(matrix(
    lower_k_grid$p_retention - upper_k_grid$p_retention, 
    nrow = 3, ncol = 3))
  row.names(effect_sizes) <- c("low", "medium", "high")
  names(effect_sizes) <- c("low", "medium", "high")
  
  # print(effect_sizes)
  # knitr::kable(effect_sizes, format = "html") %>%
  #   kableExtra::add_header_above(c("", "", "tau", ""))
  
  data.frame(
    effect_sizes = c(effect_sizes$medium[3], effect_sizes$medium[1]),
    p_retention = c(upper_k_grid$p_retention[5], lower_k_grid$p_retention[5]), 
    stringsAsFactors = FALSE)
}

rds_paths <- list("../data/lc_vollenweider.rds", #1
                  "../data/md_vollenweider.rds", #2
                  "../data/iws/cd_vollenweider.rds", #3
                  "../data/iws/sd_vollenweider.rds", #4
                  "../data/iws/bf_vollenweider.rds", #5
                  "../data/iws/ll_vollenweider.rds",
                  "../data/iws/la_vollenweider.rds",
                  "../data/iws/sr_vollenweider.rds", 
                  "../data/nws/cd_vollenweider.rds", 
                  "../data/nws/ll_vollenweider.rds", 
                  "../data/nws/la_vollenweider.rds", #11
                  "../data/nws/bf_vollenweider.rds", 
                  "../data/nws/sd_vollenweider.rds", 
                  "../data/nws/sr_vollenweider.rds")

stringr::str_extract("../data/nws/sr_vollenweider.rds", 
                     "(?<=data/)(.+)(?=/.)")

res        <- lapply(rds_paths, p_retention_by_grp)
names(res) <- make.names(rds_paths)

paste0("At median water residence time, k, and tau lakes with shorter and 
      longer link lengths respectively had a P retention of ", 
       res$...data.nws.ll_vollenweider.rds$p_retention[1], " and ", res$...data.nws.ll_vollenweider.rds$p_retention[2], ".")

res_signif  <- data.frame(path = unlist(rds_paths), stringsAsFactors = FALSE) %>%
  mutate(scale = stringr::str_extract(path, "(?<=data/)(.+)(?=/.)"), 
         metric = stringr::str_extract(path, "(.{2})(?=_)"), 
         signif = unlist(lapply(res, function(x){
           !any(abs(x$effect_sizes) < 0.009)})))

variable_key <- read.csv("../scripts/table_1.csv", 
                         stringsAsFactors = FALSE)

res_signif <- left_join(res_signif, 
                        distinct(dplyr::select(variable_key, 
                                               pnames2, parameter), 
                                 pnames2, .keep_all = TRUE), 
                             by = c("metric" = "pnames2"))

write.csv(dplyr::select(res_signif, parameter, scale, signif), 
          "model_signif.csv", row.names = FALSE)

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

global_results <- brms::bayes_R2(readRDS("../data/global_vollenweider.rds"))
hier_results <- lapply(c("../data/nws/ll_vollenweider.rds", 
                         "../data/nws/cd_vollenweider.rds", 
                         "../data/lc_vollenweider.rds", 
                         "../data/md_vollenweider.rds"), 
                       function(x) brms::bayes_R2(readRDS(x)))

range(do.call("rbind", hier_results)[,1])

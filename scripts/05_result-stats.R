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
  grp_eq <- data.frame(eq = c(
    paste0("Rp = 1 - (1 / (1 + ", 
           round(lower_k_grid$k["50%"], 2), 
           "+-", 
           round((
             round(lower_k_grid$k["75%"], 6) - 
               round(lower_k_grid$k["25%"], 6)) / 2, 2),
           "t^", 
           round(median_x, 2), 
           "))"), 
  ####
  paste0("Rp = 1 - (1 / (1 + ", 
         round(upper_k_grid$k["50%"], 2), 
         "+-", 
         round((
           round(upper_k_grid$k["75%"], 6) - 
           round(upper_k_grid$k["25%"], 6)) / 2, 2),
         "t^", 
         round(median_x, 2), 
         "))")
  ), stringsAsFactors = FALSE)
  
  data.frame(
    effect_sizes = c(effect_sizes$medium[3], effect_sizes$medium[1]),
    p_retention = c(upper_k_grid$p_retention[5], lower_k_grid$p_retention[5]),
    grp_eq = grp_eq,
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

p_retention_by_grp("../data/nws/ll_vollenweider.rds")$eq

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

rds_paths <- list("../data/lc_vollenweider.rds", #1
                  "../data/md_vollenweider.rds", #2
                  "../data/iws/cd_vollenweider.rds", #3
                  "../data/iws/sd_vollenweider.rds", #4
                  "../data/iws/bf_vollenweider.rds", #5
                  "../data/iws/ll_vollenweider.rds", #6
                  "../data/iws/la_vollenweider.rds", #7
                  "../data/iws/sr_vollenweider.rds", #8
                  "../data/nws/cd_vollenweider.rds", #9
                  "../data/nws/ll_vollenweider.rds", 
                  "../data/nws/la_vollenweider.rds", #11
                  "../data/nws/bf_vollenweider.rds", 
                  "../data/nws/sd_vollenweider.rds", 
                  "../data/nws/sr_vollenweider.rds")

global_results <- brms::bayes_R2(readRDS("../data/global_vollenweider.rds"))
hier_results <- lapply(rds_paths, 
                       function(x) brms::bayes_R2(readRDS(x)))

range(do.call("rbind", hier_results)[,1])

# ---- milstead_2013 ----
library(dplyr)
library(fulltext)

milstead_path <- "data/milstead_2013.csv"
if(!file.exists(milstead_path)){
  query <- ft_search(query = "Estimating Summer Nutrient Concentrations in
                      Northeastern Lakes from SPARROW Load Predictions and Modeled
                      Lake Depth and Volume", from = "plos")
  doi <- ft_links(query)$plos$ids
  milstead <- read.csv(ft_get_si(doi, 1))
  write.csv(milstead, milstead_path, row.names = FALSE)
}
milstead <- read.csv(milstead_path, stringsAsFactors = FALSE)
milstead$Rp <- 1 - (milstead$Poutput / milstead$Pinput)

milstead <- dplyr::filter(milstead, 
                            !is.na(NLA_ID))

# quantile(milstead$hrt, c(0.3, 0.6))
milstead$hrt_cat <- cut(milstead$hrt, c(-Inf, 0.04, 0.4, Inf), 
                        right = FALSE, include.lowest = TRUE, 
                        labels = c("short", "medium", "long"))

summary(dplyr::filter(milstead, hrt_cat == "long"))

group_by(milstead, hrt_cat) %>%
  summarize(Rp_range = range(Rp, na.rm = TRUE)[1] - range(Rp, na.rm = TRUE)[2], 
            Rp_75 = quantile(Rp, c(0.75)), 
            Rp_25 = quantile(Rp, c(0.25))) %>%
  mutate(Rp_IQR = Rp_75 - Rp_25)

boxplot(Rp ~ hrt_cat, data = milstead)

# ---- loading_vs_concentration ----
# Is loading a better predictor of retention than connectivity?
# > Not according to Fig 4 but a more thorough investigation is beyond the scope of the paper
# Is loading equivalent to concentration * wrt?

setwd("scripts"); source("01.5_loaddata.R"); source("99_utils.R")

hist(nes_nws$mean_depth)
# m2 * m = volume
test <- nes_nws %>% 
  # calculate tp_in (tp concentration) from total_inflow and p_total
  calculate_tp_in() %>% # mg / L / yr
  # calculate retention time independently from total_inflow and volume
  mutate(surface_area = surface_area * 1000000) %>% # km2 to m2
  mutate(volume = surface_area * mean_depth) %>%
  mutate(volume = volume * 1000000) %>% # m3 to cm3
  mutate(volume = volume * 0.001) %>% # ml to L
  mutate(total_inflow = (total_inflow / 0.001) * 3.154e+7) %>% # cms to L/yr
  mutate(retention_time_yr_manual = (1 / total_inflow) * volume) %>%
  # calculate p_total (p loading) from concentration and wrt
  mutate(tp_loading = p_total * 1000000) # kg to mg
  
plot(test$tp, test$tp_in)
abline(0, 1)

plot(test$retention_time_yr, test$retention_time_yr_manual)
abline(0, 1)

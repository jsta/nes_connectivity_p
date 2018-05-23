# setwd("scripts")
source("99_utils.R")
source("01_prepdata.R")

# ---- lake_characteristics_table ----
# table describing basic properties of lake population
lg                 <- lagosne_load("1.087.1")
nes_nws            <- left_join(nes_nws, dplyr::select(lg$iws, lagoslakeid, iws_ha))
nes_nws$percent_ag <- nes_iws$iws_nlcd1992_pct_81 + nes_iws$iws_nlcd1992_pct_82
  
name_key <- data.frame(
  property = c(
    "tp", "chl", "secchi", # chemistry
    "p_percent_retention", "retention_time_yr", # retention
    "maxdepth", "percent_ag", "percent_urban", # features
    "iws_ha", "nws_ha"), 
  Characteristic = c(
    "Total Phosphorus (ug/L)", "Chlorophyll (ug/L)", "Secchi Depth (m)",
    "P Retention (%)", "Residence Time (yr)",
    "Maximum Depth (m)", "Agricultural Landuse (%)", "Urban Landuse %", 
    "Inter-lake Watershed Area (ha)", "Network Watershed Area (ha)"), 
  digits = c(2, 2, 2, 
             2, 2, 
             2, 2, 2, 
             0, 0))

qs            <- function(x) quantile(nes_nws[,x], c(0.5, 0.25, 0.75), na.rm = TRUE)
summary_names <- c("tp", "chl", "secchi", 
                   "p_percent_retention", "retention_time_yr",
                   "maxdepth", "percent_ag", "iws_ha", "nws_ha")

res           <- lapply(summary_names, qs)
res           <- round(data.frame(do.call("rbind", res)), 2)
res$property  <- summary_names

# unit conversions
res[res$property == "tp", 1:3] <- 1000 * res[res$property == "tp", 1:3] # mg to ug
# organization
res <- merge(res, name_key, sort = FALSE)
# rounding
res_temp <- res
res[,c("X50.", "X25.", "X75.")] <- apply(res[,c("X50.", "X25.", "X75.")], 2, as.character)
res[,c("X50.", "X25.", "X75.")] <- t(apply(res_temp[,c("X50.", "X25.", "X75.", "digits")], 1, 
                                         function(x) as.character(round(x[-4], x[4]))))
res <- res[,c(ncol(res) - 1, 2:(ncol(res) - 2))]

res <- data.frame(res[,1], res$X50., paste0(res$X25., " - ", res$X75.))
names(res) <- c("", "Mean", "IQR")

knitr::kable(res, format = 'pandoc', align = c("lll"),
             caption = "Mean and interquartile range of selected lake characteristics.")

# ---- model_results_table ----
# table showing model results

res <- read.csv("table_1.csv", stringsAsFactors = FALSE)
res <- res[, c("parameter", "scale", "conny_type", "splits", "d_k")]

key2 <- data.frame(conny_type = unique(as.character(res$conny_type)), 
                   conny_full = c("Longitudinal", "Lateral", "-"))
res <- merge(res, key2, sort = FALSE)
res <- res[order(res$d_k, decreasing = TRUE),]
res <- res[,c(2, 3, 6, 4, 5)]
# res <- dplyr::select(res, -splits)

options(knitr.kable.NA = "-")
knitr::kable(res, 
             digits = 2, row.names = FALSE, 
             col.names = c("Metric", "Scale",
                           "Connectivity Type", "Split Value", "Delta k"), 
             caption = "Classification and ranking of connectivity metrics, lake depth, and their partition split values according to median effect size.")

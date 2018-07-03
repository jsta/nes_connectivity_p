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
    "surface_area", "maxdepth", "percent_ag", "percent_urban", # features
    "iws_ha", "nws_ha"), 
  Characteristic = c(
    "Total Phosphorus (ug/L)", "Chlorophyll (ug/L)", "Secchi Depth (m)",
    "P Retention (%)", "Residence Time (yr)",
    "Lake Area (km2)", "Maximum Depth (m)", "Agricultural Landuse (%)", "Urban Landuse %", 
    "Lake Watershed Area (km2)", "Network Watershed Area (km2)"), 
  digits = c(2, 2, 2, 
             2, 2, 
             2, 2, 2, 2, 
             0, 0))

qs            <- function(x) quantile(nes_nws[,x], c(0.5, 0.25, 0.75), na.rm = TRUE)
summary_names <- c("tp", "chl", "secchi", 
                   "p_percent_retention", "retention_time_yr",
                   "surface_area", "maxdepth", "percent_ag", "iws_ha", "nws_ha")

# range(nes_nws$p_percent_retention)
# range(nes_nws$retention_time_yr)
# 0.02 * 365

res           <- lapply(summary_names, qs)
res           <- round(data.frame(do.call("rbind", res)), 2)
res$property  <- summary_names

# unit conversions
res[res$property == "tp", 1:3] <- 1000 * res[res$property == "tp", 1:3] # mg to ug
res[res$property == "iws_ha", 1:3] <- res[res$property == "iws_ha", 1:3]  / 100 # ha to km2
res[res$property == "nws_ha", 1:3] <- res[res$property == "nws_ha", 1:3]  / 100 # ha to km2

# organization
res <- merge(res, name_key, sort = FALSE)
# rounding
res_temp <- res
res[,c("X50.", "X25.", "X75.")] <- apply(res[,c("X50.", "X25.", "X75.")], 2, as.character)
res[,c("X50.", "X25.", "X75.")] <- t(apply(res_temp[,c("X50.", "X25.", "X75.", "digits")], 1, 
                                         function(x) as.character(round(x[-4], x[4]))))
res <- res[,c(ncol(res) - 1, 2:(ncol(res) - 2))]

# space formatting 
row_nchar     <- t(apply(res[,3:4], 1, nchar))
row_nchar_max <- apply(row_nchar, 2, max)
spaces <- sapply(row_nchar_max[1] - row_nchar[,1], 
                 function(x) paste(rep(" ", x), collapse = ""))
res$X25. <- paste0(res$X25., spaces)

spaces <- sapply(row_nchar_max[2] - row_nchar[,2], 
                 function(x) paste(rep(" ", x), collapse = ""))
res$X75. <- paste0(res$X75., spaces)

res        <- data.frame(res[,1], res$X50., paste0(res$X25., " - ", res$X75.))
names(res) <- c("", "Mean", "IQR")

knitr::kable(res, format = 'pandoc', align = c("lll"),
             caption = "Mean and interquartile range of selected lake characteristics.")

# ---- model_results_table ----
# table showing model results

res <- read.csv("../scripts/table_1.csv", stringsAsFactors = FALSE)
res <- res[, c("parameter", "units", "scale", "conny_type", "splits", "d_k")]

key2 <- data.frame(conny_type = unique(as.character(res$conny_type)), 
                   conny_full = c("Longitudinal", "Lateral", "-"))
res <- merge(res, key2, sort = FALSE)
res <- res[order(res$d_k, decreasing = TRUE),]
res <- res[,c(2, 3, 4, 5, 6)]

# splits
# res <- dplyr::select(res, -splits)
res$splits[!is.na(res$splits)] <- sapply(res$splits[!is.na(res$splits)], 
                                         function(x) if(x > 100){
                                           as.character(round(x, 0))
                                           }else{as.character(round(x, 2))})

options(knitr.kable.NA = "-")
knitr::kable(res, 
             digits = 2, row.names = FALSE, 
             col.names = c("Metric", "Units", "Scale", "Split Value", "Delta k"), 
             align = c("lllcc"),
             caption = "Classification and ranking of connectivity metrics, lake depth, and their partition split values according to median effect size.")

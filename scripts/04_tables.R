# setwd("scripts")
source("99_utils.R")
source("01_prepdata.R")

# ---- lake_characteristics_table ----
# table describing basic properties of lake population
nes_iws$percent_ag <- nes_iws$iws_nlcd1992_pct_81 + nes_iws$iws_nlcd1992_pct_82
  
name_key <- data.frame(
  property = c(
    "tp", "chl", "secchi", # chemistry
    "p_percent_retention", "retention_time_yr", # retention
    "maxdepth", "percent_ag", "percent_urban"), # features
  Characteristic = c(
    "Total Phosphorus (ug/L)", "Chlorophyll (ug/L)", "Secchi Depth (m)",
    "P Retention (%)", "Residence Time (yr)",
    "Maximum Depth (m)", "Agricultural Landuse (%)", "Urban Landuse %"))

qs <- function(x) quantile(nes_iws[,x], c(0.5, 0.25, 0.75), na.rm = TRUE)
summary_names <- c("tp", "chl", "secchi", 
                   "p_percent_retention", "retention_time_yr",
                   "maxdepth", "percent_ag")
res <- lapply(summary_names, qs)
res <- round(data.frame(do.call("rbind", res)), 2)
res$property <- summary_names

# unit conversions
res[res$property == "tp", 1:3] <- 1000 * res[res$property == "tp", 1:3] # mg to ug

res <- merge(res, name_key, sort = FALSE)
res <- res[,c(ncol(res), 2:(ncol(res) - 1))]

res <- data.frame(res[,1], res$X50., paste0(res$X25., " - ", res$X75.))
names(res) <- c("", "Mean", "IQR")

knitr::kable(res, format = 'pandoc', align = c("lll"),
             caption = "Mean and interquartile range of selected study lake characteristics.")

# ---- model_results_table ----
# table showing model results

name_key <- function(){
  data.frame(
    parameter  = c("Max Depth", "Closest lake distance", 
                   "Stream order ratio", "Average Link Length", 
                   "Upstream lake area", "Baseflow", 
                   "Wetland Cover",
                   "Stream density", "Lake Connection"),
    abb        = c("Max depth", "Closest lake dist.", "S. order ratio", 
                   "Avg. link len.", "Up. lake area", "Baseflow", 
                   "W. Cover", "Stream density", "Lake Con."), 
    pnames     = c("maxdepth", "cd", 
                   "sr", "linklength", 
                   "uplakearea", "baseflow", 
                   "wetland-cover", 
                   "streamdensity", "lakeconnection"), 
    pnames2 = c("md", "cd", "sr", "ll", "la", "bf", "wc", "sd", "lc"), 
    nws_names = c("maxdepth", "closest_lake_distance", "stream_order_ratio", 
                  "avg_link_length", "upstream_lakes_4ha_area_ha", "baseflow", NA, 
                  "stream_density", "lakeconnection"),
    iws_names = c("maxdepth", "closest_lake_distance", "stream_order_ratio", 
                  "link_length", NA, "hu12_baseflowindex_mean", NA, 
                  "iws_streamdensity_streams_density_mperha", "lakeconnection"),
    conny_type = c(NA, "long", "long", "long", "long", "lat", "lat", "lat", "long")
  )
}

table_splits <- function(dir, pat = "forest.rds"){
  parts  <- list.files(dir, pattern = pat, 
                       full.names = TRUE, include.dirs = TRUE)
  pnames <- sapply(parts, function(x) strsplit(basename(x), "_")[[1]][1])
  splits <- sapply(parts, function(x) get_tree_split(gettree(readRDS(x))))
  cbind(pnames, splits)
}

misc_parts <- data.frame(table_splits("../data/"), 
                         stringsAsFactors = FALSE)
iws_parts  <- data.frame(table_splits("../data/iws/"), 
                         stringsAsFactors = FALSE)
nws_parts  <- data.frame(table_splits("../data/nws/"), 
                         stringsAsFactors = FALSE)

misc_parts$scale <- "focal"
iws_parts$scale  <- "iws"
nws_parts$scale  <- "nws"

rownames(misc_parts) <- rownames(iws_parts) <- rownames(nws_parts) <- NULL

res <- rbind(misc_parts, iws_parts, nws_parts)
res$splits <- as.numeric(res$splits)
# res <- res[,1:ncol(res)]
res <- rbind(res, data.frame(pnames = "lakeconnection", splits = NA, scale = "focal"))

res <- merge(name_key(), res)

res[res$pnames == "cd", "splits"] <- c(inv_inv_closest(-1.51, nes_rf_iws),
                                       inv_inv_closest(-1.448, nes_rf_nws))
res <- res[order(res$pnames, res$scale), ]

d_k <- lapply(list("../data/lc_vollenweider.rds",
                   "../data/md_vollenweider.rds",
                   "../data/iws/cd_vollenweider.rds",
                   "../data/iws/sd_vollenweider.rds",
                   "../data/iws/bf_vollenweider.rds",
                   "../data/iws/ll_vollenweider.rds",
                   "../data/iws/sr_vollenweider.rds", 
                   "../data/nws/cd_vollenweider.rds", 
                   "../data/nws/ll_vollenweider.rds", 
                   "../data/nws/la_vollenweider.rds", 
                   "../data/nws/bf_vollenweider.rds", 
                   "../data/nws/sd_vollenweider.rds", 
                   "../data/nws/sr_vollenweider.rds"), delta_k)

d_k <- do.call("rbind", d_k)
d_k$scale <- c("focal", "focal", rep("iws", 5), rep("nws", 6))

res <- merge(res, d_k)
res <- res[order(res$d_k, decreasing = TRUE),]

# write.csv(res[, c(1,2,4,5,6,7,8,9,10)], "../figures/table_1.csv", row.names = FALSE)

# pander::panderOptions('keep.line.breaks', TRUE)
# pander::pander(res)

res <- res[, c(1,2,4,8,9,10)]
res <- res[, c(3, 2, 4, 5, 6)]

key2 <- data.frame(conny_type = unique(as.character(res$conny_type)), 
                   conny_full = c("Longitudinal", "Lateral", NA))
res <- merge(res, key2, sort = FALSE)
res <- res[order(res$d_k, decreasing = TRUE),]
res <- res[,c(2, 3, 6, 4, 5)]
res <- dplyr::select(res, -splits)

options(knitr.kable.NA = '')
knitr::kable(res, 
             digits = 2, row.names = FALSE, 
             col.names = c("Metric", "Scale",
                           "Connectivity Type", "Delta k"), 
             caption = "Ranking of connectivity metrics according to median effect size.")

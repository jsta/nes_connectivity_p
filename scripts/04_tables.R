# setwd("scripts")
source("99_utils.R")
source("01_prepdata.R")

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
    conny_type = c(NA, "long", "long", "lat", "long", "lat", "lat", "lat", "long")
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
                   conny_full = c("Lateral", "Longitudinal", NA))
res <- merge(res, key2, sort = FALSE)
res <- res[order(res$d_k, decreasing = TRUE),]
res <- res[,c(2, 3, 6, 4, 5)]

knitr::kable(res, 
             digits = 2, row.names = FALSE, 
             col.names = c("Metric", "Scale",
                           "Connectivity Type", "Split Value", "Delta k"), 
             caption = "Results table")

# ---- lake_characteristics_table ----
# table describing basic properties of lake population
nes_iws$percent_ag <- nes_iws$iws_nlcd1992_pct_81 + nes_iws$iws_nlcd1992_pct_82
  
name_key <- data.frame(property = c("maxdepth", "percent_ag", "percent_urban", 
                                    "tp", "p_percent_retention", 
                                    "chl", "secchi", "retention_time_yr"),
                       Characteristic = c("Max Depth (m)", "Ag Landuse (%)", 
                                          "% Urban", 
                                       "TP (ug/L)", "P Retention (%)", 
                                       "Chl (ug/L)", "Secchi (m)", 
                                       "Residence Time (yr)"))

qs <- function(x) quantile(nes_iws[,x], c(0.5, 0.25, 0.75), na.rm = TRUE)
summary_names <- c("maxdepth", "percent_ag", 
                   "tp", "p_percent_retention", 
                   "secchi", "chl", "retention_time_yr")
res <- lapply(summary_names, qs)
res <- round(data.frame(do.call("rbind", res)), 2)
res$property <- summary_names

# unit conversions
res[res$property == "tp", 1:3] <- 1000 * res[res$property == "tp", 1:3] # mg to ug

res <- merge(res, name_key, sort = FALSE)
res <- res[,c(ncol(res), 2:(ncol(res) - 1))]
names(res)[2:ncol(res)] <-c("Mean", "LQ", "UQ")

knitr::kable(res, format = 'pandoc', 
             caption = "Summary of study lake characteristics")

# ---- cor_mat_table ----
# correlation matrix

library(corrr)
library(superheat)

splits      <- read.csv("../figures/table_1.csv", stringsAsFactors = FALSE)

nes_iws$p_pnt_source <-  rowSums(cbind(nes_iws$p_pnt_source_muni,
                                       nes_iws$p_pnt_source_septic,
                                       nes_iws$p_pnt_source_industrial),
                                 na.rm = TRUE)

nes_iws$p_pnt_source_pct <- nes_iws$p_pnt_source / 
  (nes_iws$p_nonpnt_source + nes_iws$p_pnt_source)

lg <- lagosne_load("1.087.1")
nes_iws <- dplyr::left_join(nes_iws, select(lg$iws, lagoslakeid, iws_ha))

# format iws data
nes_iws_sub <- dplyr::select(nes_iws, -lakeconnection)
conny_cols  <- unique(splits$iws_names)[!is.na(unique(splits$iws_names))]
conny_cols  <- conny_cols[!(conny_cols %in% "lakeconnection")]

nes_iws_sub <- nes_iws_sub[,c(conny_cols, 
                              "lake_area_ha", "p_pnt_source_pct", "iws_ha")]

iws_key <- merge(data.frame(iws_names = names(nes_iws_sub)), 
      splits[,c("iws_names", "abb")], sort = FALSE)
iws_key <- rbind(iws_key, data.frame(iws_names = c("lake_area_ha", "p_pnt_source_pct", "iws_ha"), abb = c("Lake Area", "Pnt. Src. P", "IWS area")))
iws_key <- iws_key[!duplicated(iws_key),]
names(nes_iws_sub) <- iws_key$abb

# format nws data
nes_nws_sub <- dplyr::select(nes_nws, -lakeconnection)
conny_cols  <- unique(splits$nws_names)[!is.na(unique(splits$nws_names))]
conny_cols  <- conny_cols[!(conny_cols %in% "lakeconnection")]

nes_nws_sub <- nes_nws_sub[,conny_cols]
nws_key <- merge(data.frame(nws_names = names(nes_nws_sub)), 
                 splits[,c("nws_names", "abb")], sort = FALSE)
nws_key <- nws_key[!duplicated(nws_key),]
names(nes_nws_sub) <- nws_key$abb

# combine iws and nws data
nes_sub <- dplyr::bind_rows(nes_iws_sub, nes_nws_sub)

# correlation matrix 
res <- shave(correlate(nes_sub))
res <- res[apply(res[,2:ncol(res)], 1, function(x) !all(is.na(x))),]
res <- res[,c(TRUE, apply(res[,2:ncol(res)], 2, function(x) !all(is.na(x))))]
names(res)[1] <- ""

res_names_c <- names(res)
res_f <- data.frame(res)
res_names_r <- res_f$Var.1
res_f <- res_f[,-1]
names(res_f) <- res_names_c[-1]
row.names(res_f) <- res_names_r
res_f <- res_f[rev(1:nrow(res_f)),]

# names(res) <- 
#   row.names(res) <- 
#   sapply(row.names(res), function(x) paste(strwrap(x, 10), collapse="\n "))

superheat(X = res_f, bottom.label.text.angle = 90, 
          left.label.size = 0.6, bottom.label.size = 0.6, 
          heat.pal = c("#542788", "white", "#b35806"), 
          heat.pal.values = c(0, 0.6, 1))

options(knitr.kable.NA = '')
knitr::kable(res, digits = 2, 
             caption = "Pearson correlation matrix among connectivity metrics.")

# # pander::panderOptions('keep.line.breaks', TRUE)
# pander::panderOptions("round", 2)
# pander::panderOptions("missing", '')
# # pander::pander(res, 
# #       caption = "Pearson correlation matrix among connectivity metrics", 
# #       split.cells = c(11, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8))
# pander::pandoc.table(res, style = "grid", keep.line.breaks = TRUE)

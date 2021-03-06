source("scripts/99_utils.R")
source("scripts/01_prepdata.R")

name_key <- function(){
  data.frame(
    parameter  = c("Max depth", "Closest lake distance", 
                   "Stream order ratio", "Average link length", 
                   "Upstream lake area", "Baseflow", 
                   "Wetland Cover",
                   "Stream density", "Lake connection"),
    abb        = c("Max depth", "Closest lake distance", "Stream order ratio", 
                   "Average link length", "Upstream lake area", "Baseflow", 
                   "Wetland Cover", "Stream density", "Lake connection"), 
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
                  "link_length", "upstream_lakes_4ha_area_ha", "hu12_baseflowindex_mean", NA, 
                  "iws_streamdensity_streams_density_mperha", "lakeconnection"),
    conny_type = c(NA, "long", "long", "long", "long", "lat", "lat", "lat", "long"), 
    units = c("m", "m", "-", "m", "ha", "-", "-", "-", "-")
  )
}

table_splits <- function(dir, pat = "forest.rds"){
  parts  <- list.files(dir, pattern = pat, 
                       full.names = TRUE, include.dirs = TRUE)
  pnames <- sapply(parts, function(x) strsplit(basename(x), "_")[[1]][1])
  splits <- sapply(parts, function(x) get_tree_split(gettree(readRDS(x))))
  cbind(pnames, splits)
}

misc_parts <- data.frame(table_splits("data/"), 
                         stringsAsFactors = FALSE)
iws_parts  <- data.frame(table_splits("data/iws/"), 
                         stringsAsFactors = FALSE)
nws_parts  <- data.frame(table_splits("data/nws/"), 
                         stringsAsFactors = FALSE)

misc_parts$scale <- "focal"
iws_parts$scale  <- "lws"
nws_parts$scale  <- "nws"

rownames(misc_parts) <- rownames(iws_parts) <- rownames(nws_parts) <- NULL

res        <- rbind(misc_parts, iws_parts, nws_parts)
res$splits <- as.numeric(res$splits)
# res <- res[,1:ncol(res)]
res <- rbind(res, data.frame(pnames = "lakeconnection", splits = NA, scale = "focal"))

res <- merge(name_key(), res)

res[res$pnames == "cd", "splits"] <- c(inv_inv_closest(-1.51, nes_rf_iws),
                                       inv_inv_closest(-1.448, nes_rf_nws))
res <- res[order(res$pnames, res$scale), ]

d_k <- lapply(list("data/lc_vollenweider.rds",
                   "data/md_vollenweider.rds",
                   "data/iws/cd_vollenweider.rds",
                   "data/iws/sd_vollenweider.rds",
                   "data/iws/bf_vollenweider.rds",
                   "data/iws/ll_vollenweider.rds",
                   "data/iws/la_vollenweider.rds",
                   "data/iws/sr_vollenweider.rds", 
                   "data/nws/cd_vollenweider.rds", 
                   "data/nws/ll_vollenweider.rds", 
                   "data/nws/la_vollenweider.rds", 
                   "data/nws/bf_vollenweider.rds", 
                   "data/nws/sd_vollenweider.rds", 
                   "data/nws/sr_vollenweider.rds"), delta_k)

d_k <- do.call("rbind", d_k)
d_k$scale <- c("focal", "focal", rep("lws", 6), rep("nws", 6))

res <- merge(res, d_k)
res <- res[order(res$d_k, decreasing = TRUE),]

get_sample_size <- function(x){
  if(is.na(x$splits)){
    sub <- get_sub2(nes_iws, as.character(x$iws_names), x$splits)
    as.numeric(c(sub[which.min(sub)], sub[which.max(sub)]))
  }else{
    if(x$scale == "nws"){
      sub <- get_sub2(nes_nws, as.character(x$nws_names), x$splits)
      c(nrow(sub), nrow(nes_nws) - nrow(sub))
    }else{
      sub <- get_sub2(nes_iws, as.character(x$iws_names), x$splits)
      c(nrow(sub), nrow(nes_iws) - nrow(sub))
    }
  }
}

sample_sizes <- data.frame(suppressWarnings(do.call("rbind", 
                        lapply(seq_len(nrow(res)), 
                               function(x) get_sample_size(res[x,])))))

res[,c("n_low", "n_high")] <- sample_sizes

write.csv(res[, c(1,2,4,5,6,7,8,9,10,11,12,13)], 
          "scripts/table_1.csv", row.names = FALSE)
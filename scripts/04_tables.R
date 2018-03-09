source("99_utils.R")

# ---- table_1 ----
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

misc_parts$scale <- "misc"
iws_parts$scale  <- "iws"
nws_parts$scale  <- "nws"

rownames(misc_parts) <- rownames(iws_parts) <- rownames(nws_parts) <- NULL

res <- rbind(misc_parts, iws_parts, nws_parts)
res$splits <- as.numeric(res$splits)
# res <- res[,1:ncol(res)]

name_key <- data.frame(
  parameter  = c("Max Depth", "Closest lake distance", 
                 "Stream order ratio", "Average Link Length", 
                 "Upstream lake area", "Baseflow", 
                 "Wetland Cover",
                 "Stream density"), 
  pnames = c("maxdepth", "cd", 
                "sr", "linklength", 
                "uplakearea", "baseflow", 
                "wetland-cover", 
                "streamdensity")
)

res <- merge(name_key, res)

res[res$pnames == "cd", "splits"] <- c(inv_inv_closest(-1.51, nes_rf_iws),
                                       inv_inv_closest(-1.49, nes_rf_nws))
                                       
# tidyr::spread(res, scale, splits)                                       

write.csv(res, "../figures/table_1.csv", row.names = FALSE)

knitr::kable(res[order(res$pnames, res$scale), 2:ncol(res)], 
             digits = 2, row.names = FALSE, 
             col.names = c("Metric", "Split Value", "Scale"))

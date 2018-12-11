# setwd("scripts")
# source("99_utils.R")
# source("01.5_loaddata.R")

# ---- 04_global_vollenweider ----

plot_grid(
  full_pred_plot(nes_iws, readRDS("../data/global_vollenweider.rds")),
  part_pred_plot(nes_nws, readRDS("../data/nws/ll_vollenweider.rds"),
                 8, "B.", xl = TRUE, yl = FALSE, legend_title = "Link length \n (connectivity)"), 
  rel_widths = c(1, 0.9), 
  rel_heights = c(1, 1))

# ---- recursive_partitioning_viz ----

plot_tree(gettree(readRDS("../data/maxdepth_forest.rds")), 
          title = "Max Depth")

plot_tree(gettree(readRDS("../data/iws/cd_forest.rds")), 
          title = "Distance to Closest \n Upstream Lake (standardized)")
# inv_inv_closest(-1.51, nes_rf_iws) # 3773.614
plot_tree(gettree(readRDS("../data/nws/cd_forest.rds")), 
          title = "Distance to Closest \n Upstream Lake (standardized)")
# inv_inv_closest(-1.469, nes_rf_nws) # 2681.288

plot_tree(gettree(readRDS("../data/iws/sr_forest.rds")), 
          title = "Stream Order Ratio")
plot_tree(gettree(readRDS("../data/nws/sr_forest.rds")), 
          title = "Stream Order Ratio")

plot_tree(gettree(readRDS("../data/iws/linklength_forest.rds")), 
          title = "Link Length")
plot_tree(gettree(readRDS("../data/nws/linklength_forest.rds")), 
          title = "Link Length")

plot_tree(gettree(readRDS("../data/iws/uplakearea_forest.rds")), 
          title = "Upstream Lake Area")
# plot_tree(gettree(readRDS("01_Chapter/data/nws/uplakearea_forest.rds")), 
#           title = "Upstream Lake Area")

plot_tree(gettree(readRDS("../data/iws/baseflow_forest.rds")), 
          title = "Baseflow")

plot_tree(gettree(readRDS("../data/iws/wetland-cover_forest.rds")), 
          title = "Wetland Cover")

plot_tree(gettree(readRDS("../data/iws/streamdensity_forest.rds")), 
          title = "Stream Density")

# ---- 05_partition_vollenweider ----

plot_grid(
  plot_grid(
    part_pred_plot(nes_iws, readRDS("../data/lc_vollenweider.rds"), 
                   9, "A. Effect of lake connection", xl = TRUE, rev_legend = TRUE),
    part_pred_plot(nes_nws, readRDS("../data/nws/ll_vollenweider.rds"),
                   8, "B. Effect of link length \n (NWS)", xl = TRUE, yl = FALSE), 
    rel_widths = c(1, 0.9), 
    rel_heights = c(0.6, 1)),
  NULL,
  plot_grid(
    part_pred_plot(nes_iws, readRDS("../data/md_vollenweider.rds"), 
                 8, "C. Effect of max depth", xl = TRUE, legend_title = "depth"),
    part_pred_plot(nes_nws, readRDS("../data/nws/cd_vollenweider.rds"),
                 9, "D. Effect of distance to \n closest lake (NWS)", xl = TRUE, yl = FALSE), 
    rel_widths = c(1, 0.9), 
    rel_heights = c(0.6, 1)), 
  rel_heights = c(1, 0.05, 1), ncol = 1)

# ---- graphical_exploratory_analysis ----

partition_splits <- read.csv("../figures/table_1.csv", stringsAsFactors = FALSE)

nes_iws$p_pnt_source <-  rowSums(cbind(nes_iws$p_pnt_source_muni,
                                       nes_iws$p_pnt_source_septic,
                                       nes_iws$p_pnt_source_industrial),
                                 na.rm = TRUE)

md_v_la <- ggplot(data = nes_iws) + 
  geom_point(aes(x = maxdepth, 
                 y = upstream_lakes_4ha_area_ha, color = lakeconnection)) +
  ylab("Upstream Lake Area (ha)") + xlab("Max Depth (m)") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$splits[
      partition_splits$pnames2 == "md" & partition_splits$scale == "focal"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$splits[
      partition_splits$pnames2 == "la" & partition_splits$scale == "nws"])) + 
  scale_y_log10()

md_v_ll_iws <- ggplot(data = nes_iws) + 
  geom_point(aes(x = maxdepth, 
                 y = link_length, color = lakeconnection)) +
  ylab("Link Length") + xlab("Max Depth (m)") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$splits[
      partition_splits$pnames2 == "md"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$splits[
      partition_splits$pnames2 == "ll" & partition_splits$scale == "iws"]))

md_v_ll_nws <- ggplot(data = nes_nws) + 
  geom_point(aes(x = maxdepth, 
                 y = avg_link_length, color = lakeconnection)) +
  ylab("Link Length") + xlab("Max Depth (m)") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$splits[
      partition_splits$pnames2 == "md"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$splits[
      partition_splits$pnames2 == "ll" & partition_splits$scale == "nws"]))

md_v_cd_nws <- ggplot(data = nes_nws) + 
  geom_point(aes(x = maxdepth, 
                 y = closest_lake_distance, color = lakeconnection)) +
  ylab("Closest Lake Distance") + xlab("Max Depth (m)") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$splits[
      partition_splits$pnames2 == "md"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$splits[
      partition_splits$pnames2 == "cd" & partition_splits$scale == "nws"]))

ll_v_cd_nws <- ggplot(data = nes_nws) + 
  geom_point(aes(x = avg_link_length, 
                 y = closest_lake_distance, color = lakeconnection)) +
  ylab("Closest Lake Distance") + xlab("Link Length") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$splits[
      partition_splits$pnames2 == "ll" & partition_splits$scale == "nws"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$splits[
      partition_splits$pnames2 == "cd" & partition_splits$scale == "nws"])) + 
  scale_x_log10() + scale_y_log10() + xlab("log(Link Length)") + ylab("log(cd)")

legend <- get_legend(md_v_la)

plot_grid(md_v_la + theme(legend.position = "none"), 
          md_v_ll_iws + theme(legend.position = "none"), 
          md_v_ll_nws + theme(legend.position = "none"),
          ll_v_cd_nws + theme(legend.position = "none"),
          ggplot() + geom_point(aes(nes_iws$hu12_baseflowindex_mean, nes_nws$baseflow, color = nes_iws$state %in% c("MAINE", "NEW HAMPSHIRE", "IOWA", "MISSOURI", "OHIO", "ILLINOIS"))) + 
            geom_hline(aes(yintercept = 53.43)) +
            geom_vline(aes(xintercept = 63.76)) + 
            ylab("NWS Baseflow") + xlab("IWS Baseflow") + 
            theme(legend.position = "none"),
          ggplot() + geom_point(aes(
            nes_iws$iws_streamdensity_streams_density_mperha, 
            nes_nws$stream_density, color = nes_iws$lakeconnection)) + 
            geom_hline(aes(yintercept = 10.4), color = "red") +
            geom_vline(aes(xintercept = 4.43), color = "red") +
            geom_abline(aes(slope = 1, intercept = 0)) + 
            ylab("NWS Stream Density") + xlab("IWS Stream Density") + 
            theme(legend.position = "none"),
          
          ggplot() + geom_point(aes(nes_iws$closest_lake_distance, 
                                    nes_nws$closest_lake_distance, 
                                    color = nes_iws$lakeconnection)) + 
            geom_hline(aes(yintercept = 2776.81), color = "red") +
            geom_vline(aes(xintercept = 3773.61), color = "red") +
            geom_abline(aes(slope = 1, intercept = 0)) +
            scale_x_log10() + 
            scale_y_log10() + 
            ylab("NWS Lake Distance") + xlab("IWS Lake Distance") + 
            theme(legend.position = "none"),
          ggplot(data = nes_iws) + 
            geom_point(aes(x = p_pnt_source / (p_nonpnt_source + p_pnt_source), 
                           y = p_percent_retention)), 
          ncol = 2)

# ---- graphical_exploratory_analysis_cont ----

lg <- lagosne_load("1.087.1")
nes_nws_temp <- dplyr::left_join(nes_nws, select(lg$iws, lagoslakeid, iws_ha))
nes_nws_temp$p_pnt_source <-  rowSums(cbind(nes_nws_temp$p_pnt_source_muni,
                                       nes_nws_temp$p_pnt_source_septic,
                                       nes_nws_temp$p_pnt_source_industrial),
                                 na.rm = TRUE)
nes_nws_sf    <- coordinatize(nes_nws_temp, "lat", "long")
us_states <- st_intersects(us_states(), nes_nws_sf)
us_states <- us_states()[unlist(lapply(us_states, function(x) length(x) > 0)),]

plot_grid(
ggplot() + geom_point(data = nes_nws_temp, 
                      aes(x = iws_ha, y = avg_link_length)) + 
  xlim(300, 100000) + ylab("Average Link Length") + 
  xlab("Interlake Watershed Area (ha)") + 
  geom_hline(aes(yintercept = 2380.09), color = "red"),

ggplot() + geom_point(data = nes_nws_temp, 
                      aes(x = iws_ha, y = closest_lake_distance)) + 
  xlim(300, 100000) + ylab("Closest Lake Distance") + 
  xlab("Interlake Watershed Area (ha)") + 
  geom_hline(aes(yintercept = 3273.65), color = "red"), 

ggplot() + geom_point(data = nes_nws_temp, 
                      aes(x = iws_ha, y = p_percent_retention)) + 
  xlim(300, 100000) + ylab("P Retention") + 
  xlab("Interlake Watershed Area (ha)"),

ggplot() + geom_sf(data = us_states) + 
  geom_sf(data = nes_nws_sf, aes(color = p_pnt_source / (p_nonpnt_source + p_pnt_source))) + theme(legend.direction = "horizontal", legend.position = "bottom") + ggtitle("% pnt source P"),

ggplot() + geom_sf(data = us_states) + 
  geom_sf(data = nes_nws_sf, aes(color = log(maxdepth))) + theme(legend.direction = "horizontal", legend.position = "bottom") + ggtitle("max depth"), 
ggplot() + geom_sf(data = us_states) + 
  geom_sf(data = nes_nws_sf, aes(color = log(iws_ha))) + theme(legend.direction = "horizontal", legend.position = "bottom") + ggtitle("iws area"), 
ncol = 2)

# ---- graphical_exploratory_analysis_cont_cont ----

lg <- lagosne_load("1.087.1")

nes_nws_temp <- nes_nws
nes_nws_temp$p_pnt_source <-  rowSums(cbind(nes_nws_temp$p_pnt_source_muni,
                                            nes_nws_temp$p_pnt_source_septic,
                                            nes_nws_temp$p_pnt_source_industrial),
                                      na.rm = TRUE)
nes_nws_sf    <- coordinatize(nes_nws_temp, "lat", "long")
us_states <- st_intersects(us_states(), nes_nws_sf)
us_states <- us_states()[unlist(lapply(us_states, function(x) length(x) > 0)),]

plot_grid(
  ggplot() + geom_point(data = nes_nws_temp, 
                        aes(x = nws_ha, y = upstream_lakes_4ha_area_ha)) + 
    xlim(300, 100000) + ylab("Up. Lk. Area") + 
    xlab("Network Watershed Area (ha)"),
  ggplot() + geom_point(data = nes_nws_temp, 
                        aes(x = nws_ha, y = stream_order_ratio)) + 
    xlim(300, 100000) + ylab("Stream Order Ratio") + 
    xlab("Network Watershed Area (ha)"), 
  
  ggplot() + geom_point(data = nes_nws_temp, 
                        aes(x = nws_ha, y = p_percent_retention)) + 
    xlim(300, 100000) + ylab("P Retention") + 
    xlab("Network Watershed Area (ha)"),
  
  ggplot() + geom_point(data = nes_nws_temp, 
                        aes(x = nws_ha, 
                            y = p_pnt_source / (p_nonpnt_source + p_pnt_source))) + 
    xlim(300, 100000) + ylab("% Pnt Source") + 
    xlab("Network Watershed Area (ha)"), 
  ncol = 2)

# ---- 07_cor_mat_hmap ----
# correlation matrix heatmap
# setwd("scripts")
library(corrr)

splits      <- read.csv("../figures/table_1.csv", stringsAsFactors = FALSE)

nes_iws$p_pnt_source <-  rowSums(cbind(nes_iws$p_pnt_source_muni,
                                       nes_iws$p_pnt_source_septic,
                                       nes_iws$p_pnt_source_industrial),
                                 na.rm = TRUE)

nes_iws$p_pnt_source_pct <- nes_iws$p_pnt_source / 
  (nes_iws$p_nonpnt_source + nes_iws$p_pnt_source)

lg      <- lagosne_load("1.087.1")
nes_iws <- dplyr::left_join(nes_iws, select(lg$iws, lagoslakeid, iws_ha))

nes_nws$p_pnt_source <-  rowSums(cbind(nes_nws$p_pnt_source_muni,
                                       nes_nws$p_pnt_source_septic,
                                       nes_nws$p_pnt_source_industrial),
                                 na.rm = TRUE)

nes_nws$p_pnt_source_pct <- nes_nws$p_pnt_source / 
  (nes_nws$p_nonpnt_source + nes_nws$p_pnt_source)
nes_nws                  <- dplyr::left_join(nes_nws, 
                                             select(lg$iws, 
                                                    lagoslakeid, iws_ha))

# format iws data
nes_iws_sub <- dplyr::select(nes_iws, -lakeconnection)
conny_cols  <- unique(splits$iws_names)[!is.na(unique(splits$iws_names))]
conny_cols  <- conny_cols[!(conny_cols %in% "lakeconnection")]

nes_iws_sub <- nes_iws_sub[,c(conny_cols, 
                              "lake_area_ha", "p_pnt_source_pct", "iws_ha")]

iws_key <- merge(data.frame(iws_names = names(nes_iws_sub)), 
                 splits[,c("iws_names", "abb")], sort = FALSE)
iws_key <- rbind(iws_key, 
                 data.frame(iws_names = c("lake_area_ha", 
                                          "p_pnt_source_pct", "iws_ha"), 
                            abb = c("Lake Area", "Point Source P", "LWS area")))
iws_key <- iws_key[!duplicated(iws_key),]
names(nes_iws_sub) <- iws_key$abb

# format nws data
nes_nws_sub <- dplyr::select(nes_nws, -lakeconnection)
conny_cols  <- unique(splits$nws_names)[!is.na(unique(splits$nws_names))]
conny_cols  <- conny_cols[!(conny_cols %in% "lakeconnection")]


nes_nws_sub <- nes_nws_sub[,c(conny_cols, 
                              "lake_area_ha", "p_total", 
                              "p_pnt_source_pct", "iws_ha", 
                              "nws_ha", "retention_time_yr")]
nws_key <- merge(data.frame(nws_names = names(nes_nws_sub)), 
                 splits[,c("nws_names", "abb")], sort = FALSE)
nws_key <- rbind(nws_key, data.frame(nws_names = c("lake_area_ha", "p_total", 
                                                   "p_pnt_source_pct", "iws_ha", 
                                                   "nws_ha", "retention_time_yr"), 
                                     abb = c("Lake Area", "P Loading", 
                                             "Point Source P", "LWS area", 
                                             "NWS area", "Water Residence Time")))
nws_key <- nws_key[!duplicated(nws_key),]
names(nes_nws_sub) <- nws_key$abb

# combine iws and nws data
nes_sub <- dplyr::bind_rows(nes_iws_sub, nes_nws_sub)
# nes_sub <- nes_iws_sub
# nes_sub <- nes_nws_sub

# p_values
get_corrr_p_values <- function(d){
  var_pairs <- t(combn(names(d), 2)) %>%
    as_data_frame() %>% 
    setNames(c("x", "y"))
  
  p_values <- var_pairs %>% 
    dplyr::mutate(r.test = purrr::map2(x, y, ~ stats::cor.test(d[[.x]], d[[.y]])),
                  r.test = purrr::map(r.test, broom::tidy)) %>%
    tidyr::unnest(r.test)
  
  res <- matrix(NA, 
                nrow = length(names(d)), 
                ncol = length(names(d)))
  res[lower.tri(res)] <- p_values$p.value
  res <- data.frame(res, stringsAsFactors = FALSE, row.names = names(d))
  setNames(res, names(d))
}

p_values <- get_corrr_p_values(nes_sub)

# correlation matrix 
res <- shave(correlate(nes_sub))
sig_positions <- which(p_values > 0.05, TRUE)
for(i in seq_len(nrow(sig_positions))){
  res[sig_positions[i, 1], sig_positions[i, 2] + 1] <- NA
}

res <- res[apply(res[,2:ncol(res)], 1, function(x) !all(is.na(x))),]
res <- res[,c(TRUE, apply(res[,2:ncol(res)], 2, function(x) !all(is.na(x))))]
names(res)[1] <- ""

res_names_c      <- names(res)
res_f            <- data.frame(res)
res_names_r      <- res_f$Var.1
res_f            <- res_f[,-1]
names(res_f)     <- res_names_c[-1]
row.names(res_f) <- res_names_r
# res_f            <- res_f[rev(1:nrow(res_f)),]

# names(res) <- 
#   row.names(res) <- 
#   sapply(row.names(res), function(x) paste(strwrap(x, 10), collapse="\n "))

fnc <- function(var, decimal.places) {
  var <- sprintf(paste0("%1.", decimal.places, "f"), var)
  var[var=="NA"] <- ""
  var <- gsub("-", "", var)
  var
}

pheatmap::pheatmap(res_f, na_col = "grey", 
                   cluster_cols = FALSE, cluster_rows = FALSE, 
                   color = colorRampPalette(
                     rev(
                       RColorBrewer::brewer.pal(n = 7, name ="RdBu")[2:6]
                      ))(100), 
                   display_numbers = apply(res_f, 2, function(x) fnc(x, 2)), 
                   fontsize_number = 6)

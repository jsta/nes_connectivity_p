# setwd("scripts")
# source("99_utils.R")
# source("01_prepdata.R")

# ---- global_vollenweider_viz ----

full_pred_plot(nes_iws, readRDS("../data/global_vollenweider.rds"))

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

# ---- partition_vollenweider_viz ----

plot_grid(
  plot_grid(
    part_pred_plot(nes_iws, readRDS("../data/lc_vollenweider.rds"), 
                   9, "Lake connection", xl = FALSE),
    part_pred_plot(nes_nws, readRDS("../data/nws/ll_vollenweider.rds"),
                   8, "Link length", xl = FALSE, yl = FALSE), 
    rel_widths = c(1, 0.93), labels = c("IWS", "NWS"), 
    hjust = c(-4, -2.7), vjust = 2.5),
  NULL,
  plot_grid(
    part_pred_plot(nes_iws, readRDS("../data/md_vollenweider.rds"), 
                 8, "Max depth", xl = TRUE),
    part_pred_plot(nes_nws, readRDS("../data/nws/cd_vollenweider.rds"),
                 9, "Distance to \n closest lake", xl = TRUE, yl = FALSE), 
    rel_widths = c(1, 0.93)), 
  rel_heights = c(0.91, 0.1, 1), ncol = 1)

# ---- k_viz ----

fit_df_iws <- rbind(
  extract_coefs(readRDS("../data/global_vollenweider.rds"), "k"), 
  extract_coefs2(readRDS("../data/lc_vollenweider.rds"), 
                 c("k1_lc", "k2_lc")),
  extract_coefs2(readRDS("../data/md_vollenweider.rds"), 
                 c("k1_md", "k2_md")),
  extract_coefs2(readRDS("../data/iws/cd_vollenweider.rds"),
                 c("k1_cd", "k2_cd")),
  extract_coefs2(readRDS("../data/iws/ll_vollenweider.rds"),
                 c("k1_ll", "k2_ll")),
  extract_coefs2(readRDS("../data/iws/bf_vollenweider.rds"),
                 c("k1_bf", "k2_bf")),
  extract_coefs2(readRDS("../data/iws/sd_vollenweider.rds"),
                 c("k1_sd", "k2_sd")),
  extract_coefs2(readRDS("../data/iws/sr_vollenweider.rds"), 
                 c("k1_sr", "k2_sr")))
  #               extract_coefs2(fit_wc, c("k1_wc", "k2_wc")), 

# controls ordering of ridge plot
fit_df_iws$key <- factor(fit_df_iws$key, levels = c("k", 
                                            "k1_sd", "k2_sd",
                                            "k1_cd", "k2_cd",
                                            "k1_sr", "k2_sr",
                                            "k1_bf", "k2_bf",
                                            "k1_ll", "k2_ll",
                                            "k1_md", "k2_md",
                                            "k1_lc", "k2_lc"))
#                                             "k1_wc", "k2_wc",

k_viz_iws <- ggplot(
  fit_df_iws, 
  aes(x = coef, y = key, fill = key)) + 
    geom_density_ridges(size = 1, rel_min_height = 0.03) + 
    scale_fill_manual(values = c("grey", 
                            rep(c(viridis::viridis(1, begin = 0.5), 
                                  viridis::viridis(1, begin = 0)), 7))) + 
    cowplot::theme_cowplot() + theme(axis.text = element_text(size = 10),
                                     axis.title = element_text(size = 12),
                               plot.title = element_text(size = 12, face = "plain"),
                                     axis.line.x = element_line(size = 1, 
                                     colour = "black"),
                                     axis.line.y = element_line(size = 1, 
                                                       colour = "black"), 
                                     legend.position = "none", 
                                     plot.margin = margin(0.5,0,0,0, "cm")) + 
  xlab("Coefficient \n Value") + ylab("") + 
  scale_y_discrete(breaks = levels(fit_df_iws$key)[
    c(1, seq(3, length(levels(fit_df_iws$key)), by = 2))], 
                   labels = rev(c("Lake \n connection",  
                                                 "Max depth",
                                                 "Link length",
                                                 "Baseflow",
                                                 "Stream \n order ratio",
                                                 "Closest \n lake distance",
                                                 "Stream \n density",
                                                 "Global")))

fit_df_nws <- rbind(
  extract_coefs(readRDS("../data/global_vollenweider.rds"), "k"), 
  extract_coefs2(readRDS("../data/lc_vollenweider.rds"), 
                 c("k1_lc", "k2_lc")),
  extract_coefs2(readRDS("../data/nws/cd_vollenweider.rds"), 
                 c("k1_cd", "k2_cd")),
  extract_coefs2(readRDS("../data/nws/ll_vollenweider.rds"), 
                 c("k1_ll", "k2_ll")),
  extract_coefs2(readRDS("../data/nws/la_vollenweider.rds"), 
                 c("k1_la", "k2_la")),
  extract_coefs2(readRDS("../data/md_vollenweider.rds"), 
                 c("k1_md", "k2_md")),
  extract_coefs2(readRDS("../data/nws/bf_vollenweider.rds"), 
                 c("k1_bf", "k2_bf")),
  extract_coefs2(readRDS("../data/nws/sd_vollenweider.rds"), 
                 c("k1_sd", "k2_sd")),
  extract_coefs2(readRDS("../data/nws/sr_vollenweider.rds"), 
                 c("k1_sr", "k2_sr")))

fit_df_nws$key <- factor(fit_df_nws$key, levels = c("k",
                                                    "k1_sr", "k2_sr",
                                                    "k1_bf", "k2_bf",
                                                    "k1_md", "k2_md",
                                                    "k1_la", "k2_la",
                                                    "k1_lc", "k2_lc",
                                                    "k1_sd", "k2_sd",
                                                    "k1_cd", "k2_cd",
                                                    "k1_ll", "k2_ll"))

k_viz_nws <- ggplot(
  fit_df_nws, 
  aes(x = coef, y = key, fill = key)) + 
  geom_density_ridges(size = 1, rel_min_height = 0.03) + 
  scale_fill_manual(values = c("grey", 
                               rep(c(viridis::viridis(1, begin = 0.5), 
                                     viridis::viridis(1, begin = 0)), 8))) + 
  cowplot::theme_cowplot() + theme(axis.text = element_text(size = 10),
                                   axis.title = element_text(size = 12),
                               plot.title = element_text(size = 12, face = "plain"),
                                   axis.line.x = element_line(size = 1, 
                                                  colour = "black"),
                                   axis.line.y = element_line(size = 1, 
                                                  colour = "black"), 
                                   legend.position = "none", 
                                   plot.margin = margin(0.5,0,0,-0.7, "cm")) + 
  xlab("Coefficient \n Value") +  ylab("") +
  scale_y_discrete(breaks = levels(fit_df_nws$key)[
    c(1, seq(3, length(levels(fit_df_nws$key)), by = 2))], 
    labels = rev(c("Link length",
                   "Closest \n lake distance",
                   "Stream \n density",
                   "Lake \n connection",  
                   "Upstream \n lake area",
                   "Max depth",
                   "Baseflow",
                   "Stream \n order ratio",
                   "Global")))

plot_grid(
  k_viz_iws + ggtitle("IWS scale") + xlim(0.8, 1.6),
  k_viz_nws + ggtitle("NWS scale") + xlim(0.8, 1.6), ncol = 2,
  align = "v")

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
us_states <- st_intersects(us_states(), nes_sf)
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

# ---- maps ----

partition_splits <- read.csv("../figures/table_1.csv")
partition_splits <- dplyr::filter(partition_splits, !is.na(splits), 
                                  scale %in% c("nws", "iws"))
partition_splits <- droplevels(partition_splits)
partition_splits$pnames2 <- factor(partition_splits$pnames2, 
                              levels = c("cd", "ll", "la", "sd", "bf", "sr"))

nes_sf    <- coordinatize(nes_nws, "lat", "long")
us_states <- st_intersects(us_states(), nes_sf)
us_states <- us_states()[unlist(lapply(us_states, function(x) length(x) > 0)),]

get_sub <- function(dt, col_name, split_value){
  coordinatize(
    dplyr::filter(dt, UQ(rlang::sym(as.character(col_name))) <= split_value), 
    "lat", "long")
}

lower_iws_maps <- lapply(which(partition_splits$scale == "iws"), 
                         function(x) {
                           nes_sub <- get_sub(nes_iws, 
                                            partition_splits[x, "iws_names"],
                                            partition_splits[x, "splits"])
                           
                           ggplot() + 
                             geom_sf(data = us_states) +
                             geom_sf(data = nes_sf, 
                                     color = viridis::viridis(1, begin = 0),
                                     size = 0.6) + 
                             geom_sf(data = nes_sub, 
                                color = viridis::viridis(1, begin = 0.5), 
                                size = 0.6) + 
                             coord_sf(datum = NA) + 
                             ggtitle(partition_splits[x, "abb"]) + 
                             theme(plot.margin = unit(c(0, -0.1, 0, -0.1), "cm"),
                               plot.title = element_text(face = "plain", size = 12))
                         })

lower_nws_maps <- lapply(which(partition_splits$scale == "nws"), 
                         function(x) {
                           nes_sub <- get_sub(nes_nws, 
                                              partition_splits[x, "nws_names"],
                                              partition_splits[x, "splits"])
                           
                           nes_sf_temp <- nes_sf
                           st_geometry(nes_sf_temp) <- NULL
                           na_rows <- !is.na(
                    nes_sf_temp[,as.character(partition_splits[x, "nws_names"])])
                           
                           ggplot() + 
                             geom_sf(data = us_states) +
                             geom_sf(data = nes_sf[na_rows,], 
                                     color = viridis::viridis(1, begin = 0),
                                     size = 0.6) + 
                             geom_sf(data = nes_sub,
                                     color = viridis::viridis(1, begin = 0.5), 
                                     size = 0.6) + 
                             coord_sf(datum = NA) + 
                             ggtitle(partition_splits[x, "abb"]) + 
                             theme(plot.margin = unit(c(0, -0.1, 0, -0.1), "cm"), 
                               plot.title = element_text(face = "plain", size  = 12))
                         })

lower_iws_maps <- lower_iws_maps[
  match(levels(partition_splits$pnames2), 
        partition_splits$pnames2[which(partition_splits$scale == "iws")])
  ]

# lower_iws_maps[1][[1]] <- lower_iws_maps[1][[1]] + 
#   theme(plot.margin = unit(c(1,0,0,0), "lines"))

lower_nws_maps <- lower_nws_maps[
  match(levels(partition_splits$pnames2), 
        partition_splits$pnames2[which(partition_splits$scale == "nws")])
  ]

plot_grid(
  plot_grid(NULL,
            NULL, 
            NULL,
            NULL,
            NULL, 
            nrow = 4, ncol = 1,
            labels = c("IWS:", "", "NWS:", "", ""), 
            label_colour = "gray", vjust = c(1.2, 0, 2, 0, 0), 
            hjust = c(-0.1, 0, -0.12, 0, 0)),
  plot_grid(
    plot_grid(plotlist = lower_iws_maps, 
                   ncol = 3, nrow = 2, align = "v"),
    plot_grid(ggplot() + geom_blank()),
    plot_grid(plotlist = lower_nws_maps,
                   ncol = 3, nrow = 2, align = "v"), 
    ncol = 1, rel_heights = c(1, 0.1, 1)), 
rel_widths = c(0.2, 1), ncol = 2)

# ggplot() + geom_point(data = nes_nws, aes(x = lg_long, y = lg_lat, color = baseflow))

# ---- cor_mat_hmap ----
# correlation matrix heatmap

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

nes_nws$p_pnt_source <-  rowSums(cbind(nes_nws$p_pnt_source_muni,
                                       nes_nws$p_pnt_source_septic,
                                       nes_nws$p_pnt_source_industrial),
                                 na.rm = TRUE)

nes_nws$p_pnt_source_pct <- nes_nws$p_pnt_source / 
  (nes_nws$p_nonpnt_source + nes_nws$p_pnt_source)
nes_nws <- dplyr::left_join(nes_nws, select(lg$iws, lagoslakeid, iws_ha))

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

nes_nws_sub <- nes_nws_sub[,c(conny_cols, 
                              "lake_area_ha", "p_pnt_source_pct", "iws_ha")]
nws_key <- merge(data.frame(nws_names = names(nes_nws_sub)), 
                 splits[,c("nws_names", "abb")], sort = FALSE)
nws_key <- rbind(nws_key, data.frame(nws_names = c("lake_area_ha", "p_pnt_source_pct", "iws_ha"), abb = c("Lake Area", "Pnt. Src. P", "IWS area")))
nws_key <- nws_key[!duplicated(nws_key),]
names(nes_nws_sub) <- nws_key$abb

# combine iws and nws data
nes_sub <- dplyr::bind_rows(nes_iws_sub, nes_nws_sub)
# nes_sub <- nes_iws_sub
# nes_sub <- nes_nws_sub

# correlation matrix 
res <- shave(correlate(nes_sub))
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

pheatmap::pheatmap(res_f, na_col = "grey", 
                   cluster_cols = FALSE, cluster_rows = FALSE, 
                   color = colorRampPalette(
                     rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100))

# test <- superheat(
#           X = res_f, bottom.label.text.angle = 90, 
#           left.label.size = 1.4, bottom.label.size = 1.4, 
#           heat.pal = c("#542788", "white", "#b35806"), 
#           heat.pal.values = c(0, 0.6, 1))

# options(knitr.kable.NA = '')
# knitr::kable(res, digits = 2, 
#              caption = "Pearson correlation matrix among connectivity metrics.")

# # pander::panderOptions('keep.line.breaks', TRUE)
# pander::panderOptions("round", 2)
# pander::panderOptions("missing", '')
# # pander::pander(res, 
# #       caption = "Pearson correlation matrix among connectivity metrics", 
# #       split.cells = c(11, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8))
# pander::pandoc.table(res, style = "grid", keep.line.breaks = TRUE)

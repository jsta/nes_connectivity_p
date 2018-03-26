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
  part_pred_plot(nes_iws, readRDS("../data/lc_vollenweider.rds"), 
                 9, "Lake Connection", xl = FALSE),
  part_pred_plot(nes_iws, readRDS("../data/nws/cd_vollenweider.rds"),
                 9, "Distance to Closest Lake", xl = FALSE, yl = FALSE), 
  part_pred_plot(nes_iws, readRDS("../data/md_vollenweider.rds"), 
                 8, "Max Depth", xl = FALSE),
  part_pred_plot(nes_nws, readRDS("../data/nws/ll_vollenweider.rds"),
                 8, "Link Length", xl = FALSE, yl = FALSE),
  part_pred_plot(nes_iws, readRDS("../data/iws/ll_vollenweider.rds"),
                 7, "Link Length", xl = FALSE), 
  part_pred_plot(nes_iws, readRDS("../data/lc_vollenweider.rds"), 
                 7, "Lake Connection", yl = FALSE, xl = FALSE),
  part_pred_plot(nes_iws, readRDS("../data/iws/bf_vollenweider.rds"),
                 6, "Baseflow"),
  part_pred_plot(nes_nws, readRDS("../data/nws/la_vollenweider.rds"),
                 6, "Upstream Lake Area", yl = FALSE), 
  nrow = 4, ncol = 2, rel_widths = c(1, 0.9), 
  rel_heights = c(0.75, 0.75, 0.75, 1))

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
    geom_density_ridges(size = 2) + ylab("") + 
    scale_fill_manual(values = c(viridis::viridis(10)[1], 
                                 rep(viridis::viridis(10)[2:10], each = 2))) + 
    cowplot::theme_cowplot() + theme(text = element_text(size = 20), 
                                     axis.text = element_text(size = 10), 
                                     axis.line.x = element_line(size = 2, 
                                     colour = "black"),
                                     axis.line.y = element_line(size = 2, 
                                                       colour = "black"), 
                      legend.position = "none") + xlab("Coefficient \n Value")


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
                                                    "k1_bf", "k2_bf",
                                                    "k1_sr", "k2_sr",
                                                    "k1_sd", "k2_sd",
                                                    "k1_md", "k2_md",
                                                    "k1_la", "k2_la",
                                                    "k1_lc", "k2_lc",
                                                    "k1_ll", "k2_ll",
                                                    "k1_cd", "k2_cd"))

k_viz_nws <- ggplot(
  fit_df_nws, 
  aes(x = coef, y = key, fill = key)) + 
  geom_density_ridges(size = 2) + ylab("") + 
  scale_fill_manual(values = c(viridis::viridis(10)[1], 
                               rep(viridis::viridis(10)[2:10], each = 2))) + 
  cowplot::theme_cowplot() + theme(text = element_text(size = 20), 
                                   axis.text = element_text(size = 10), 
                                   axis.line.x = element_line(size = 2, 
                                                  colour = "black"),
                                   axis.line.y = element_line(size = 2, 
                                                  colour = "black"), 
                                   legend.position = "none") + 
  xlab("Coefficient \n Value")

plot_grid(
  k_viz_iws + ggtitle("IWS Scale") + xlim(0.75, 1.75),
  k_viz_nws + ggtitle("NWS Scale") + xlim(0.75, 1.75))

# ---- graphical_exploratory_analysis ----

partition_splits <- read.csv("../figures/table_1.csv", stringsAsFactors = FALSE)

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
          
          legend, ncol = 2)

# nes_iws$p_pnt_source <-  rowSums(cbind(nes_iws$p_pnt_source_muni,
#                                        nes_iws$p_pnt_source_septic,
#                                        nes_iws$p_pnt_source_industrial), 
#                                  na.rm = TRUE)
# 
# plot(nes_iws$p_pnt_source / 
#        nes_iws$p_nonpnt_source, nes_iws$p_percent_retention)

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
                             geom_sf(data = nes_sf, size = 0.6) + 
                             geom_sf(data = nes_sub, color = "red", size = 0.6) + 
                             coord_sf(datum = NA) + 
                             ggtitle(partition_splits[x, "pnames2"]) + 
                             theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
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
                                 size = 0.6) + 
                             geom_sf(data = nes_sub, color = "red", size = 0.6) + 
                             coord_sf(datum = NA) + 
                             ggtitle(partition_splits[x, "pnames2"]) + 
                             theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
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
  NULL,  plot_grid(plotlist = lower_iws_maps),
  NULL,  plot_grid(plotlist = lower_nws_maps), 
  ncol = 1, labels = c(NA, "iws", NA, "nws"), 
  rel_heights = c(0.2, 1, 0.2, 1), vjust = -0.6, label_colour = "gray")

# ggplot() + geom_point(data = nes_nws, aes(x = lg_long, y = lg_lat, color = baseflow))
source("99_utils.R")
source("01_prepdata.R")

# ---- load_packages ----

library(partykit)
library(data.tree)

library(forcats)
library(rstan)
library(nesRdata)
library(bayesplot)
library(patchwork)
library(tidyr)
library(cowplot)
library(ggsn)
library(tidybayes)
library(akima)

# ---- global_vollenweider ----

# fit <- full_model(nes_iws)
# saveRDS(fit, "01_Chapter/data/global_vollenweider.rds")

# ---- recursive_partioning ----

# max depth
# tree3 <- cforest(nes_rf$p_perce ~ ., data = dplyr::select(nes_rf,
#                                                           maxdept),
#                control = ctree_control(mincriterion = 0.8,
#                                        minsplit = 2))
# saveRDS(tree3, "01_Chapter/data/maxdepth_forest.rds")

# iws closest lake distance
# nes_rf_iws$inv_closest <- 1 / scale(nes_rf_iws$closest)
# tree3 <- cforest(nes_rf$p_perce ~ ., data = dplyr::select(nes_rf,
#                                                 inv_closest),
#                control = ctree_control(mincriterion = 0.55,
#                                        minsplit = 3))
# saveRDS(tree3, "01_Chapter/data/iws/cd_forest.rds")
# 
# nes_rf_nws$inv_closest <- 1 / scale(nes_rf_nws$closest)
# tree3 <- cforest(nes_rf_nws$p_perce ~ ., data = dplyr::select(nes_rf_nws,
#                                                 inv_closest),
#                control = ctree_control(mincriterion = 0.45,
#                                        minsplit = 4))
# # plot(gettree(tree3))
# saveRDS(tree3, "01_Chapter/data/nws/cd_forest.rds")
# inv_inv_closest(-1.51, nes_rf_iws) # 3773.614
#inv_inv_closest(-1.469, nes_rf_nws) # 2681.288

# stream order ratio
# tree3 <- cforest(nes_rf$p_perce ~ ., 
#         data = dplyr::select(nes_rf, 
#                              stream_), 
#         control = ctree_control(mincriterion = 0.1))
# saveRDS(tree3, "01_Chapter/data/iws/sr_forest.rds")
#
# tree3 <- cforest(nes_rf_nws$p_perce ~ .,
#         data = dplyr::select(nes_rf_nws,
#                              stream_),
#         control = ctree_control(mincriterion = 0.45))
# # plot(gettree(tree3))
# saveRDS(tree3, "01_Chapter/data/nws/sr_forest.rds")

# average link length
# tree3 <- cforest(nes_rf$p_perce ~ .,
#                  data = dplyr::select(nes_rf,
#                                       link_le),
#                  control = ctree_control(mincriterion = 0.3))
# saveRDS(tree3, "01_Chapter/data/iws/linklength_forest.rds")
#
# tree3 <- cforest(nes_rf_nws$p_perce ~ .,
#                  data = dplyr::select(nes_rf_nws,
#                                       link_le),
#                  control = ctree_control(mincriterion = 0.29))
# plot(gettree(tree3))
# saveRDS(tree3, "01_Chapter/data/nws/linklength_forest.rds")

# upstream lake area
# tree3 <- cforest(nes_rf$p_perce ~ .,
#                  data = dplyr::select(nes_rf_iws,
#                                                 upstrea),
#                  control = ctree_control(mincriterion = 0.7))
# saveRDS(tree3, "01_Chapter/data/iws/uplakearea_forest.rds")
# 
# tree3 <- cforest(nes_rf_nws$p_perce ~ .,
#                  data = dplyr::select(nes_rf_nws,
#                                                 upstrea),
#                  control = ctree_control(mincriterion = 0.7))
# plot(gettree(tree3))
# saveRDS(tree3, "01_Chapter/data/nws/uplakearea_forest.rds")

# Baseflow
# tree3 <- cforest(nes_rf$p_perce ~ .,
#                  data = dplyr::select(nes_rf,
#                                       hu12_ba),
#                  control = ctree_control(mincriterion = 0.2))
# saveRDS(tree3, "01_Chapter/baseflow_forest.rds")

# wetland cover
# tree3 <- cforest(nes_rf$p_perce ~ .,
#                  data = dplyr::select(nes_rf,
#                                       iws_wl_),
#                  control = ctree_control(mincriterion = 0.68))
# saveRDS(tree3, "01_Chapter/wetland-cover_forest.rds")

# stream density
# tree3 <- cforest(nes_rf$p_perce ~ .,
#                  data = dplyr::select(nes_rf,
#                                       iws_str),
#                  control = ctree_control(mincriterion = 0.4))
# saveRDS(tree3, "01_Chapter/streamdensity_forest.rds")

# ---- model_on_partitions ----

# Global model
# nes_sub <- make_standata(
#   formula = bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
#                k + x ~ 1, nl = TRUE),
#   data    = dplyr::select(nes, 
#                           p_percent_retention, retention_time_yr),
#   prior   = c(prior(normal(1.3, 0.3),  lb = 0, nlpar = k), 
#               prior(normal(0.45, 0.3), lb = 0, nlpar = x)))
# stan_mod <- stan_model(file = "01_Chapter/vollenweider.stan")
# m        <- sampling(stan_mod, nes_sub)
# saveRDS(m, "01_Chapter/global_vollenweider.rds")

m <- readRDS("../data/global_vollenweider.rds")

# Fit on partitions

# Max Depth
# fit <- part_model(nes_iws, "maxdepth", 0, 19.8, Inf, 1)
# saveRDS(fit, "../data/md_vollenweider.rds")

# lake connnection
# fit <- part_model(nes_iws, fac = "lakeconnection")
# saveRDS(fit, "../data/lc_vollenweider.rds")

# Closest lake distance
# fit <- part_model(nes_iws, "closest_lake_distance", 0, 3773.614, Inf, 2)
# saveRDS(fit, "../data/iws/cd_vollenweider.rds")
# fit <- part_model(nes_nws, "closest_lake_distance", 0, 2681.288, Inf, 2)
# saveRDS(fit, "../data/nws/cd_vollenweider.rds")

# Stream order ratio
# fit <- part_model(nes_iws, "stream_order_ratio", 0, 0.7, Inf)
# saveRDS(fit, "../data/iws/sr_vollenweider.rds")
# fit <- part_model(nes_nws, "stream_order_ratio", 0, 0.5, Inf)
# saveRDS(fit, "../data/nws/sr_vollenweider.rds")

# Link length
# fit <- part_model(nes_iws, "link_length", 0, 2177.1, Inf)
# saveRDS(fit, "../data/iws/ll_vollenweider.rds")
# fit <- part_model(nes_nws, "link_length", 0, 2237.34, Inf)
# saveRDS(fit, "../data/nws/ll_vollenweider.rds")

# not used currently 
# nes$part <- forcats::fct_recode(nes$part,
#                                 "1" = "2",
#                                 "2" = "1")

# Upstream lake area
# fit <- part_model(nes_nws, "upstream_lakes_4ha_area_ha", 0, 153.5, Inf, na_var = 1)
# saveRDS(fit, "../data/nws/la_vollenweider.rds")

# Baseflow
nes$part <- factor(as.numeric(cut(nes$hu12_baseflowindex_mean, 
                                  c(0, 63.8, Inf), include.highest = TRUE)))
# fit <- part_model(nes)
# saveRDS(fit, "01_Chapter/bf_vollenweider.rds")

fit_bf <- readRDS("01_Chapter/bf_vollenweider.rds")

gg_bf <- part_pred_plot(nes, fit_bf)

# Wetland cover

nes$part <- factor(as.numeric(cut(
  nes$iws_wl_allwetlandsdissolved_overlapping_area_pct, 
                                  c(0, 11.7, Inf), include.highest = TRUE)))
# fit <- part_model(nes)
# saveRDS(fit, "01_Chapter/wc_vollenweider.rds")

fit_wc <- readRDS("01_Chapter/wc_vollenweider.rds")

gg_wc <- part_pred_plot(nes, fit_wc)

# Stream density
# fit <- part_model(nes_iws, "iws_streamdensity_streams_density_mperha", 0, 4.4, Inf)
# saveRDS(fit, "../data/iws/sd_vollenweider.rds")

fit_df <- rbind(extract_coefs(m, "k"), 
                extract_coefs2(fit_md, c("k1_md", "k2_md")), 
                extract_coefs2(fit_cd, c("k1_cd", "k2_cd")), 
                extract_coefs2(fit_sr, c("k1_sr", "k2_sr")),
                extract_coefs2(fit_sd, c("k1_sd", "k2_sd")),
                extract_coefs2(fit_ll, c("k1_ll", "k2_ll")), 
                extract_coefs2(fit_la, c("k1_la", "k2_la")), 
                extract_coefs2(fit_bf, c("k1_bf", "k2_bf")), 
                extract_coefs2(fit_wc, c("k1_wc", "k2_wc")), 
                extract_coefs2(fit_lc, c("k1_lc", "k2_lc")))

fit_df$key <- factor(fit_df$key, levels = c("k", 
                                            "k1_wc", "k2_wc",
                                            "k1_sd", "k2_sd",
                                            "k1_cd", "k2_cd",
                                            "k2_sr", "k1_sr",
                                            "k1_bf", "k2_bf",
                                            "k2_ll", "k1_ll",
                                            "k1_la", "k2_la",
                                            "k1_lc", "k2_lc",
                                            "k1_md", "k2_md"))

# Compare coefficients

(gg_k <- ggplot(
  fit_df, 
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
                                     legend.position = "none")) + 
  xlab("Coefficient Value")

# ggplot(fit_df) + 
#   stat_boxplot(geom = 'errorbar', aes(x = coef,
#                                                    y = key,
#                                                    group = key))
# 
# ggplot(fit_df) + 
#   geom_boxplot(aes(x = key, y = coef))

# Compare predictions

gg_md_pred <- ggplot() +
  gg_md[1] + gg_md[2] + theme(legend.position = "na") + 
  geom_point(data = nes, aes(x = retention_time_yr, y = p_percent_retention)) + 
  scale_fill_viridis_d(begin = (1 / 9) * 9) + 
  ylab("P Retention (%)") +  xlab("Residence Time (yr)") + 
  ggtitle("Max Depth")

gg_md_pred + gg_md_pred + scale_x_log10() + ylab("")

gg_lc_pred <- ggplot() +
  gg_lc[1] + gg_lc[2] + theme(legend.position = "na") + 
  geom_point(data = nes, aes(x = retention_time_yr, y = p_percent_retention)) + 
  scale_fill_viridis_d(begin = (1 / 9) * 8) + 
  ylab("P Retention (%)") +  xlab("Residence Time (yr)") + 
  ggtitle("Lake Connection") 

gg_lc_pred + gg_lc_pred + scale_x_log10() + ylab("")

gg_la_pred <- ggplot() +
  gg_la[1] + gg_la[2] + theme(legend.position = "na") + 
  geom_point(data = nes, aes(x = retention_time_yr, y = p_percent_retention)) + 
  scale_fill_viridis_d(begin = (1 / 9) * 7) + 
  ylab("P Retention (%)") +  xlab("Residence Time (yr)") + 
  ggtitle("Upstream Lake Area") 

gg_la_pred + gg_la_pred + scale_x_log10() + ylab("")

gg_sr_pred <- ggplot() +
  gg_sr[1] + gg_sr[2] + theme(legend.position = "na") + 
  geom_point(data = nes, aes(x = retention_time_yr, y = p_percent_retention)) + 
  scale_fill_viridis_d(begin = (1 / 9) * 5) + 
  ylab("P Retention (%)") +  xlab("Residence Time (yr)") + 
  ggtitle("Stream Order Ratio") 

gg_sr_pred + gg_sr_pred + scale_x_log10() + ylab("")

gg_cd_pred <- ggplot() +
  gg_cd[1] + gg_cd[2] + theme(legend.position = "na") + 
  geom_point(data = nes, aes(x = retention_time_yr, y = p_percent_retention)) + 
  scale_fill_viridis_d(begin = (1 / 9) * 4) + 
  ylab("P Retention (%)") +  xlab("Residence Time (yr)") + 
  ggtitle("Distance to Closest \n Upstream Lake") 

gg_cd_pred + gg_cd_pred + scale_x_log10() + ylab("")

gg_wc_pred <- ggplot() +
  gg_wc[1] + gg_wc[2] + theme(legend.position = "na") + 
  geom_point(data = nes, aes(x = retention_time_yr, y = p_percent_retention)) + 
  scale_fill_viridis_d(begin = (1 / 9) * 2) + 
  ylab("P Retention (%)") +  xlab("Residence Time (yr)") + 
  ggtitle("Wetland Cover (%)")

gg_wc_pred + gg_wc_pred + scale_x_log10() + ylab("")

# ---- wind_roses ----

nes_sub <- dplyr::select(nes, 
                         stream_order_ratio:closest_lake_distance, 
                         hu12_baseflowindex_mean, 
                         iws_streamdensity_streams_density_mperha, 
                         upstream_lakes_4ha_area_ha, 
                         maxdepth, 
                         link_length, 
                         iws_wl_allwetlandsdissolved_overlapping_area_pct)

wind_rose_levels <- c("maxdepth", 
  "hu12_baseflowindex_mean",
  "upstream_lakes_4ha_area_ha", 
  "iws_wl_allwetlandsdissolved_overlapping_area_pct", 
  "stream_order_ratio", 
  "iws_streamdensity_streams_density_mperha", 
  "closest_lake_distance", 
  "link_length")

# global
nes_global <- nes_sub
nes_global$upstream_lakes_4ha_area_ha <- log(nes_global$upstream_lakes_4ha_area_ha)
nes_global <- dplyr::filter(nes_global, !is.na(closest_lake_distance))

test <- data.frame(apply(nes_global, 2, 
                         function(x) (x - mean(x, na.rm = TRUE)) / 
                           sd(x, na.rm = TRUE)))
test <- tidyr::gather(test, "key")
test$key <- factor(test$key, levels = wind_rose_levels)
test$quantile <- cut_number(test$value, 7)

gg <- ggplot(test, aes(x = key, fill = quantile)) 
gg_global <- gg + geom_bar(width = 1) + coord_polar() + 
  theme(axis.text.x = element_text(
angle = 360 / (2 * pi) * rev(pi / 2 + seq(pi / ncol(nes_sub), 2 * pi - pi / ncol(nes_sub), len = ncol(nes_sub))))) +  
  scale_fill_viridis_d(direction = 1, 
        labels = levels(cut(seq(0, 1, by = 0.1), 7, include.lowest = TRUE, 
                            dig.lab = 1))) + 
  theme(axis.title.x = element_blank(), 
        legend.box.margin = unit(c(0, 0, 0, 0), "cm"), 
        legend.direction = "horizontal", 
        legend.position = "bottom", 
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# partitioned

part_windrose_plot <- function(nes_sub, nes_column, partition, wind_rose_levels){

  if(!is.na(partition)){
    nes_sub$part <- factor(as.numeric(cut(nes_sub[,nes_column], 
                                    c(0, partition, Inf), include.highest = TRUE)))
    test1 <- nes_sub
    test1$upstream_lakes_4ha_area_ha <- log(test1$upstream_lakes_4ha_area_ha + 
                                              0.01)
    test1 <- dplyr::filter(test1, !is.na(closest_lake_distance))
  }else{ #lakeconnection
    nes_sub$part <- factor(as.numeric(nes_sub$lakeconnection))
    nes_sub <- dplyr::select(nes_sub, -lakeconnection)
    test1 <- nes_sub
    test1$upstream_lakes_4ha_area_ha <- log(test1$upstream_lakes_4ha_area_ha + 
                                              0.01)
    test1 <- rbind(test1[test1$part == 2,], 
                   test1[test1$part == 1 & !is.na(test1$closest_lake_distance),])
    # test1 <- test1[test1$part == 2 | 
    #                  (!is.na(test1$closest_lake_distance) & test1$part == 1),]
  }
  
  test1[,1:(ncol(test1) - 1)]  <- data.frame(
    apply(test1[,1:(ncol(test1) - 1)], 2, 
          function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
  
  test1 <- tidyr::gather(test1, "key", "value", -part)
  test1$key <- factor(test1$key, levels = wind_rose_levels)
  test1$quantile <- cut_number(test1$value, 7)
  
  gg   <- ggplot(dplyr::filter(test1, part == 1), aes(x = key, fill = quantile)) 
  gg_1 <- gg + geom_bar(width = 1) + coord_polar() + 
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    scale_fill_viridis_d(direction = 1)
  
  gg   <- ggplot(dplyr::filter(test1, part == 2), aes(x = key, fill = quantile)) 
  gg_2 <- gg + geom_bar(width = 1) + coord_polar() + 
    theme(legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
          scale_fill_viridis_d(direction = 1)
  
  list(gg_1, gg_2)
}

gg_md <- part_windrose_plot(nes_sub, "maxdepth", 20, wind_rose_levels)
gg_lc <- part_windrose_plot(dplyr::select(nes,
                                      stream_order_ratio:closest_lake_distance,
                                          hu12_baseflowindex_mean,
                                      iws_streamdensity_streams_density_mperha,
                                          upstream_lakes_4ha_area_ha,
                                          maxdepth,
                                          link_length,
                              iws_wl_allwetlandsdissolved_overlapping_area_pct,
                                        lakeconnection), "lakeconnection", NA, 
                            wind_rose_levels)
gg_la <- part_windrose_plot(nes_sub, "upstream_lakes_4ha_area_ha", 15584, 
                            wind_rose_levels)
gg_sr <- part_windrose_plot(nes_sub, "stream_order_ratio", 0.7, wind_rose_levels)
gg_cd <- part_windrose_plot(nes_sub, "closest_lake_distance", 9148.6, 
                            wind_rose_levels)
gg_ll <- part_windrose_plot(nes_sub, "link_length", 2463.9, wind_rose_levels)
gg_wc <- part_windrose_plot(nes_sub, "iws_wl_allwetlandsdissolved_overlapping_area_pct", 12.1, wind_rose_levels)
gg_bf <- part_windrose_plot(nes_sub, "hu12_baseflowindex_mean", 57.5, wind_rose_levels)

## plotting

gg_global

gg_md[1][[1]] + labs(subtitle = "Partition on max depth") +
  gg_md[2][[1]] + plot_layout(ncol = 2) +
gg_lc[1][[1]] + labs(subtitle = "Partition on lake connection") +
  gg_lc[2][[1]] + plot_layout(ncol = 2)
gg_la[1][[1]] + labs(subtitle = "Partition on upstream lake area") +
  gg_la[2][[1]] + plot_layout(ncol = 2) +
gg_sr[1][[1]] + labs(subtitle = "Partition on stream order ratio") +
  gg_sr[2][[1]] + plot_layout(ncol = 2)
gg_cd[1][[1]] + labs(subtitle = "Partition on closest lake distance") +
  gg_cd[2][[1]] + plot_layout(ncol = 2) +
gg_ll[1][[1]] + labs(subtitle = "Partition on link length") +
  gg_ll[2][[1]] + plot_layout(ncol = 2) 
gg_bf[1][[1]] + labs(subtitle = "Partition on baseflow") +
  gg_bf[2][[1]] + plot_layout(ncol = 2) +
gg_wc[1][[1]] + labs(subtitle = "Partition on wetland cover") +
  gg_wc[2][[1]] + plot_layout(ncol = 2)
   

## dont clip labels
# library(grid)
# 
# gt <- ggplot_gtable(ggplot_build(gg_2))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)

# ---- post_hoc_plots ----

md_v_la <- ggplot(data = nes) + 
  geom_point(aes(x = maxdepth, 
                 y = upstream_lakes_4ha_area_ha, color = lakeconnection)) +
  ylab("Upstream Lake Area (ha)") + xlab("Max Depth (m)") + 
    geom_vline(data = partition_splits, aes(
      xintercept = partition_splits$split_value[
                    partition_splits$col_names == "maxdepth"])) + 
    geom_hline(data = partition_splits, aes(
      yintercept = partition_splits$split_value[
              partition_splits$col_names == "upstream_lakes_4ha_area_ha"])) + 
  scale_y_log10()

md_v_ll <- ggplot(data = nes) + 
  geom_point(aes(x = maxdepth, 
                 y = link_length, color = lakeconnection)) +
  ylab("Link Length") + xlab("Max Depth (m)") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$split_value[
      partition_splits$col_names == "maxdepth"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$split_value[
      partition_splits$col_names == "link_length"]))

la_v_ll <- ggplot(data = nes) + 
  geom_point(aes(x = link_length, 
                 y = upstream_lakes_4ha_area_ha, color = lakeconnection)) +
  ylab("Upstream Lake Area (ha)") + xlab("Link Length") + 
  geom_vline(data = partition_splits, aes(
    xintercept = partition_splits$split_value[
      partition_splits$col_names == "link_length"])) + 
  geom_hline(data = partition_splits, aes(
    yintercept = partition_splits$split_value[
      partition_splits$col_names == "upstream_lakes_4ha_area_ha"])) + 
  scale_y_log10()

legend <- get_legend(md_v_la)

plot_grid(md_v_la + theme(legend.position = "none"), 
          md_v_ll + theme(legend.position = "none"), 
          la_v_ll + theme(legend.position = "none"), 
          legend, nrow = 2, rel_heights = c(1, 1))

# ---- connectivity_index ----

conn_index <- function(nes){
  res <- data.frame(conn_index   = matrix(NA, nrow = nrow(nes)),
                    lake_index   = matrix(NA, nrow = nrow(nes)),
                    stream_index = matrix(NA, nrow = nrow(nes)))
  res$la <- boot::inv.logit(scale(nes$upstream_lakes_4ha_area_ha))
  # res$la[res$la == min(res$la)] <- NA
  res$lc <- ((as.numeric(fct_rev(nes$lakeconnection)) - 1) / 3) + 0.33333
  res$ll <- boot::inv.logit(scale(nes$link_length))
  res$sr <- boot::inv.logit(scale(
    (1 + max(nes$stream_order_ratio, na.rm = TRUE)) - nes$stream_order_ratio))
  
  res$conn_index <- rowSums(res[,c(4, 6, 7)], na.rm = TRUE) /
    max(rowSums(res[,c(4, 6, 7)], na.rm = TRUE))
  res$lake_index <- rowSums(res[,4], na.rm = TRUE) /
    max(rowSums(res[,4], na.rm = TRUE))
  res$stream_index <- rowSums(res[,6:7], na.rm = TRUE) /
    max(rowSums(res[,6:7], na.rm = TRUE))
  res
}
 
nes$conn_index   <- conn_index(nes)[,1]
nes$lake_index   <- conn_index(nes)[,2]
nes$stream_index <- conn_index(nes)[,3]

lg  <- lagosne_load("1.087.1")
nes <- left_join(nes, dplyr::select(lg$iws, lagoslakeid, iws_ha))

# plot(nes$iws_ha, nes$conn_index)
# plot(nes$retention_time_yr, nes$conn_index)

nes_sub <- select(nes, conn_index, stream_index, lake_index, 
                  p_percent_retention, link_length, lat, retention_time_yr)

gg_ci_vs_p <- ggplot(nes_sub, aes(x = conn_index, y = p_percent_retention)) + 
  geom_point() + stat_smooth(method = lm, se = FALSE) + 
  ylab("P Retention (%)") + xlab("Connectivity Index")

gg_si_vs_p <- ggplot(nes_sub, aes(x = stream_index, y = p_percent_retention)) + 
    geom_point() + stat_smooth(method = lm, se = FALSE) +
    ylab("P Retention (%)") + xlab("Stream Connectivity Index")

gg_li_vs_p <- ggplot(nes_sub, aes(x = lake_index, y = p_percent_retention)) + 
    geom_point() + # add segmented regression results
    ylab("P Retention (%)") + xlab("Lake Connectivity Index")

gg_si_vs_li <- ggplot(nes_sub, aes(x = stream_index, y = lake_index)) + 
  geom_point() +
  ylab("Lake Connectivity Index") + xlab("Stream Connectivity Index")

gg_si_vs_lat <- ggplot(nes_sub, aes(x = stream_index, y = lat)) + 
  geom_point() +
  ylab("Lattitude") + xlab("Stream Connectivity Index")

gg_li_vs_lat <- ggplot(nes_sub, aes(x = lake_index, y = lat)) + 
  geom_point() +
  ylab("Lattitude") + xlab("Lake Connectivity Index")

gg_si_vs_ret <- ggplot(nes_sub, aes(x = stream_index, y = retention_time_yr)) + 
  geom_point() +
  ylab("Residence Time (yr)") + xlab("Stream Connectivity Index")

gg_li_vs_ret <- ggplot(nes_sub, aes(x = lake_index, y = retention_time_yr)) + 
    geom_point() +
    ylab("Residence Time (yr)") + xlab("Lake Connectivity Index")

gg_blank <- ggplot() + geom_blank()

plot_grid(gg_ci_vs_p, 
          gg_si_vs_p,
          gg_li_vs_p,
          gg_si_vs_li,
          gg_si_vs_lat,
          gg_li_vs_lat,
          geom_blank,
          gg_si_vs_ret,
          gg_li_vs_ret,
          nrow = 3)

# plot(nes$iws_ha, nes$lake_index)
# plot(nes$retention_time_yr, nes$lake_index)

# plot(nes$iws_ha, nes$stream_index)
# plot(nes$retention_time_yr, nes$stream_index)

# df <- tidyr::complete(dplyr::count(nes, lake_index, stream_index), lake_index, stream_index, fill = list(n = 0))

df <- dplyr::count(nes, lake_index, stream_index)

# resolution <- 0.02
# image(interp(x = df$stream_index, y = df$lake_index, z = df$n,
#             xo = seq(min(df$stream_index), max(df$stream_index),
#                      by = resolution),
#             yo = seq(min(df$lake_index), max(df$lake_index), by = resolution),  duplicate = "mean"))

# ggplot(df,
#        aes(stream_index, lake_index)) +
#   stat_density2d(aes(fill = ..density..), geom = "tile", contour = FALSE) +
#   scale_fill_gradient2(low = "red", high = "white")

nes$lake_type <- factor(nes$lake_type)

plot_grid(ggplot(nes, aes(x = lake_type, y = lake_index)) +
  geom_boxplot(outlier.shape = NA) + ylab("Lake Connectivity") + xlab("") + 
  scale_y_continuous(limits = c(0.4, 0.5)),
ggplot(nes, aes(x = lake_type, y = stream_index)) +
  geom_boxplot(outlier.shape = NA) + ylab("Stream Connectivity") + xlab(""))

#

# rbind(
#   head(dplyr::filter(nes[order(nes$conn_index),],
#                    iws_ha < 15000, retention_time_yr < 2,
#                    lake_area_ha < 1000, lake_area_ha > 200,
#                    p_percent_retention > 0.6))[1:3,],
# head(dplyr::filter(nes[order(nes$conn_index, decreasing = TRUE),],
#                    iws_ha < 10000, retention_time_yr < 2,
#                    lake_area_ha < 1000,
#                    p_percent_retention < 0.4))[1,])
# head(nes[order(nes$conn_index, decreasing = TRUE),])

# ggplot(data = nes) + 
#   geom_point(aes(x = log(retention_time_yr), 
#                  y = p_percent_retention, 
#                  color = conn_index))




# nes$lakeconnection      <- factor(tolower(substr(nes$lakeconnection, 4, 
#                                nchar(as.character(nes$lakeconnection)))), 
#                                levels = c("stream", "lakestream"))


# ggplot(data = nes) + 
#   geom_point(aes(x = log(retention_time_yr), 
#                  y = log(p_percent_retention), 
#                  color = lakeconnection))
# 
# ggplot(data = nes) + 
#   geom_point(aes(x = log(retention_time_yr), 
#                  y = log(p_percent_retention), 
#                  color = closest_lake_distance))
# 
# ggplot(data = nes) + 
#   geom_point(aes(x = log(retention_time_yr), 
#                  y = log(p_percent_retention), 
#                  color = stream_order_ratio))
# 
# 
# ggplot(data = nes) + 
#   geom_point(aes(x = log(retention_time_yr), 
#                  y = log(p_percent_retention), 
#                  color = log(upstream_lakes_4ha_area_ha)))
# 
# ggplot(nes, aes(x = log(retention_time_yr), y = lakeconnection)) +
#   geom_density_ridges()
# 
# ggplot(nes, aes(x = log(p_percent_retention), y = lakeconnection)) +
#   geom_density_ridges()
  
# ---- brms ----

# simple vollenweider
# fit <- brm(bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), k + x ~ 1, nl = TRUE), data = nes, 
#            prior = c(prior(normal(1.3, 0.1), lb = 0, nlpar = k), 
#                      prior(normal(0.45, 0.1), lb = 0, nlpar = x)), 
#            iter = 4000)

nes$p_percent_retention <- boot::logit(nes$p_percent_retention)
nes_sub <- make_standata(
  formula = bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
               k + x ~ 1, nl = TRUE),
  data    = dplyr::select(nes, 
                          p_percent_retention, retention_time_yr),
  prior   = c(prior(normal(1.3, 0.3),  lb = 0, nlpar = k), 
              prior(normal(0.45, 0.3), lb = 0, nlpar = x)))
stan_mod <- stan_model(file = "01_Chapter/vollenweider.stan")
m        <- sampling(stan_mod, nes_sub)

theta  <- extract(m)
ppc_dens_overlay(boot::inv.logit(nes$p_percent_retention), 
                 t(apply(theta$Y_pred[1:50,], 1, boot::inv.logit)))
test <- lapply(ncol(theta$Y_orig), function(x) boot::inv.logit(theta$Y_pred) / boot::inv.logit(theta$Y_orig))
hist(unlist(test), xlim = c(-2, 10))
abline(v = 1)
abline(v = median(unlist(test)), col = "red")

# https://stats.stackexchange.com/questions/86273/calculate-log-likelihood-by-hand-for-generalized-nonlinear-least-squares-regre



test2 <- median(extract(m)$sigma) * sqrt((nrow(nes) - 2)/ nrow(nes))
sum(dnorm(nes$p_percent_retention, 
          mean = theta$Y_pred[1,], 
          sd = test2, log = TRUE))

log_lik_1 <- extract_log_lik(m)
loo_1 <- loo(log_lik_1)

# group-level "random" effect of lake connection
# fit <- brm(bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), x ~ 1, k ~ 1 + (1|lakeconnection), nl = TRUE), data = nes, 
#            prior = c(prior(normal(1.3, 0.1), lb = 0, nlpar = k), 
#                      prior(normal(0.45, 0.1), lb = 0, nlpar = x)), 
#            iter = 4000)

# ---- hierarchical_on_lake-connection ----
# nes$p_percent_retention <- nes$p_percent_retention - 0.01
fit <- brm(bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
              x ~ 1, k ~ 0 + lakeconnection, nl = TRUE), 
           data = nes, 
           prior = c(prior(normal(1.3, 0.1),  lb = 0, nlpar = k), 
                     prior(normal(0.45, 0.1), lb = 0, nlpar = x)), 
           iter = 4000)

# # try Beta
# fit <- brm(bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
#               x ~ 1, k ~ 0 + lakeconnection, nl = TRUE), 
#            data = nes, 
#            prior = c(prior(gamma(8, 1), class = "phi"), 
#                      prior(normal(1.3, 0.1),  lb = 0, nlpar = k), 
#                      prior(normal(0.45, 0.1), lb = 0, nlpar = x)), 
#            iter = 4000, family = Beta())

# from external stan file
library(rstan)
library(tidybayes)
# nes_sub2 <- compose_data(
#   dplyr::select(nes, p_percent_retention, retention_time_yr, lakeconnection))

nes_sub <- make_standata(
  formula = bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
               x ~ 1, k ~ 0 + lakeconnection, nl = TRUE),
  data    = dplyr::select(nes, 
                       p_percent_retention, retention_time_yr, lakeconnection),
  prior   = c(prior(normal(1.3, 0.1),  lb = 0, nlpar = k), 
              prior(normal(0.45, 0.1), lb = 0, nlpar = x)))
stan_mod <- stan_model(file = "01_Chapter/hier_lakeconn.stan")
m <- sampling(stan_mod, nes_sub)

theta2  <- extract(m)
ppc_dens_overlay(nes$p_percent_retention, theta2$Y_pred[1:50,])
test2 <- lapply(ncol(theta2$Y_orig), function(x) theta2$Y_pred / theta2$Y_orig)

hist(unlist(test2), xlim = c(-1, 3))
abline(v = 1)
abline(v = median(unlist(test2)), col = "red")

#### 

nes_sub <- make_standata(
  formula = bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
               x ~ 1, k ~ 0 + maxdepth + link_length, nl = TRUE),
  data    = dplyr::select(nes, 
                 p_percent_retention, retention_time_yr, maxdepth, link_length),
  prior   = c(prior(normal(1.3, 0.1),  lb = 0, nlpar = k), 
              prior(normal(0.45, 0.1), lb = 0, nlpar = x)))
stan_mod <- stan_model(file = "01_Chapter/hier_lakeconn.stan")
m <- sampling(stan_mod, nes_sub)

theta3  <- extract(m)
ppc_dens_overlay(theta3$Y_orig[1,], theta3$Y_pred[1:50,])
test3 <- lapply(ncol(theta3$Y_orig), function(x) theta3$Y_pred / theta3$Y_orig)

library(ggridges)
test4 <- data.frame(theta_ratio = c(unlist(test), 
                                    unlist(test2), 
                                    unlist(test3)), 
            model = c(rep("base", length(unlist(test)), 
                      rep("hier", length(unlist(test2)), 
                      rep("hier2", length(unlist(test3)))))))

ggplot(test4, aes(x = theta_ratio, y = model)) + 
  geom_density_ridges() + 
  xlim(-1, 5) + geom_vline(aes(xintercept = 1))



# saveRDS(fit, "hier_lakeconn.rds")

# ---- hierarchical_on_lake-connection_upstream-lake-area ----
fit <- brm(bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), x ~ 1, k ~ lakeconnection + upstream_lakes_4ha_area_ha, nl = TRUE), data = nes, 
           prior = c(prior(normal(1.3, 0.1), lb = 0, nlpar = k), 
                     prior(normal(0.45, 0.1), lb = 0, nlpar = x)), 
           iter = 4000)

# hierarchical model of tp_lake on lake connection
# fit <- brm(bf(tp ~ 0 + lakeconnection), data = nes,
#            iter = 2000, chains = 1, family = lognormal())
# 
# fit <- brm(bf(tp ~ tp_in / (1 + k * retention_time_yr), k ~ 0 + lakeconnection, nl = TRUE), 
#            prior = c(prior(normal(1.3, 1e6), lb = 0, nlpar = k)), 
#            iter = 2000, chains = 1, data = nes, family = lognormal())

# hist(nes$p_percent_retention)
# hist(rlnorm(nrow(nes), ))

# ---- model_plotting ----
fit <- readRDS("hier_lakeconn.rds")

fit_df <- as.data.frame(fit)
names(fit_df) <- c("x", "k (lakestream)", "k (stream)", "sigma", "lp")
fit_df <- tidyr::gather(fit_df, value = "removal rate coeffcient (k)") %>% 
  filter(key %in% c("k (lakestream)", "k (stream)"))

lakeconnection2 <- data.frame(lakeconnection = c("DR_Stream", "DR_LakeStream"), 
                              lakeconnection2 = c("stream", "lakestream"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(3)[2:3]

nes <- left_join(nes, lakeconnection2, by = "lakeconnection")

(gg_ret_time <- ggplot(nes, 
  aes(x = log(retention_time_yr), y = lakeconnection2, 
      fill = lakeconnection2)) +
  geom_density_ridges(size = 2) + ylab("") + 
  scale_fill_cyclical(values = c(pal)) + 
  cowplot::theme_cowplot() + theme(text = element_text(size = 40), 
                                         axis.text = element_text(size = 30), 
                                         axis.line.x = element_line(size = 2, 
                                                          colour = "black"),
                                         axis.line.y = element_line(size = 2, 
                                                        colour = "black")))
gg_ret_p <- ggplot(nes, 
  aes(x = log(p_percent_retention), y = lakeconnection2, 
      fill = lakeconnection2)) +
  geom_density_ridges(size = 2) + ylab("") + 
  scale_fill_cyclical(values = c(pal)) + 
  cowplot::theme_cowplot() + theme(text = element_text(size = 40), 
                                     axis.text = element_text(size = 30), 
                                     axis.line.x = element_line(size = 2, 
                                                           colour = "black"),
                                     axis.line.y = element_line(size = 2, 
                                                            colour = "black"))

(gg_k <- ggplot(
  fit_df, 
  aes(x = `removal rate coeffcient (k)`, y = key, fill = key)) + 
  geom_density_ridges(size = 2) + ylab("") + 
    scale_fill_cyclical(values = c(pal)) + 
    cowplot::theme_cowplot() + theme(text = element_text(size = 40), 
                                     axis.text = element_text(size = 30), 
                                     axis.line.x = element_line(size = 2, 
                                                          colour = "black"),
                                     axis.line.y = element_line(size = 2, 
                                                            colour = "black")
                                     ))

cowplot::plot_grid(gg_ret_time, gg_ret_p, 
                   rel_widths = c(0.5, 0.5), ncol = 1)

gg_k

# ---- model_checking ----

plot(fit)
plot(marginal_effects(fit), points = TRUE)
test <- data.frame(lakeconnection = unique(nes$lakeconnection))
rownames(test) <- unique(nes$lakeconnection)
plot(marginal_effects(fit, conditions = test, 
                      re_formula = NULL, method = "predict"), points = TRUE)

test <- data.frame(a = rep(seq(1, 4, by = 1), each = 2))
test$b <- rep(c("g", "h"), length.out = nrow(test))

new_data <- data.frame(retention_time_yr = 
                         exp(rep(seq(
                           -10, 
                           log(max(nes$retention_time_yr)), by = 0.1), 
                             each = 2)))
new_data$lakeconnection <- rep(c("DR_LakeStream", "DR_Stream"), 
                               length.out = nrow(new_data))
new_data <- cbind(new_data, predict(fit, newdata = new_data))

plot(log(new_data$retention_time_yr), log(new_data$Estimate), pch = 19, 
       col = rep(c("red", "green"), length.out = nrow(new_data)))
points(log(nes$retention_time_yr), log(nes$p_percent_retention))

# See https://rpubs.com/bbolker/3423

# lme ####




# non-linear mle ####
# http://www.magesblog.com/2015/10/non-linear-growth-curves-with-stan.html

nes_stream <- dplyr::filter(nes, lakeconnection == "DR_Stream")
nes_lakestream <- dplyr::filter(nes, lakeconnection == "DR_LakeStream")

nlm <- nls(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
           data = nes, start = list(k = 0.1, x = 0.1))

nlm_resid <- nes$p_percent_retention - predict(nlm)

nlm_stream <- nls(
  p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
                  data = nes_stream, start = list(k = 0.1, x = 0.1))
  
nlm_lakestream <- nls(
  p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)), 
  data = nes_lakestream, start = list(k = 0.1, x = 0.1))


# plot(fitted(nlm), nes$p_percent_retention)
# hist(residuals(nlm))
# qqplot(residuals(nlm), nes$p_percent_retention)
# abline(0.4, 1)
# 
# confint(nlm)
# confint(nlm_stream)
# confint(nlm_lakestream)
# 
# hist(fitted(nlm))
# hist(nes$p_percent_retention)
# 
# plot(log(nes$retention_time_yr), nes$p_percent_retention)
# points(log(nes$retention_time_yr), fitted(nlm), col = "red", pch = 19)
# points(log(nes_stream$retention_time_yr), fitted(nlm_stream), col = "blue", pch = 19)
# points(log(nes_lakestream$retention_time_yr), fitted(nlm_lakestream), col = "darkgreen", pch = 19)

cor.test(nes$p_percent_retention, fitted(nlm)) # pseudo R^2

get_coefs <- function(md){
  res <- cbind(broom::tidy(confint(md)), broom::tidy(md))
  # res$modname <- deparse(substitute(md))
  res
}

dt <- list(nlm_stream, nlm_lakestream)
dt <- lapply(dt, get_coefs)
dt <- do.call("rbind", dt)
dt$modname <- factor(rep(c("stream", "lakestream"), each = 2), 
                     levels = c("stream", "lakestream"))

ggplot(dt) + geom_boxplot(aes(x = modname, 
                              ymin = X2.5., lower = estimate - std.error, 
                              middle = estimate, 
                              upper = estimate + std.error, ymax = X97.5.), 
                          stat = "identity") + 
  facet_wrap(~term, scales = "free") + 
  cowplot::theme_cowplot() + xlab("Lake Connection")
ggsave("01_Chapter/figures/nls_coefs.pdf")


gg_retention <- ggplot(nes) + 
  geom_boxplot(aes(x = lakeconnection, y = p_percent_retention)) + 
  cowplot::theme_cowplot() + ylab("P Retention (%)") + xlab("Lake Connection")

# time labels ####

yr_labels <- c("1e-07", "1e-05", 
               as.character(10/365), 
               as.character(2/12), as.character(6/12), 
               "2", "10")

format_labels <- function(x, mult_factor){
  gsub("\\.+$", "\\1", gsub("0+$", "\\1",
                            trimws(format(round(as.numeric(x), 6) * 
                            mult_factor, scientific = FALSE)), perl = TRUE))
}

day_labels   <- ceiling(as.numeric(
  format_labels(yr_labels, 365)) / 10) * 10
minute_labels  <- ceiling(as.numeric(
  format_labels(yr_labels, 365 * 24 * 60)) / 10) * 10
month_labels  <- round(as.numeric(
  format_labels(yr_labels, 12)), 0)

mixed_labels <- paste0(
  c(minute_labels[1:2], day_labels[3], month_labels[4:5], yr_labels[6:7]), 
  c(" min", " min", " days", " months", " months", " years", " years"))

gg_residence <- ggplot(nes) + 
  geom_boxplot(aes(x = lakeconnection, y = retention_time_yr), 
               outlier.shape = NA) + 
  scale_y_log10(labels = mixed_labels, breaks = as.numeric(yr_labels), 
                limits = c(1 / 365, 11)) + 
  ylab("Residence time log(yr)") + xlab("Lake Connection") + 
  cowplot::theme_cowplot()

cowplot::plot_grid(gg_retention, gg_residence)
ggsave("01_Chapter/figures/p_retention_boxplot.pdf")

# non-linear mixed effects
library(lme4)

nform <- ~ 1 - (1 / (1 + k * retention_time_yr ^ x))
nfun <- deriv(nform, namevec = c("k", "x"), 
              function.arg = c("retention_time_yr", "k", "x"))
fit <- nlmer(p_percent_retention ~ nfun(retention_time_yr, k, x) ~ 
               (k|lakeconnection), 
      start = c(k = 1.2, x = 0.3), data = nes)

nes$predict <- predict(fit)
# nlmer(p_percent_retention ~ nfun(retention_time_yr, k, x) ~ k|hu2, 
#       start = c(k = 1.2, x = 0.3), data = nes)

ggplot(data = nes) + 
  geom_point(aes(x = log(retention_time_yr), 
                 y = log(p_percent_retention), 
                 color = lakeconnection)) + 
  geom_point(aes(x = log(retention_time_yr), 
                 y = log(predict), 
                 color = lakeconnection))

# beta regression ####

library(lme4)
library(bbmle)
library(lattice)
library(ggplot2)

library(betareg)

mod0 <- betareg(p_percent_retention ~ log(retention_time_yr) + 
                  iws_nlcd1992_pct_81 + 
                  iws_nlcd1992_pct_82 + 
                  iws_streamdensity_streams_density_mperha +
                  hu12_baseflowindex_mean + 
                  scale(upstream_lakes_4ha_area_ha), data = nes)


library(glmmTMB)

mod1 <- glmmTMB(p_percent_retention ~ scale(retention_time_yr) + 
                  iws_nlcd1992_pct_81 + 
                  iws_nlcd1992_pct_82 + 
                  iws_streamdensity_streams_density_mperha +
                  hu12_baseflowindex_mean + 
                  scale(upstream_lakes_4ha_area_ha) + 
                  (1 | hu2), family = list(family = "beta", link = "logit"), 
                data = nes)

mod2 <- glmmTMB(p_percent_retention ~ scale(retention_time_yr) + 
                  iws_nlcd1992_pct_81 + 
                  iws_nlcd1992_pct_82 + 
                  iws_streamdensity_streams_density_mperha +
                  hu12_baseflowindex_mean + 
                  scale(upstream_lakes_4ha_area_ha) + 
                  (1 | edu), family = list(family = "beta", link = "logit"), 
                data = nes)

mod3 <- glmmTMB(p_percent_retention ~ scale(retention_time_yr) + 
                  iws_nlcd1992_pct_81 + 
                  iws_nlcd1992_pct_82 + 
                  iws_streamdensity_streams_density_mperha +
                  hu12_baseflowindex_mean + 
                  scale(upstream_lakes_4ha_area_ha) + 
                  (1 | hu4), family = list(family = "beta", link = "logit"), 
                data = nes)

AICtab(mod1, mod2, mod3)

# linear non-Vollenweider ####
mod1 <- lm(car::logit(p_percent_retention) ~
            scale(retention_time_yr) + 
            iws_nlcd1992_pct_81 + 
            iws_nlcd1992_pct_82 + 
            iws_streamdensity_streams_density_mperha +
            hu12_baseflowindex_mean + 
            scale(upstream_lakes_4ha_area_ha), data = nes)

plot(fitted(mod1), residuals(mod1))
plot(car::logit(nes$p_percent_retention), fitted(mod1)); abline(0, 1)

mod1_aic <- MASS::stepAIC(mod1)
plot(fitted(mod1_aic), residuals(mod1_aic))
plot(car::logit(nes$p_percent_retention), fitted(mod1_aic)); abline(0, 1)
plot(scale(nes$retention_time_yr), car::logit(nes$p_percent_retention)); abline(mod1_aic)

mod1_poly <- lm(car::logit(p_percent_retention) ~
                poly(as.numeric(scale(nes$retention_time_yr)), 2) + 
                scale(upstream_lakes_4ha_area_ha), data = nes)

summary(mod1_poly)
plot(mod1_poly)
plot(car::logit(nes$p_percent_retention), fitted(mod1_poly)); abline(0, 1)
plot(scale(nes$retention_time_yr), car::logit(nes$p_percent_retention)) 
nes_predict <- data.frame(x = scale(nes$retention_time_yr), y = predict(mod1_poly))
nes_predict <- nes_predict[order(nes_predict$x),]
lines(nes_predict)

mod2 <- lmer(p_percent_retention ~
               log(retention_time_yr) + 
               iws_nlcd1992_pct_81 + 
               iws_nlcd1992_pct_82 + 
               iws_streamdensity_streams_density_mperha + 
               hu12_baseflowindex_mean +
               scale(upstream_lakes_4ha_area_ha) + 
               (1 | edu), data = nes, REML = FALSE)

mod3 <- lmer(p_percent_retention ~
             log(retention_time_yr) + 
             iws_nlcd1992_pct_81 + 
             iws_nlcd1992_pct_82 + 
             iws_streamdensity_streams_density_mperha + 
             hu12_baseflowindex_mean +
             scale(upstream_lakes_4ha_area_ha) + 
               (1 + log(retention_time_yr)| edu), data = nes, REML = FALSE)

mod4 <- lmer(p_percent_retention ~
               log(retention_time_yr) + 
               iws_nlcd1992_pct_81 + 
               iws_nlcd1992_pct_82 + 
               iws_streamdensity_streams_density_mperha + 
               hu12_baseflowindex_mean +
               scale(upstream_lakes_4ha_area_ha) + 
               (1 | hu4_zoneid), data = nes, REML = FALSE)

mod5 <- lmer(p_percent_retention ~
               log(retention_time_yr) + 
               iws_nlcd1992_pct_81 + 
               iws_nlcd1992_pct_82 + 
               iws_streamdensity_streams_density_mperha + 
               hu12_baseflowindex_mean +
               scale(upstream_lakes_4ha_area_ha) + 
               (1 + log(retention_time_yr)| hu4_zoneid), 
             data = nes, REML = FALSE)

mod6 <- lmer(p_percent_retention ~
               log(retention_time_yr) + 
               iws_nlcd1992_pct_81 + 
               iws_nlcd1992_pct_82 + 
               iws_streamdensity_streams_density_mperha + 
               hu12_baseflowindex_mean +
               scale(upstream_lakes_4ha_area_ha) + 
               (log(retention_time_yr) || hu4_zoneid), 
             data = nes, REML = FALSE)

mod7 <- lmer(p_percent_retention ~
               log(retention_time_yr) + 
               iws_nlcd1992_pct_81 + 
               iws_nlcd1992_pct_82 + 
               iws_streamdensity_streams_density_mperha + 
               hu12_baseflowindex_mean +
               scale(upstream_lakes_4ha_area_ha) + 
               (1 | hu2), data = nes, REML = FALSE)

bbmle::AICtab(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

summary(mod7)

mod7 <- lmer(p_percent_retention ~
               log(retention_time_yr) + 
               iws_nlcd1992_pct_82 + 
               hu12_baseflowindex_mean +
               (1 | hu2), data = nes, REML = FALSE)

plot(mod6)
plot(resid(mod6) ~ nes$retention_time_yr)

plot(resid(mod6) ~ fitted(mod6))

plot(fitted(mod6) ~ nes$p_percent_retention, 
     xlim = c(0, 100), 
     ylim = c(0, 100))
abline(a = 0, b = 1)

table(nes$edu)

xyplot(log(p_percent_retention) ~ log(retention_time_yr) | edu, nes, aspect = "xy",
       layout = c(2, 2), type = c("g", "p", "r"),
       index.cond = function(x, y) coef(lm(y ~ x))[2],
       xlab = "log(Retention Time (yr))",
       ylab = "log(P Retention (%))",
       as.table = TRUE)

nes$flag <- NA

my_pal <- rep("#000000", nrow(nes))
my_pal[nes$edu %in% c("75", "73", "10", "25", "61")] <- "#ce4d1e"
ggplot(data = nes) + 
  geom_point(aes(x = log(retention_time_yr), 
                 y = log(p_percent_retention), 
                 color = edu)) + 
  scale_color_manual(values = my_pal)


fit <- lmer(log(p_percent_retention) ~ 
              log(retention_time_yr) + lake_type + 
              (log(retention_time_yr) | edu), data = nes)

library(merTools)

shinyMer(fit, simData = nes)

library(rstan)

N <- nrow(nes)
fit <- stan("01_Chapter/vollenweider.stan", 
            data = list(N = N, rp = nes$p_percent_retention, tw = nes$retention_time_yr), iter = 1000, chains = 4)

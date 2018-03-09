source("01_Chapter/scripts/99_utils.R")
source("01_Chapter/scripts/01_prepdata.R")

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
  part_pred_plot(nes_iws, readRDS("../data/md_vollenweider.rds"), 
                 9, "Max Depth", xl = FALSE), ggplot() + geom_blank(),
  part_pred_plot(nes_iws, readRDS("../data/lc_vollenweider.rds"), 
                 8, "Lake Connection", xl = FALSE), ggplot() + geom_blank(),
  part_pred_plot(nes_iws, readRDS("../data/iws/cd_vollenweider.rds"),
                 7, "Closest Lake Distance", xl = FALSE), 
  part_pred_plot(nes_nws, readRDS("../data/nws/cd_vollenweider.rds"),
                 7, "Closest Lake Distance", xl = FALSE, yl = FALSE),
  part_pred_plot(nes_iws, readRDS("../data/iws/sr_vollenweider.rds"),
                 6, "Stream Order Ratio"), 
  part_pred_plot(nes_nws, readRDS("../data/nws/sr_vollenweider.rds"),
                 6, "Stream Order Ratio", yl = FALSE),
  nrow = 4, ncol = 2)

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
  extract_coefs2(readRDS("../data/iws/ll_vollenweider.rds"),
                 c("k1_sd", "k2_sd")),
  extract_coefs2(readRDS("../data/iws/sr_vollenweider.rds"), 
                 c("k1_sr", "k2_sr")))
 
  #               extract_coefs2(fit_sd, c("k1_sd", "k2_sd")),
  #               extract_coefs2(fit_wc, c("k1_wc", "k2_wc")), 


fit_df_iws$key <- factor(fit_df_iws$key, levels = c("k", 
                                            "k1_cd", "k2_cd",
                                            "k1_sr", "k2_sr",
                                            "k1_ll", "k2_ll",
                                            "k1_bf", "k2_bf",
                                            "k1_sd", "k2_sd",
                                            "k1_lc", "k2_lc",
                                            "k1_md", "k2_md"))
#                                             "k1_wc", "k2_wc",
#                                             "k1_bf", "k2_bf",

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
  extract_coefs2(readRDS("../data/nws/cd_vollenweider.rds"), 
                 c("k1_ll", "k2_ll")),
  extract_coefs2(readRDS("../data/nws/la_vollenweider.rds"), 
                 c("k1_la", "k2_la")),
  extract_coefs2(readRDS("../data/md_vollenweider.rds"), 
                 c("k1_md", "k2_md")),
  extract_coefs2(readRDS("../data/nws/sr_vollenweider.rds"), 
                 c("k1_sr", "k2_sr")))

fit_df_nws$key <- factor(fit_df_nws$key, levels = c("k",
                                                    "k1_sr", "k2_sr",
                                                    "k1_md", "k2_md",
                                                    "k1_lc", "k2_lc",
                                                    "k1_la", "k2_la",
                                                    "k1_ll", "k2_ll",
                                                    "k1_cd", "k2_cd"))
#                                             "k1_sd", "k2_sd",

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

ggplot(nes) + geom_point(aes(x = retention_time_yr,
                             y = p_percent_retention,
                             color = maxdepth)) + 
  geom_vline(aes(xintercept = 0.482)) + 
  geom_vline(aes(xintercept = 3.6)) 

ggplot(nes) + geom_point(aes(x = retention_time_yr,
                             y = p_percent_retention,
                             color = link_length))

ggplot(nes) + geom_point(aes(x = retention_time_yr,
                             y = p_percent_retention,
                             color = iws_wl_allwetlandsdissolved_overlapping_area_pct)) + 
  theme(legend.position = "NA")

ggplot(dplyr::filter(nes, iws_lakes_overlapping_area_pct < 1.9)) + 
  geom_point(aes(x = retention_time_yr,
                 y = p_percent_retention,
                 color = log(iws_lakes_overlapping_area_pct + 
                               iws_wl_allwetlandsdissolved_overlapping_area_pct))) + 
  theme(legend.position = "NA")

ggplot(nes) + geom_point(aes(x = retention_time_yr,
                             y = p_percent_retention,
                             color = upstream_lakes_4ha_area_ha))

ggplot(nes) + geom_point(aes(x = closest_lake_distance,
                             y = p_percent_retention,
                             color = stream_order_ratio))

nes[is.na(nes$link_length),"link_length"] <- max(nes$link_length, na.rm = TRUE)
ggplot(nes) + geom_point(aes(x = retention_time_yr,
                             y = p_percent_retention,
                             color = 1/link_length))

plot(nes$retention_time_yr, nes$stream_order_ratio)
plot(nes$iws_streamdensity_streams_density_mperha, nes$link_length)

library(randomForest)
library(randomForestExplainer)

RF1 = randomForest(nes_rf$p_perce ~ ., data = nes_rf
                   , importance = TRUE
                   , na.action = na.omit 
                   , ntree = 1000)

varImpPlot(RF1, scale = FALSE)
# partialPlot(x = RF1, pred.data = test, x.var = maxdept, n.pt = 50)

plot_predict_interaction(RF1, test, "retenti", "link_le")
# explain_forest(RF1, interactions = TRUE, data = test)

library(heatR)
cormat <- cor(dplyr::select(test, -lakecon), use = 'pair')
corrheat(cormat, psychOptions = list(fm = 'ml', rotate = 'promax'))

plot(nes_rf$retenti, nes_rf$p_perce)
points(nes_rf$retenti, predict(ctree(nes_rf$p_perce ~ ., data = nes_rf, 
                                     controls = ctree_control( 
                                       minsplit = 5, 
                                       mincriterion = 0.7, 
                                       mtry = 0))), col = "red", pch = 21)

test <- data.frame(link_le = nes_rf$link_le, 
                   retention_resid = unlist(predict(
                     ctree(nes_rf$p_perce ~ retenti, 
                           data = nes_rf, 
                           controls = ctree_control( 
                             minsplit = 5, 
                             mincriterion = 0.7, 
                             mtry = 0)))) - nes_rf$p_perce)
fit <- lm(test$nes_rf.p_perce ~ test$link_le)
plot(test)
abline(fit)

# 
# + 
#   geom_point(data = nes[c(3,  11,  22,  29,  46,  62,  85, 136, 147, 160, 181, 184, 204, 207, 233),], 
#              aes(x = closest_lake_distance, 
#                  y = p_percent_retention), 
#              color = "red")
# plot(nes$closest_lake_distance, nes$p_percent_retention)
# identify(nes$closest_lake_distance, nes$p_percent_retention)
# c(3,  11,  22,  29,  46,  62,  85, 136, 147, 160, 181, 184, 204, 207, 233)
# 
# test <- nes[c(3,  11,  22,  29,  46,  62,  85, 136, 147, 160, 181, 184, 204, 207, 233),]
# library(sf)
# 
# nes_sf <- st_as_sf(nes, coords = c("lg_long", "lg_lat"), crs = 4326)
# test_sf <- st_as_sf(test, coords = c("lg_long", "lg_lat"), crs = 4326)
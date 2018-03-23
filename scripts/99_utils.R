
# ---- source_utils ----

library(LAGOSNE)
library(dplyr)
library(brms)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(data.tree)
library(partykit)
library(modelr)
library(tidyr)
library(ggridges)
library(USAboundaries)
library(sf)

prep_full_nes <- function(nes_x_lagos, conny_metrics_path, join_lagos_gis = FALSE){
  
  nes         <- read.csv(system.file("extdata/nes.csv", package = "nesRdata"),
                          stringsAsFactors = FALSE)
  nes_x_lagos <- read.csv(nes_x_lagos, stringsAsFactors = FALSE)
  
  if(!("storet_code" %in% names(nes_x_lagos))){ # if lagos-ne
    lg  <- lagosne_load("1.087.1")
    nes <- dplyr::left_join(nes_x_lagos, nes, by = c("nesid" = "storet_code"))
    nes <- dplyr::left_join(nes, select(                # lakes.geo
      lg$lakes.geo,
      lakeconnection, lagoslakeid, hu12_zoneid, upstream_lakes_4ha_area_ha))
    nes <- filter(nes, upstream_lakes_4ha_area_ha < 4000000)
    nes <- dplyr::left_join(nes, select(lg$lakes_limno, # lakes_limno
                                        lagoslakeid, meandepth, maxdepth))
    nes <- dplyr::left_join(nes, select(lg$iws.conn,    # iws.conn
                                        lagoslakeid, iws_streamdensity_streams_density_mperha, 
                                        iws_wl_allwetlandsdissolved_overlapping_area_pct,
                                        iws_lakes_overlapping_area_pct
    ))
    nes <- dplyr::left_join(nes, select(lg$hu12.chag,   # hu12.chag
                                        hu12_zoneid, hu12_baseflowindex_mean))
    nes <- dplyr::left_join(nes, select(lg$iws.lulc,    # iws.lulc
                                        lagoslakeid, iws_nlcd1992_pct_81, iws_nlcd1992_pct_82))
    # nes <- dplyr::left_join(nes, select(lg$iws,         # iws
    #                                     lagoslakeid, iws_ha))
    nes <- dplyr::left_join(nes, select(lg$locus,       # locus
                                        lagoslakeid, lake_area_ha))
    
    connectivity_metrics <- read.csv(conny_metrics_path)
    nes <- dplyr::left_join(nes, connectivity_metrics)
    
    nes <- dplyr::filter(nes, 
                         lagoslakeid != 6208) # strange outlier
    
    nes <- dplyr::filter(nes, 
                         !(nesid %in% c("1909", "1907", "1839", 
                                        "2314", "2314", "4224")))
    
  }else{ # lagos-us
    nes <- dplyr::left_join(nes, nes_x_lagos)
  }
  
  nes <- dplyr::filter(nes,
                       !is.na(retention_time),
                       !is.na(retention_time_units),
                       !is.na(p_percent_retention), 
                       p_percent_retention > 3,
                       surface_area < 1000, 
                       mean_depth < 70)
  nes$retention_time_yr <- nes$retention_time
  nes$retention_time_yr[nes$retention_time_units == "days"] <-
    nes$retention_time_yr[nes$retention_time_units == "days"] / 365
  nes <- dplyr::filter(nes, retention_time_yr < 20)
  nes <- dplyr::filter(nes, retention_time_yr > 0.02739726) # 10 days
  if("LakeConnectivity" %in% names(nes)){
    nes$lakeconnection <- nes$LakeConnectivity
  }
  nes <- dplyr::filter(nes, lakeconnection != "Isolated",
                       lakeconnection != "Headwater")
  nes$lakeconnection <- droplevels(nes$lakeconnection)
  nes$p_percent_retention <- nes$p_percent_retention / 100
  
  nes$tp_in <- nesRdata:::calculate_tp_in(nes)
  
  if(join_lagos_gis){
    nes
  }
  nes
}

theme_pred <- function(){
  theme(legend.position = "na")
}

full_model <- function(nes){
  brm(bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)),
         k + x ~ 1, nl = TRUE),
      data = nes,
      prior = c(prior(normal(1.3, 0.1),  lb = 0, nlpar = k),
                prior(normal(0.45, 0.1), lb = 0, nlpar = x)),
      iter = 4000)
}

full_pred_plot <- function(nes, fit){
  test <- tidybayes::add_fitted_samples(nes, fit)
  
  ggplot() + 
    stat_lineribbon(data = test, 
                             aes(x = retention_time_yr, 
                                 y = estimate), .prob = c(0.95)) + 
    geom_point(data = nes, aes(x = retention_time_yr, 
                               y = p_percent_retention)) + 
    scale_x_log10() + ylim(0, 1) +
    scale_fill_discrete("grey") + 
    theme_pred() + 
    ylab("P Retention (%)") + xlab("Residence Time (yr)")
}

part_model <- function(nes, part_var = NA, lower = NA, mid = NA, 
                       upper = NA, na_var = NA, fac = NA){
  
  if(!is.na(fac)){
    nes$part <- as.factor(as.numeric(nes[, fac]))
  }else{
    nes$part <- factor(as.numeric(cut(nes[, part_var], 
                                    c(lower, mid, upper), 
                                    include.highest = TRUE)))
    nes$part[is.na(nes$part)] <- na_var
  }
  
  res <- brm(
    bf(p_percent_retention ~ 1 - (1 / (1 + k * retention_time_yr ^ x)),
         x ~ 1, k ~ 0 + part, nl = TRUE),
      data = nes,
      prior = c(prior(normal(1.3, 0.1),  lb = 0, nlpar = k),
                prior(normal(0.45, 0.1), lb = 0, nlpar = x)),
      iter = 4000)
  res$part <- nes$part
  res
}

part_pred_plot <- function(nes, fit, ind, title, xl = TRUE, yl = TRUE){
  nes$part <- fit$part
  test     <- data_grid(nes, part, 
                     retention_time_yr = seq_range(retention_time_yr, n = 100))
  test     <- tidybayes::add_fitted_samples(dplyr::filter(test, !is.na(part)), 
                                            fit)
  test     <- split(test, f = test$part)
  test     <- lapply(test, function(x) group_by(x, retention_time_yr))
  
  browser()
  
  gg_format <- function(gg){
    gg <- gg + geom_point(data = nes, aes(x = retention_time_yr, 
                               y = p_percent_retention)) +
      theme(legend.position = "none") + 
      ylim(0, 1) + scale_x_log10() + ggtitle(title) +
      ylab("P Retention (%)") + xlab("Residence Time (yr)")
    
    if(!xl){
      gg <- gg + theme(axis.title.x = element_blank())
    }
    
    if(!yl){
      gg <- gg + theme(axis.title.y = element_blank())
    }
  }
  
  # https://stackoverflow.com/a/37623958/3362993
  gg_1 <- ggplot() + 
    stat_summary(data = test[1][[1]], geom = "ribbon", 
                 aes(x = retention_time_yr, y = estimate, 
                     fill = forcats::fct_rev(ordered(...prob..))), 
                 fun.data = median_qi, fun.args = list(.prob = 0.95), alpha = 0.3) +
    stat_summary(data = test[1][[1]], aes(x = retention_time_yr, y = estimate), 
                 fun.y = median, geom = "line", color = viridis::viridis(1, begin = 1), 
                 size = 0.5) +
    scale_fill_viridis_d(begin = 1)
  
  gg_2 <- ggplot() +
    stat_summary(data = test[2][[1]], geom = "ribbon", 
                 aes(x = retention_time_yr, y = estimate, 
                     fill = forcats::fct_rev(ordered(...prob..))), 
                 fun.data = median_qi, fun.args = list(.prob = c(0.95)), alpha = 0.3) + 
    stat_summary(data = test[2][[1]], aes(x = retention_time_yr, y = estimate), 
                 fun.y = median, geom = "line", color = viridis::viridis(1, begin = 0.5), 
                 size = 0.5) +
    scale_fill_viridis_d(begin = 0.5) + 
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(
            color = "transparent", 
            fill = "transparent"),
          panel.border = element_blank())
    
  gg_1 <- gg_format(gg_1)
  gg_2 <- gg_format(gg_2)
  g1 <- ggplotGrob(gg_1)
  g2 <- ggplotGrob(gg_2)
  g2 <- gtable::gtable_filter(g2, "panel")
  
  gg <- gtable::gtable_add_grob(g1, g2, t = 6, l = 4)
  
  grid::grid.draw(gg)
}

pprint_equalities <- function(node, digits = 1){
  node_raw <- node$parent$splitlevels[node$position]
  node_raw <- strsplit(node_raw,  " ")[[1]]
  node_raw[[2]] <- as.character(round(as.numeric(node_raw[[2]]), digits = digits))
  node$splitLevel <- paste(node_raw, collapse = " ")
}

plot_tree <- function(tree, title){
  tree2 <- as.Node(tree)
  SetNodeStyle(tree2,
               label = function(node) paste0(node$name, ": ", node$splitname),
               tooltip = function(node) paste0(nrow(node$data), " observations"),
               fontname = "helvetica")
  
  tree2$Do(function(node) pprint_equalities(node), filterFun = isNotRoot)
  
  SetEdgeStyle(tree2,
               arrowhead = "none",
               label = function(node) node$splitLevel,
               fontname = "helvetica",
               penwidth = function(node){
                 12 * nrow(node$data)/nrow(node$root$data)
               })
  
  graph <- ToDiagrammeRGraph(tree2, "climb", NULL)
  DiagrammeR::render_graph(graph, output = "graph", title = title)
}

make_rf_dt <- function(nes){
  nes_rf <- dplyr::select(nes, 
                          lakeconnection:retention_time_yr, 
                          -meandepth, -hu12_zoneid, p_percent_retention)
  # nes_rf$retention_time_yr <- log(nes_rf$retention_time_yr)
  nes_rf          <- nes_rf[complete.cases(nes_rf),]
  nes_rf_name_key <- cbind(names(nes_rf), 
                           make.names(substring(names(nes_rf), 0, 7), 
                                      unique = TRUE))
  names(nes_rf)   <- make.names(substring(names(nes_rf), 0, 7), 
                                unique = TRUE)
  
  list(nes_rf = nes_rf, 
       nes_rf_name_key = nes_rf_name_key)
}

inv_inv_closest <- function(x, nes_rf) {
  ((( 1 / x) * sd(nes_rf$closest, na.rm = TRUE)) +
               mean(nes_rf$closest, na.rm = TRUE))
}

get_tree_split <- function(tree){
  tree2 <- as.Node(tree)
  SetNodeStyle(tree2,
               label = function(node) paste0(node$name, ": ", node$splitname),
               tooltip = function(node) paste0(nrow(node$data), " observations"),
               fontname = "helvetica")
  res <- capture.output(tree2$Get(function(node) print(node$parent$splitlevels), 
                   filterFun = isLeaf, simplify = "regular"))
  as.numeric(gsub('"', '', 
                  gsub("\\\\", "", strsplit(res[1], " ")[[1]][3], 
                                fixed = TRUE)))
}

extract_coefs <- function(fit = NA, col_names){
  if(any(!is.na(fit))){
    fit_part_df <- as.data.frame(fit)
    fit_part_df <- dplyr::select(fit_part_df, `b_k_Intercept`)
    names(fit_part_df) <- col_names
    tidyr::gather(fit_part_df, value = "coef")
  }else{
    data.frame(key = col_names, coef = rep(NA, length(col_names)))
  }
}

extract_coefs2 <- function(fit = NA, col_names){
  if(any(!is.na(fit))){
    fit_part_df <- as.data.frame(fit)
    fit_part_df <- dplyr::select(fit_part_df, b_k_part1, b_k_part2)
    names(fit_part_df) <- col_names
    tidyr::gather(fit_part_df, value = "coef")
  }else{ 
    data.frame(key = col_names, coef = rep(NA, length(col_names)))
  }
}

delta_k <- function(rds_path){
  dt <- readRDS(rds_path)
  dt <- extract_coefs2(dt, c("k1", "k2"))
  data.frame(pnames2 = substring(basename(rds_path), 0, 2), d_k = 
  abs(median(dplyr::filter(dt, key == "k2")$coef) - 
    median(dplyr::filter(dt, key == "k1")$coef)))
}

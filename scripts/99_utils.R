
# ---- source_utils ----

library(nlaR)
library(pheatmap)
library(superheat)
library(corrr)
library(LAGOSNE)
library(nesRdata)
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
library(magrittr)
library(purrr)

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
  nes$lakeconnection      <- droplevels(nes$lakeconnection)
  nes$p_percent_retention <- nes$p_percent_retention / 100
  
  nes$tp_in <- nesRdata:::calculate_tp_in(nes)
  
  if(join_lagos_gis){
    nes
  }
  nes
}

theme_pred <- function(){
  theme(legend.position = "na", 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 10, face = "bold", color = "black", hjust = 0))
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
                                 y = estimate), .prob = c(0.95), 
                    fill = "grey", color = "black", alpha = 0.5) + 
    geom_point(data = nes, aes(x = retention_time_yr, 
                               y = p_percent_retention), 
               size = 0.6, color = "black") + 
    scale_x_log10() + ylim(0, 1) +
    scale_fill_discrete("grey") + 
    theme_pred() + 
    ylab("P retention") + xlab("Residence time (yr)") + ggtitle("A.")
}

# na_var: value assigned to rows with NAs
# fac: group assignment for non-numeric variables like lakeconnection?
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

part_pred_plot <- function(nes, fit, ind, title, xl = TRUE, yl = TRUE, add_legend = TRUE, 
                           rev_legend = FALSE, legend_title = "connectivity"){
  nes$part <- fit$part
  test     <- data_grid(nes, part, 
                     retention_time_yr = seq_range(retention_time_yr, n = 100))
  test     <- tidybayes::add_fitted_samples(dplyr::filter(test, !is.na(part)), 
                                            fit)
  test     <- split(test, f = test$part)
  test     <- lapply(test, function(x) group_by(x, retention_time_yr))
  
  gg_format <- function(gg){
    gg <- gg + 
      theme(
        legend.position = "none", 
        plot.title = element_text(size = 10, face = "bold", color = "black", hjust = 0), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) + 
      ylim(0, 1) + scale_x_log10() + ggtitle(title) +
      ylab("P retention") + xlab("Residence time (yr)")
    
    if(!xl){
      gg <- gg + theme(axis.title.x = element_blank(), 
                       plot.margin = margin(t = 1.2, r = 0.2, l = 0.2, unit = "cm"))
    }
    
    if(!yl){
      gg <- gg + theme(axis.title.y = element_blank())
    }
    gg
  }
  # https://stackoverflow.com/a/37623958/3362993
  gg_1 <- ggplot() + 
    geom_point(data = nes, aes(x = retention_time_yr, 
                               y = p_percent_retention, color = part, shape = part), 
               size = 0.8) +
    scale_color_manual(values = c(viridis::viridis(1, begin = 0.5), 
                                  viridis::viridis(1, begin = 0))) +
    stat_summary(data = test[1][[1]], geom = "ribbon", 
                 aes(x = retention_time_yr, y = estimate, 
                     fill = forcats::fct_rev(ordered(...prob..))), 
                 fun.data = median_qi, fun.args = list(.prob = 0.95), 
                 alpha = 0.25) +
    stat_summary(data = test[1][[1]], aes(x = retention_time_yr, y = estimate),
                 fun.y = median, geom = "line", 
                 color = viridis::viridis(1, begin = 0.5),
                 size = 0.5) +
    scale_fill_viridis_d(begin = 0.5)
  
  gg_2 <- ggplot() +
    stat_summary(data = test[2][[1]], geom = "ribbon", 
                 aes(x = retention_time_yr, y = estimate, 
                     fill = forcats::fct_rev(ordered(...prob..))), 
                 fun.data = median_qi, fun.args = list(.prob = c(0.95)), alpha = 0.25) + 
    stat_summary(data = test[2][[1]], aes(x = retention_time_yr, y = estimate),
                 fun.y = median, geom = "line", color = viridis::viridis(1, begin = 0),
                 size = 0.5, fun.args = list(alpha = 0.15)) +
    scale_fill_viridis_d(begin = 0) + 
    theme(panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(
            color = "transparent", 
            fill = "transparent"),
          panel.border = element_blank())
    
  gg_1 <- gg_format(gg_1)
  gg_2 <- gg_format(gg_2)
  
  # create dummy plot to pull legend
  if(add_legend){
    dummy_df <- data.frame(connectivity = c("low", "high"), 
                           foobar = c(1,2), stringsAsFactors = FALSE)
    if(rev_legend){
      dummy_df <- data.frame(connectivity = c("high", "low"), 
                             foobar = c(1,2), stringsAsFactors = FALSE)
      dummy_df$connectivity <- factor(dummy_df$connectivity, levels = c("low", "high"))
    }

    dummy_legend <- ggplot(dummy_df) + 
      geom_line(aes(x = connectivity, y = foobar, color = connectivity)) + 
      scale_color_manual(values = c(viridis::viridis(1, begin = 0), 
                                    viridis::viridis(1, begin = 0.5))) + 
      theme(legend.position = c(0.7, 0.95), 
            legend.title = element_text(size = 7, face = "bold"), 
            legend.text = element_text(size = 7), 
            legend.background = element_rect(fill = "white", 
                                                 color = "white", 
                                                 size = 0), 
            legend.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"),
            legend.box.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"), 
            legend.spacing = unit(0.01, "cm")) + 
      labs(color = legend_title)
    
    dl <- ggplotGrob(dummy_legend)
    legend <- gtable::gtable_filter(dl, "guide-box")
    gg_1 <- gg_1 + annotation_custom(grob = legend, 
                                     ymax = 0.15, 
                                     xmax = 0.89)
  }
  
  g1 <- ggplotGrob(gg_1)
  g2 <- ggplotGrob(gg_2)
  # gtable::gtable_show_layout(g1)
  # print(g2)
  g2 <- gtable::gtable_filter(g2, "panel")
  
  gg <- gtable::gtable_add_grob(g1, g2, t = 7, l = 5)
  
  # grid::grid.newpage()
  # grid::grid.draw(gg)
  
  gg
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

break_word <- function(x, max.len){
  
  x_split <- stringr::str_split(x, " ")[[1]]
  n_chars <- sapply(x_split, nchar)
  
  first_line  <- cumsum(n_chars + 1) < max.len
  second_line <- !first_line
  
  res <- paste(x_split[first_line], collapse = " ")
  paste(res, "\n", paste(x_split[second_line], collapse = " "))
}

get_sub <- function(dt, col_name, split_value, higher_connectivity){
  if(higher_connectivity == "higher"){
    coordinatize(
      dplyr::filter(dt, UQ(rlang::sym(as.character(col_name))) <= split_value), 
      "lat", "long")
  }else{
    coordinatize(
      dplyr::filter(dt, UQ(rlang::sym(as.character(col_name))) > split_value), 
      "lat", "long")
  }
}

get_sub2 <- function(dt, col_name, split_value){
  if(is.na(split_value)){
    table(dt[,col_name]) 
  }else{
    dplyr::filter(dt, UQ(rlang::sym(as.character(col_name))) <= split_value)
  }
}

signif_star <- function(x){
  if(!is.na(x)){
    if(x){
      "*"
    }else{
      ""
    }
  }else{
    ""
  }
}

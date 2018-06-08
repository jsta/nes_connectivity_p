cmdargs <- commandArgs(trailingOnly = TRUE)

# ---- load_packages ----

# cran
library(concaveman)
library(sf)
library(dplyr)
library(mapview)
library(sp)
library(ggsn)

# github
library(vapour)
library(nhdR)
library(LAGOSNE)
library(streamnet)

# ---- source_functions ----

pull_iws <- function(lagoslakeid, maxsteps = 15){
  gdb_path <- path.expand("~/.local/share/LAGOS-GIS/lagos-ne_gis.gpkg")
  crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  
  layer_name <- "LAGOS_NE_All_Lakes_4ha"
  query      <- paste0("SELECT * FROM ", layer_name, 
                       " WHERE lagoslakeid = '", 
                       lagoslakeid, "'")
  
  lakepoly_geom <- st_as_sfc(vapour_read_geometry(gdb_path, sql = query))
  lakepoly      <- st_sf(vapour_read_attributes(gdb_path, sql = query), 
                         geometry = lakepoly_geom, crs = crs)
  lakepoly <- st_zm(lakepoly)
  
  layer_name <- "IWS"
  query      <- paste0("SELECT * FROM ", layer_name, 
                       " WHERE lagoslakeid = '", 
                       lagoslakeid, "'")
  
  poly_geom <- st_as_sfc(vapour_read_geometry(gdb_path, sql = query))
  poly      <- st_sf(vapour_read_attributes(gdb_path, sql = query), 
                     geometry = poly_geom, crs = crs)
  poly      <- st_zm(poly)
  poly      <- st_union(poly, lakepoly)
  
  lines <- nhd_plus_query(poly = poly, dsn = "NHDFlowLine", 
                            approve_all_dl = TRUE)$sp$NHDFlowLine
  
  lg        <- lagosne_load("1.087.1")
  wb_coords <- as.numeric(
    lake_info(lg, lagoslakeid = lagoslakeid)[,c("nhd_long", "nhd_lat")])

  names(lines) <- tolower(names(lines))
  extract_network(wb_coords[1], wb_coords[2], lines = lines, maxsteps = maxsteps)
}

pull_nws <- function(lagoslakeid, map = FALSE, maxsteps = 15, buffer_dist = 0.01){
    print(lagoslakeid)
    lg        <- lagosne_load("1.087.1")
    wb_coords <- as.numeric(
      LAGOSNE::lake_info(lg, lagoslakeid = lagoslakeid)[,c("nhd_long", "nhd_lat")])
    
    extract_network(wb_coords[1], wb_coords[2], maxsteps = maxsteps, 
                    buffer_dist = buffer_dist, approve_all_dl = TRUE)
}

pull_lakes <- function(lines){
  poly  <- concaveman::concaveman(st_cast(lines, "POINT"))$polygons
  nhd_plus_query(poly = poly, 
                 dsn = "NHDWaterbody", approve_all_dl = TRUE)$sp$NHDWaterbody
}
  
run_lake <- function(lagoslakeid, lines, lakes){
    print(lagoslakeid)
    calc_metrics(lines, lakes)
}

pull_nws_metrics <- function(lakes){
  
  if(nrow(lakes) == 0){
    list(baseflow = NA, 
         stream_density = NA)
  }else{
    # pull llids
    lg <- lagosne_load("1.087.1")
    poly                <- concaveman::concaveman(
                                st_cast(lakes, "POINT"))$polygons
    lg_sf               <- st_transform(coordinatize(lg$locus), st_crs(poly))
    poly_lake_intersect <- st_intersects(lg_sf, poly)
    res                 <- lg_sf[which(unlist(lapply(poly_lake_intersect, 
                          function(x) length(x) > 0))),]
    res <- data.frame(lagoslakeid = res$lagoslakeid, stringsAsFactors = FALSE)
    
    # join metrics
    
    # baseflowindex - area-weighted average
    hu12s <- dplyr::filter(lg$locus, lagoslakeid %in% res$lagoslakeid)
    
    hu12_ids <- unique(hu12s$hu12_zoneid)
    hu12_areas <- dplyr::filter(lg$hu12, hu12_zoneid %in% hu12_ids)$hu12_ha
    hu12_baseflow <- dplyr::filter(
      dplyr::select(lg$hu12.chag, hu12_zoneid, hu12_baseflowindex_mean), 
      hu12_zoneid %in% hu12_ids)$hu12_baseflowindex_mean 
      
    baseflow <- sum(hu12_areas * hu12_baseflow) / sum(hu12_areas)
    
    # stream density - area-weighted average
    res <- dplyr::left_join(res, dplyr::select(lg$iws, lagoslakeid, iws_ha))
    res <- dplyr::left_join(res, dplyr::select(lg$iws.conn, lagoslakeid, 
                                       iws_streamdensity_streams_density_mperha))
    res <- res[!is.na(res$iws_ha),]
    
    nws_ha <- sum(res$iws_ha)
    
    stream_density <- 
      sum(res$iws_ha * res$iws_streamdensity_streams_density_mperha) / nws_ha
    
    list(baseflow = baseflow, 
         stream_density = stream_density, 
         nws_ha = nws_ha)
    }
  }

# ---- execute commands ----
flist <-c("data/nes_x_lagos-ne.csv", 
          "../data/nes_x_lagos-ne.csv", 
          "nes_x_lagos-ne.csv")
nes_x_lagos <- read.csv(flist[which.max(sapply(flist, file.exists))], 
                        stringsAsFactors = FALSE)

# all(nes_x_lagos$lagoslakeid %in% lg$locus$lagoslakeid)

# lagoslakeid <- cmdargs <- 5514                         # small iws
# lagoslakeid <- nes_x_lagos$lagoslakeid[30]  # very large iws
# lagoslakeid <- 4351                         # na conny value in prev analysis
# lagoslakeid <- 4566                         # all na in prev analysis
# lagoslakeid <- 7195                         # straddles utm zone boudary
# lagoslakeid <- 4686
# lagoslakeid <- 4314
# lagoslakeid <- 4819
# lagoslakeid <- 6753
# lagoslakeid <- 2639
# lagoslakeid <- 2227
# lagoslakeid <- 2197
# lagoslakeid <- 3523
# lagoslakeid <- 6302
# lagoslakeid <- cmdargs <- 1906
# lagoslakeid <- cmdargs <- 8016
# lagoslakeid <- cmdargs <- x <- 6076
# iws <- pull_iws(lagoslakeid, maxsteps = Inf)
# nws <- pull_nws(lagoslakeid, maxsteps = Inf)
# iws_lakes <- pull_lakes(iws)
# nws_lakes <- pull_lakes(nws)
# nws_avgs  <- pull_nws_metrics(nws_lakes)
#  
# mapview(nws) + mapview(nws_lakes) + mapview(iws, color = "green")
#  
# iws_metrics <- run_lake(lagoslakeid, lines = iws, lakes = iws_lakes)
# nws_metrics <- run_lake(lagoslakeid, lines = nws, lakes = nws_lakes)
# tryCatch({
#   nws_metrics <- run_lake(lagoslakeid, lines = nws, lakes = nws_lakes)
# }, error = function(e) {
#   list(avg_link_length = NA, 
#        stream_order_ration = NA,
#        closest_lake_distance = NA,
#        num_up_lakes = NA)
# })
# cbind(nws_metrics, iws_metrics)

if(length(cmdargs) > 0){
  nes_x_lagos <- dplyr::filter(nes_x_lagos, lagoslakeid == cmdargs)
}else{
  # don't run some strange lakes
  nes_x_lagos <- dplyr::filter(nes_x_lagos, 
                               lagoslakeid != 2639,
                               lagoslakeid != 2227,
                               lagoslakeid != 2197,
			                         lagoslakeid != 4930)
}
  
res   <- lapply(nes_x_lagos$lagoslakeid, function(x) {
  # add a test to exclude coastal lakes?
  tryCatch({
    lines <- pull_nws(x, maxsteps = Inf)
    }, error = function(e){
    lines <- NA
  })
      
  if(any(!is.na(lines)) & any(class(lines) != "function")){
    print(1)
    lakes    <- pull_lakes(lines)
    nws_avgs <- pull_nws_metrics(lakes)
    
    tryCatch({
      res <- run_lake(x, lines = lines, lakes = lakes)
      print(2)
      return(c(list(lagoslakeid = x), res, nws_avgs))
    }, error = function(e) {
      print(3)
      c(list(
           lagoslakeid = x,
           avg_link_length = NA, 
           stream_order_ratio = NA,
           closest_lake_distance = NA,
           num_up_lakes = NA,
           lake_area    = NA), 
        nws_avgs)
    })
  }else{
    print(4)
    list(lagoslakeid = x,
         avg_link_length = NA, 
         stream_order_ratio = NA,
         closest_lake_distance = NA,
         num_up_lakes = NA,
         lake_area    = NA,
         baseflow     = NA,
         stream_density = NA, 
         nws_ha = NA)
  }

})

if(length(cmdargs) == 0){
  res <- do.call("rbind", res)
  write.table(res, file = "connectivity_metrics_nws.csv", append = TRUE, 
              sep = ",", row.names = FALSE, col.names = TRUE)
}
  print(res)

# state_code <- "ME.rds"
# state_code <- "NH.rds"
# state_code <- "VT.rds"
# state_code <- "MA.rds"
# state_code <- "CT.rds"
# state_code <- "RI.rds"
# state_code <- "NY.rds"
state_code <- commandArgs(trailingOnly = TRUE)

library(nesRdata)
library(LAGOSNE)
library(sf)
library(dplyr)
library(mapview)
library(nhdR)

lg_raw     <- lagos_load("1.087.1")
states_raw <- sf::st_as_sf(map("state", plot = FALSE, fill = TRUE))
# nes_get("1")
nes_raw <- nes_load(1)$nes_data

state_codes  <- tolower(state.name[state.abb == gsub(".rds", "", state_code)])
states       <- states_raw[states_raw$ID %in% state_codes,]

nes_raw$nes_long <- nes_raw$long
nes_raw$nes_lat <- nes_raw$lat
nes <- st_as_sf(nes_raw, coords = c("long", "lat"), crs = 4326)
nes <- filter(nes, state == toupper(state_codes))

lg_epi <- lg_raw$epi_nutr[
            !duplicated(lg_raw$epi_nutr$lagoslakeid),]
lg     <- left_join(lg_epi, lg_raw$locus)
lg$lg_long <- lg$nhd_long
lg$lg_lat  <- lg$nhd_lat
lg     <- st_as_sf(lg, coords = c("nhd_long", "nhd_lat"), crs = 4326)

lg <- select(lg, lagoslakeid, nhdid, gnis_name, lg_long, lg_lat)
nes <- select(nes, storet_code, nes_long, nes_lat)

lg     <- st_transform(lg, 5070)
nes    <- st_transform(nes, 5070)
states <- st_transform(states, 5070)

# select nhd water bodies that intersect nes points ####

# pull waterbodies for given state
wb_raw <- nhdR::nhd_load(state.abb[which(tolower(state.name) %in% states$ID)], 
                     dsn = "NHDWaterbody", approve_all_dl = TRUE)
wb_raw <- st_transform(wb_raw, 5070)
wb <- st_buffer(wb_raw, 200)
wb <- select(wb, Permanent_Identifier, GNIS_Name)

contains_nes          <- sapply(st_intersects(wb, nes), 
                                 function(x){length(x) > 0})
intersects_nes        <- st_intersects(wb, nes)[contains_nes]
names(intersects_nes) <- which(contains_nes)

key        <- reshape2::melt(intersects_nes)
names(key) <- c("b", "a")

wb_sub <- wb[key$a,]
key$nhdid <- wb_sub$Permanent_Identifier
key$nesid <- nes$storet_code[key$b]

# in cases with multiple overlap because of buffer select largest
if(any(duplicated(key$b))){
  key$area <- st_area(wb_sub)
}
key <- group_by(key, b)
key <- data.frame(top_n(key, n = 1))
wb_sub <- wb[key$a,]

wb_sub <- as.data.frame(wb_sub)
wb_sub <- select(wb_sub, -Shape)

#####
res <- inner_join(key, data.frame(lg))
res <- select(res, -geometry, -a, -b)
res <- left_join(res, data.frame(nes), by = c("nesid" = "storet_code"))
if(any(names(res) == "area")){
  res <- select(res, -area)
}
res <- select(res, -geometry)

# test <- st_as_sf(res, coords = c("lg_long", "lg_lat"), crs = 4326)
# test2 <- st_as_sf(res, coords = c("nes_long", "nes_lat"), crs = 4326)

# mapview(test) + mapview(test2, color = "red")

# mapview(wb_sub) + mapview(nes) + mapview(lg_sub, color = "red") + mapview(lg, color = "green")

# filter(lg_raw$lakes.geo, lakes_nhdid %in% lg_sub$nhdid)$lakeconnection

saveRDS(res, state_code)

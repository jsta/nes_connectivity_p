# setwd("scripts")
source("99_utils.R")

# ---- prep_nes_lagos ----

# pre_path <- "/home/jose/Documents/Science/Dissertation/Analysis/nes/"
pre_path <- "../data"

nes_iws <- prep_full_nes(file.path(pre_path, "nes_x_lagos-ne.csv"), 
                         file.path(pre_path, "connectivity_metrics_iws.csv"))

pre_path <- "../data"
nes_nws <- prep_full_nes(file.path(pre_path, "nes_x_lagos-ne.csv"), 
                         file.path(pre_path, "connectivity_metrics_nws.csv"))

nes_rf_iws <- make_rf_dt(nes_iws)$nes_rf
# make_rf_dt(nes_iws)$nes_rf_name_key

nes_rf_nws <- make_rf_dt(nes_nws)$nes_rf
# make_rf_dt(nes_nws)$nes_rf_name_key
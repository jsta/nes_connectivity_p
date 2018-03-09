source("scripts/99_utils.R")

pre_path <- "/home/jose/Documents/Science/Dissertation/Analysis/nes/"

# ---- prep_nes_lagos ----

nes_iws <- prep_full_nes(file.path(pre_path, "nes_x_lagos-ne.csv"), 
                         file.path(pre_path, "connectivity_metrics.csv"))

nes_nws <- prep_full_nes(file.path(pre_path, "nes_x_lagos-ne.csv"), 
                         file.path(pre_path, "connectivity_metrics2.csv"))

nes_rf_iws <- make_rf_dt(nes_iws)$nes_rf

nes_rf_nws <- make_rf_dt(nes_nws)$nes_rf
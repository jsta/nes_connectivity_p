# setwd("scripts")
if(file.exists("99_utils.R")){
  source("99_utils.R")
}else{
  if(file.exists("../scripts/99_utils.R")){
    source("../scripts/99_utils.R")  
  }else{
    source("scripts/99_utils.R")  
  }
}

# ---- prep_nes_lagos ----

# Base data
if(file.exists("../data")){
  pre_path <- "../data"
}else{
  pre_path <- "data"
}

# if(!all(sapply(c("nes_iws", "nes_nws", "nes_rf_iws", "nes_rf_nws"), exists))){
  
  nes_iws <- prep_full_nes(file.path(pre_path, "nes_x_lagos-ne.csv"), 
                         file.path(pre_path, "connectivity_metrics_iws.csv"))
  
  

  nes_nws <- prep_full_nes(file.path(pre_path, "nes_x_lagos-ne.csv"), 
                         file.path(pre_path, "connectivity_metrics_nws.csv"))
  
  saveRDS(nes_iws, "nes_iws.rds")
  saveRDS(nes_nws, "nes_nws.rds")
  
  # RandomForest data
  nes_rf_iws <- make_rf_dt(nes_iws)$nes_rf
  # make_rf_dt(nes_iws)$nes_rf_name_key

  nes_rf_nws <- make_rf_dt(nes_nws)$nes_rf
  # make_rf_dt(nes_nws)$nes_rf_name_key
# }

  saveRDS(nes_rf_iws, "nes_rf_iws.rds")
  saveRDS(nes_rf_nws, "nes_rf_nws.rds")
  
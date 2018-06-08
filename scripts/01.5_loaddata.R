
if(file.exists("../data/nes_nws.rds")){
  pre_path <- "../"
}else{
  pre_path <- ""
}

nes_nws <- readRDS(paste0(pre_path, "data/nes_nws.rds"))
nes_iws <- readRDS(paste0(pre_path, "data/nes_iws.rds"))
# rds <- dir("Analysis/nes/", pattern = "rds", full.names = TRUE, include.dirs = TRUE)
args <- commandArgs(trailingOnly = TRUE)
rds <- dir(".", pattern = "rds")

dt <- lapply(rds, readRDS)
dt <- do.call("rbind", dt)

write.csv(dt, args[1], row.names = FALSE)

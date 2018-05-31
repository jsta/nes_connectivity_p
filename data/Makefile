all: connectivity_metrics_nws.csv nes_iws.rds nes_nws.rds nes_rf_iws.rds nes_rf_nws.rds
	
connectivity_metrics_nws.csv: connectivity_metrics.R
	Rscript connectivity_metrics.R

lagos:
	Rscript -e "LAGOSNE::lagosne_get('1.087.1')"

nes_iws.rds: ../scripts/01_prepdata.R ../scripts/99_utils.R
	Rscript $<

nes_nws.rds: ../scripts/01_prepdata.R ../scripts/99_utils.R
	Rscript $<

nes_rf_iws.rds: ../scripts/01_prepdata.R ../scripts/99_utils.R
	Rscript $<

nes_rf_nws.rds: ../scripts/01_prepdata.R ../scripts/99_utils.R
	Rscript $<
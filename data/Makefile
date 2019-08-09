STATES=ME CT IL IN IA MN MO NH VT MA MI NJ NY OH PA RI WI
rds_FILES=$(foreach state, $(STATES), $(state).rds)
# STATES_US=AL AZ AR CA CO CT DE FL GA ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY
# rds_FILES_US=$(foreach state, $(STATES_US), $(state).rds)

rdss: $(rds_FILES)
# rdsus: $(rds_FILES_US)

%.rds: $(patsubst %, %.rds, $(STATES))
	Rscript nes_x_lagos.R $@

# %.rds: $(patsubst %, %.rds, $(STATES_US))
#	Rscript nes_x_lagos-us.R $@

# ne and us are combined in a single target
# need to recompute us rdss manually?
nes_x_lagos-ne.csv: rdss
	Rscript nes_x_lagos_merge.R nes_x_lagos.csv
	-rm $(rds_FILES)

all: connectivity_metrics_nws.csv nes_iws.rds nes_nws.rds nes_rf_iws.rds nes_rf_nws.rds nes_x_lagos-ne.csv
	
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

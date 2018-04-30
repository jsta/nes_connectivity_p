all: connectivity_metrics_nws.csv
	
connectivity_metrics_nws.csv: connectivity_metrics.R
	Rscript connectivity_metrics.R

lagos:
	Rscript -e "LAGOSNE::lagosne_get('1.087.1')"
	
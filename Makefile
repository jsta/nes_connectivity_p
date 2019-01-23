.PHONY: figures data tables all

all: manuscript/manuscript.pdf figures tables scripts/table_1.csv

data: 
	cd data && make all

Dockerfile: manuscript/manuscript.Rmd
	Rscript -e "liftr::lift('$<', output_dir = '.')"

manuscript/manuscript.pdf: manuscript/manuscript.Rmd
	Rscript -e "rmarkdown::render('$<')"

figures: manuscript/figures.pdf manuscript/appendix.pdf manuscript/grayscale_figures.pdf

manuscript/appendix.pdf: manuscript/appendix.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	# -pdftk manuscript/appendix.pdf cat 2-end output manuscript/appendix2.pdf
	# -mv manuscript/appendix2.pdf manuscript/appendix.pdf

manuscript/figures.pdf: manuscript/figures.Rmd figures/01_conceptual_p-cycle_pt2.pdf figures/02_conny_metric_key_crop_small.pdf figures/03_iws_nws.pdf figures/04_global_vollenweider-1.pdf figures/05_partition_vollenweider-1.pdf figures/06_k-1.pdf figures/07_cor_mat_hmap-1.pdf figures/08_maps-1.pdf
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	-pdftk manuscript/figures.pdf cat 2-end output manuscript/figures2.pdf
	-mv manuscript/figures2.pdf manuscript/figures.pdf
	
manuscript/grayscale_figures.pdf: manuscript/grayscale_figures.Rmd manuscript/figures.pdf
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	-pdftk manuscript/grayscale_figures.pdf cat 2-end output manuscript/figures2.pdf
	-mv manuscript/figures2.pdf manuscript/grayscale_figures.pdf

figures/03_iws_nws.pdf: figures/03_iws_nws.tex figures/beamer-lake-fig/beamer-lake-fig.tex
	cd figures && pdflatex 03_iws_nws.tex
	pdfcrop figures/03_iws_nws.pdf figures/03_iws_nws.pdf
	
figures/02_conny_metric_key.pdf: figures/02_conny_metric_key.tex figures/beamer-lake-fig/beamer-lake-fig.tex
	cd figures && pdflatex 02_conny_metric_key.tex
	
figures/02_conny_metric_key_crop_small.pdf: figures/02_conny_metric_key.pdf
	pdfcrop $< figures/02_conny_metric_key_crop.pdf
	convert figures/02_conny_metric_key_crop.pdf -units pixelsperinch -density 300 -page ArchC figures/02_conny_metric_key_crop_small.pdf
	rm figures/02_conny_metric_key_crop.pdf
	
figures/04_global_vollenweider-1.pdf: figures/04_global_vollenweider.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/04_global_vollenweider.pdf

figures/04_global_vollenweider.Rmd: scripts/03_viz.R

figures/05_partition_vollenweider-1.pdf: figures/05_partition_vollenweider.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/05_partition_vollenweider.pdf

figures/05_partition_vollenweider.Rmd: scripts/03_viz.R 

figures/06_k-1.pdf: figures/06_k.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/06_k.pdf
	pdfcrop figures/06_k-1.pdf figures/06_k-1.pdf

figures/07_cor_mat_hmap-1.pdf: figures/07_cor_mat_hmap.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/07_cor_mat_hmap.pdf
	
figures/07_cor_mat_hmap.Rmd: scripts/03_viz.R 

figures/08_maps-1.pdf: figures/08_maps.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	pdfcrop $@ $@
	rm figures/08_maps.pdf

tables: manuscript/tables.pdf

manuscript/tables.pdf: tables/01_lake_characteristics_table.pdf tables/02_model_results_table.pdf
	pdftk $^ cat output manuscript/tables.pdf

tables/01_lake_characteristics_table.pdf: tables/01_lake_characteristics_table.Rmd scripts/04_tables.R
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	
tables/02_model_results_table.pdf: tables/02_model_results_table.Rmd scripts/table_1.csv
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	
scripts/table_1.csv: scripts/table_1.R
	Rscript scripts/table_1.R

test: 
	Rscript -e 'dir.create(rappdirs::user_data_dir(), recursive = TRUE); example("calc_metrics", package = "streamnet", run.dontrun = TRUE)'
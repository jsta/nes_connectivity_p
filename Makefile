.PHONY: figures data tables all

data: 
	cd data && make all

all: manuscript/draft.pdf figures tables scripts/table_1.csv

manuscript/draft.pdf: manuscript/draft.md
	pandoc manuscript/draft.md -o manuscript/draft.pdf

figures: manuscript/figures.pdf manuscript/appendix.pdf

manuscript/appendix.pdf: manuscript/appendix.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	-pdftk manuscript/appendix.pdf cat 2-end output manuscript/appendix2.pdf
	-mv manuscript/appendix2.pdf manuscript/appendix.pdf

manuscript/figures.pdf: manuscript/figures.Rmd figures/iws_nws.pdf figures/conny_metric_key.pdf figures/04_global_vollenweider-1.pdf figures/05_partition_vollenweider-1.pdf figures/06_k-1.pdf figures/07_cor_mat_hmap-1.pdf figures/08_maps-1.pdf
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	-pdftk manuscript/figures.pdf cat 2-end output manuscript/figures2.pdf
	-mv manuscript/figures2.pdf manuscript/figures.pdf

figures/iws_nws.pdf: figures/iws_nws.tex figures/beamer-lake-fig/beamer-lake-fig.tex
	pdflatex figures/iws_nws.tex
	pdfcrop figures/iws_nws.pdf figures/iws_nws.pdf
	
figures/conny_metric_key.pdf: figures/conny_metric_key.tex
	pdflatex figures/conny_metric_key.tex
	pdfcrop figures/conny_metric_key.pdf figures/conny_metric_key_crop.pdf
	
figures/04_global_vollenweider-1.pdf: figures/04_global_vollenweider.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/04_global_vollenweider.pdf

figures/05_partition_vollenweider-1.pdf: figures/05_partition_vollenweider.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/05_partition_vollenweider.pdf

figures/06_k-1.pdf: figures/06_k.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/06_k.pdf

figures/07_cor_mat_hmap-1.pdf: figures/07_cor_mat_hmap.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/07_cor_mat_hmap.pdf

figures/08_maps-1.pdf: figures/08_maps.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm figures/08_maps.pdf

tables: manuscript/tables.pdf

manuscript/tables.pdf: tables/01_lake_characteristics_table.pdf tables/02_model_results_table.pdf
	pdftk $^ cat output manuscript/tables.pdf

tables/01_lake_characteristics_table.pdf: tables/01_lake_characteristics_table.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	
tables/02_model_results_table.pdf: tables/02_model_results_table.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	
scripts/table_1.csv: scripts/table_1.R
	Rscript scripts/table_1.R

.PHONY: all

all: 01_lake_characteristics_table.pdf 02_model_results_table.pdf

01_lake_characteristics_table.pdf: 01_lake_characteristics_table.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	
02_model_results_table.pdf: 02_model_results_table.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
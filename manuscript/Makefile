
.PHONY: all

all: figures.pdf appendix.pdf tables.pdf

draft.pdf: all draft.md
	pandoc draft.md -o draft.pdf

figures.pdf: figures.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	-pdftk figures.pdf cat 2-end output figures2.pdf
	-mv figures2.pdf figures.pdf
	
appendix.pdf: appendix.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	
tables.pdf: ../tables/01_lake_characteristics_table.pdf ../tables/02_model_results_table.pdf
	pdftk $^ cat output tables.pdf

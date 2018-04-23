analysis.pdf: analysis.Rmd 01_prepdata.R 03_viz.R 04_tables.R 99_utils.R
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	pdftk analysis.pdf cat 8-end output appendix.pdf
	pdftk analysis.pdf cat 1-7 output analysis2.pdf
	mv analysis2.pdf analysis.pdf

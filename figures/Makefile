.PHONY: all

all: iws_nws.pdf conny_metric_key.pdf 04_global_vollenweider-1.pdf 05_partition_vollenweider-1.pdf 06_k-1.pdf 07_cor_mat_hmap-1.pdf 08_maps-1.pdf
	mogrify -density 300 -format png *.pdf
	echo done

iws_nws.pdf: iws_nws.tex beamer-lake-fig/beamer-lake-fig.tex
	pdflatex iws_nws.tex
	pdfcrop iws_nws.pdf iws_nws.pdf
	
conny_metric_key.pdf: conny_metric_key.tex
	pdflatex conny_metric_key.tex
	pdfcrop conny_metric_key.pdf conny_metric_key_crop.pdf
	# convert -density 400 -trim conny_metric_key.pdf conny_metric_key.pdf
	#-rm conny_metric_key.pdf

04_global_vollenweider-1.pdf: 04_global_vollenweider.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm 04_global_vollenweider.pdf

05_partition_vollenweider-1.pdf: 05_partition_vollenweider.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm 05_partition_vollenweider.pdf

06_k-1.pdf: 06_k.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm 06_k.pdf

07_cor_mat_hmap-1.pdf: 07_cor_mat_hmap.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm 07_cor_mat_hmap.pdf

08_maps-1.pdf: 08_maps.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'pdf_document')"
	rm 08_maps.pdf



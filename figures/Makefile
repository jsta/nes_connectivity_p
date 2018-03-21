.PHONY: all

all: iws_nws.pdf conny_metric_key.pdf
	echo done

iws_nws.pdf: iws_nws.tex
	pdflatex iws_nws.tex
	pdfcrop iws_nws.pdf iws_nws.pdf
	
conny_metric_key.pdf: conny_metric_key.tex
	pdflatex conny_metric_key.tex
	pdfcrop conny_metric_key.pdf conny_metric_key_crop.pdf
	# convert -density 400 -trim conny_metric_key.pdf conny_metric_key.pdf
	#-rm conny_metric_key.pdf

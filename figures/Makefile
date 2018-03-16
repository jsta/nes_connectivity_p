.PHONY: all

all: iws_nws.png conny_metric_key.png
	echo done

iws_nws.png: iws_nws.tex
	pdflatex iws_nws.tex
	convert -density 400 iws_nws.pdf -trim iws_nws.png
	-rm iws_nws.pdf

conny_metric_key.png: conny_metric_key.tex
	pdflatex conny_metric_key.tex
	convert -density 400 conny_metric_key.pdf conny_metric_key.png
	#-rm conny_metric_key.pdf

crop/01_lake_characteristics_table.png: 01_lake_characteristics_table.pdf
	pdfcrop $< $(@F).pdf
	convert -density 300 -trim $(@F).pdf -quality 100 $@
	rm $(@F).pdf
	
crop/02_model_results_table.png: 02_model_results_table.pdf
	pdfcrop $< $(@F).pdf
	convert -density 300 -trim $(@F).pdf -quality 100 $@
	rm $(@F).pdf

.PHONY: png

png: crop/01_lake_characteristics_table.png crop/02_model_results_table.png

.PHONY: figures analysis data

data: 
	cd data && make all

draft.pdf: draft.md scripts/analysis.pdf
	pandoc draft.md -o draft.pdf

analysis: analysis.pdf

analysis.pdf: scripts/analysis.Rmd
	cd scripts && make analysis.pdf

figures:  
	cd figures && make all

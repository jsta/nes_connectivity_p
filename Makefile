.PHONY: figures analysis

draft.pdf: draft.md scripts/analysis.pdf
	pandoc draft.md -o draft.pdf

analysis: scripts/analysis.Rmd
	cd scripts && make analysis

figures:  
	cd figures && make all

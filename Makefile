draft.pdf: draft.md
	pandoc draft.md -o draft.pdf

analysis: scripts/analysis.Rmd
	cd scripts && make analysis

figures: analysis
	cd figures && make

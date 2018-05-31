.PHONY: figures analysis data tables

data: 
	cd data && make all

draft.pdf: draft.md scripts/analysis.pdf
	pandoc draft.md -o draft.pdf

all: analysis

analysis: analysis.pdf

analysis.pdf: scripts/analysis.Rmd scripts/table_1.csv
	cd scripts && make all

figures:  
	cd figures && make all

tables:
	cd tables && make all

scripts/table_1.csv: scripts/table_1.R
	Rscript scripts/table_1.R

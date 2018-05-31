.PHONY: figures data tables all

data: 
	cd data && make all

all: manuscript/draft.pdf figures tables scripts/table_1.csv

manuscript/draft.pdf:
	cd manuscript && make draft.pdf

figures:  
	cd figures && make all

tables:
	cd tables && make all

scripts/table_1.csv: scripts/table_1.R
	Rscript scripts/table_1.R

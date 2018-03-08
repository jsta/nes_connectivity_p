draft.pdf: draft.md
	pandoc draft.md -o draft.pdf

figures:
	cd figures && make


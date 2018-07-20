.PHONY: all

all:
	-mogrify -density 300 -set colorspace Gray *.pdf
	for file in *.pdf; do \
		mv $$file "gray_$$file"; \
	done
	-rm *~

clean:
	-rm *.pdf
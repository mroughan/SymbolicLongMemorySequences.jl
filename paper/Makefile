NAME = s5_annotated_bibliography
PDF = $(NAME).pdf
TEX = $(NAME).tex
AUX = $(NAME).aux $(NAME).log $(NAME).out

.PHONY: all clean distclean

all: $(PDF)

$(PDF): $(TEX)
	pdflatex -interaction=nonstopmode -halt-on-error $(TEX)

clean:
	rm -f $(AUX)

distclean: clean
	rm -f $(PDF)

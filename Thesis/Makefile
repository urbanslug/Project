all: main.pdf

main.pdf: main.tex main.bbl figures/*
	pdflatex -interaction nonstopmode main
	pdflatex -interaction nonstopmode main

main.bbl: library.bib
	pdflatex -interaction nonstopmode main
	bibtex main </dev/null

clean:
	rm -f *.aux *.toc *.pdf *.log *.bbl *.out *.lof *.lot *.blg

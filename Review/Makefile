all: review.pdf

review.pdf: main.tex review.bbl figures/*
	pdflatex  --jobname="review" -interaction nonstopmode main
	pdflatex  --jobname="review" -interaction nonstopmode main

review.bbl: ../library.bib
	pdflatex --jobname="review" -interaction nonstopmode main
	bibtex review </dev/null

clean:
	rm -f *.aux *.toc *.pdf *.log *.bbl *.out *.lof *.lot *.blg

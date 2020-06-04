all: review.pdf

review.pdf: review.tex review.bbl figures/*
	pdflatex -interaction nonstopmode review
	pdflatex -interaction nonstopmode review

review.bbl: library.bib
	pdflatex -interaction nonstopmode review
	bibtex review </dev/null

clean:
	rm -f *.aux *.toc *.pdf *.log *.bbl *.out *.lof *.lot *.blg

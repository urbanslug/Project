all: paper

paper: review.tex
	pdflatex --jobname="review" -shell-escape review.tex

clean:
	rm -f *.aux *.toc *.pdf *.log *.bbl *.out *.tex*

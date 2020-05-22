all: paper

paper: paper.tex
	pdflatex --jobname="review" -shell-escape paper.tex

clean:
	rm -f *.aux *.toc *.pdf *.log *.bbl *.out

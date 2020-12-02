latex: report.tex
	pdflatex report.tex > /dev/null;\
	bibtex report> /dev/null;\
	pdflatex report.tex> /dev/null;\
	pdflatex report.tex> /dev/null;

clean:
	$(RM) *.aux *.pdf *run *.blg *.toc *.xml *blx.bib *.bbl *.log

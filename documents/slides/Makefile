#file=Presentazione\ PACS.tex
file=Presentazione\ finale\ NAPDE.tex
all:
	pdflatex $(file) > /dev/null ;\
	bibtex $(file)> /dev/null;\
	pdflatex  $(file)> /dev/null;\
	pdflatex  $(file)> /dev/null;

clean:
	rm -f *.aux *.nav *.out *.snm *.toc *.pdf *.blg *.log *.bbl *.run.xml *.bib *.vrb
	

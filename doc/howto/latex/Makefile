NOTE          =

all doc: pdf


pdf:
	@pdflatex $(NOTE)
	@pdflatex $(NOTE)
	@makeindex $(NOTE).idx -s monindex.ist
	@pdflatex $(NOTE)
	@pdflatex $(NOTE)

clean:
	@\rm -f .tex .log *.dvi *.out *.log *.aux *.toc *.lof *.lot \
	 *.ilg *.idx *.bbl *.blg *.ind *.sym *.nsy 

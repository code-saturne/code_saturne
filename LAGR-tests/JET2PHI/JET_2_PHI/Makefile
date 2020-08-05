DOC = jet2phi

STYLE = ../STYLE/
REF = write-up

# Processor
LATEX = latex
PDF = pdflatex
BIBTEX = bibtex

# Error message
BADREF = "There were undefined references"

# Default values
MODE = --interaction batchmode
MODE =

pdf: $(STYLE)$(REF).tex
	make clean
	make include
	count=0 ;\
	while ( ($(PDF) $(MODE) $(STYLE)$(REF).tex ;\
		make bibtex ;\
		$(PDF) $(MODE) $(STYLE)$(REF).tex ;\
		$(PDF) $(MODE) $(STYLE)$(REF).tex ;\
		grep $(BADREF) $(REF).log > /dev/null) ) do \
		count=`expr $$count + 1` ;\
		if [ "$${count}" = "5" ] ;\
		  then break ;\
		fi ;\
	done
	make clean

include:
	for name in $(DOC) ; do \
		echo "\include{$$name}" >> infiles.tex ;\
	done

bibtex:
	for name in $(DOC) ; do \
		$(BIBTEX) $$name.aux ;\
	done

clean:
	rm -f *~
	rm -f *.log
	rm -f *.aux
	rm -f *.out
	rm -f *.toc
	rm -f *.cb*
	rm -f *.blg
	rm -f *.bbl
	rm -f infiles.tex

cleanall:
	rm -f *.pdf

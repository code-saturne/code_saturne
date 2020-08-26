DOC = lagr_point_source_dispersion

STYLE = ../STYLE/
REF = write-up

# Processor
LATEX = latex
PDF = pdflatex

# Error message
BADREF = "There were undefined references"

# Default values
MODE = --interaction batchmode
MODE =

pdf: $(STYLE)$(REF).tex
	make clean
	make input
	count=0 ;\
	while ( ($(PDF) $(MODE) $(STYLE)$(REF).tex ;\
		grep $(BADREF) $(REF).log > /dev/null) ) do \
		count=`expr $$count + 1` ;\
		if [ "$${count}" = "5" ] ;\
		  then break ;\
		fi ;\
	done
	make clean

input:
	for name in $(DOC) ; do \
		echo "\input{$$name}" >> infiles.tex ;\
	done

clean:
	rm -f *~
	rm -f *.log
	rm -f *.aux
	rm -f *.out
	rm -f *.toc
	rm -f *.cb*
	rm -f infiles.tex

cleanall:
	rm -f *.pdf

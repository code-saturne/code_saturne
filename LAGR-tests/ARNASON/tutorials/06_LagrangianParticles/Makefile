# Main LaTeX document

TEXMAIN = 06_LagrangianParticles.tex

TEXINPUTFILES = input

# Environment variables for LaTeX

PDFLATEX = pdflatex
BIBTEX = bibtex

TEXINPUTS = ../style/:
BIBINPUTS = ../style/:
BSTINPUTS = ../style/:

# LaTeX commands

PDFLATEXCMD = export TEXINPUTS="$(TEXINPUTS)"; "$(PDFLATEX)"
BIBTEXCMD = export BIBINPUTS="$(BIBINPUTS)"; "$(BIBTEX)"

# Generated documents (pdf and index files)

PDFMAIN = $(TEXMAIN:.tex=.pdf)

# Targets

all doc: pdf

pdf:
	@$(PDFLATEXCMD) $(TEXMAIN)
	@$(PDFLATEXCMD) $(TEXMAIN)
	@$(PDFLATEXCMD) $(TEXMAIN)

install:
	@mv *.pdf ..

clean:
	@\rm -f *.dvi *.out *.loe *.log *.aux *.toc *.lof *.los *.lot \
	 *.ilg *.idx *.bbl *.blg *.ind *.sym *.nsy

distclean:
	@\rm -f *.pdf *~
	make clean

#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main LaTeX document

TEXMAIN = user.tex

# Subdirectories (order is important)

SUBDIRS = graphics .

# Environment variables for LaTeX

TEXINPUTS = $(top_srcdir)/doc/style/:$(srcdir)//:$(builddir)//:
BIBINPUTS = $(top_srcdir)/doc/style/:
BSTINPUTS = $(top_srcdir)/doc/style/:
INDEXSTYLE = $(top_srcdir)/doc/style/:

# LaTeX commands

PDFLATEX = export TEXINPUTS="$(TEXINPUTS)"; pdflatex
BIBTEX = export BIBINPUTS="$(BIBINPUTS)"; bibtex
MAKEINDEX = export INDEXSTYLE="$(INDEXSTYLE)"; makeindex

# Generated documents (pdf and index files)

PDFMAIN = $(TEXMAIN:.tex=.pdf)
IDXMAIN = $(TEXMAIN:.tex=.idx)

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local:
	@$(PDFLATEX) $(TEXMAIN)
	@$(PDFLATEX) $(TEXMAIN)
	@$(MAKEINDEX) $(IDXMAIN) -s csindex.sty
	@$(PDFLATEX) $(TEXMAIN)
	@$(PDFLATEX) $(TEXMAIN)

clean-local:
	-rm -f *.out *.log *.aux *.toc *.lof *.lot \
	 *.ilg *.idx *.bbl *.blg *.ind *.sym *.nsy 

install-pdf-local:
	$(MKDIR_P) "$(DESTDIR)$(pdfdir)"
	cp -p $(PDFMAIN) $(DESTDIR)$(pdfdir)

distclean-local: clean
	-rm -f $(PDFMAIN)

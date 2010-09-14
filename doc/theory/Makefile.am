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
#   along with the Code_Saturne Kernel; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main LaTeX document

TEXMAIN = theory.tex

# Other LaTex document

TEXOTHERS = \
bilsc2.tex \
cfbase.tex \
cfener.tex \
cfmsvl.tex \
cfqdmv.tex \
cfxtcl.tex \
clptrg.tex \
clptur.tex \
clsyvt.tex \
codits.tex \
combbase.tex \
combcharb.tex \
combfuel.tex \
combgaz.tex \
condli.tex \
covofi.tex \
elbase.tex \
gradmc.tex \
gradrc.tex \
inimas.tex \
introd.tex \
itrmas.tex \
matrix.tex \
navsto.tex \
preduv.tex \
recvmc.tex \
resolp.tex \
thermo.tex \
summary.tex \
turbke.tex \
turrij.tex \
viscfa.tex \
visort.tex \
vissec.tex \
vortex.tex

TEXINPUTFILES = input

# Distributed files

EXTRA_DIST = $(TEXMAIN) $(TEXOTHERS) $(TEXINPUTFILES)

# Subdirectories (order is important)

SUBDIRS = graphics .

# Environment variables for LaTeX

TEXINPUTS = $(top_srcdir)/doc/style/:$(top_builddir)/doc/style/:$(srcdir)//:$(builddir)//:
BIBINPUTS = $(top_srcdir)/doc/style/:
BSTINPUTS = $(top_srcdir)/doc/style/:

# LaTeX commands

PDFLATEX = export TEXINPUTS="$(TEXINPUTS)"; pdflatex
BIBTEX = export BIBINPUTS="$(BIBINPUTS)"; bibtex

# Generated documents (pdf and index files)

PDFMAIN = $(TEXMAIN:.tex=.pdf)

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local:
	@$(PDFLATEX) $(TEXMAIN)
	@$(PDFLATEX) $(TEXMAIN)
	@$(PDFLATEX) $(TEXMAIN)

clean-local:
	-rm -f *.out *.log *.aux *.toc *.lof *.lot \
	 *.ilg *.idx *.bbl *.blg *.ind *.sym *.nsy 

install-pdf-local:
	$(MKDIR_P) "$(DESTDIR)$(pdfdir)"
	cp -p $(PDFMAIN) $(DESTDIR)$(pdfdir)

distclean-local: clean-local
	-rm -f $(PDFMAIN)

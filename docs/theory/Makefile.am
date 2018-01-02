## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

# Main LaTeX document

TEXMAIN = theory.tex

# Other LaTex documents

TEXOTHERS = \
algebr.tex \
apintr.tex \
bilsc2.tex \
biblio.tex \
bndcnd.tex \
cavitation.tex \
cfbase.tex \
cfener.tex \
cfmsvl.tex \
cfqdmv.tex \
cfxtcl.tex \
clptrg.tex \
clptur.tex \
clsyvt.tex \
codits.tex \
combustion_coal.tex \
combustion_gas.tex \
combustion_intro.tex \
combustion_thermo.tex \
comprs.tex \
condli.tex \
covofi.tex \
elbase.tex \
goveqn.tex \
gradmc.tex \
gradrc.tex \
ground_water.tex \
inimas.tex \
introd.tex \
itrmas.tex \
matrix.tex \
mesh_algo.tex \
mesh_quality.tex \
navstv.tex \
operat.tex \
predvv.tex \
resopv.tex \
spadis.tex \
timstp.tex \
turbul.tex \
summary.tex \
turbke.tex \
turrij.tex \
viscfa.tex \
visort.tex \
visecv.tex \
vortex.tex


TEXINPUTFILES = input

# Distributed files

EXTRA_DIST = $(TEXMAIN) $(TEXOTHERS) $(TEXINPUTFILES) \
graphics/clthermique.pdf \
graphics/loilog.pdf

# Generated pdf document

PDFMAIN = $(TEXMAIN:.tex=.pdf)

if HAVE_LATEX

# Subdirectories (order is important)

SUBDIRS = graphics .

# Environment variables for LaTeX

TEXINPUTS = $(top_srcdir)/docs/style/$(cs_tex_path_sep)$(top_builddir)/docs/style/$(cs_tex_path_sep)$(srcdir)//$(cs_tex_path_sep)$(builddir)//$(cs_tex_path_end)
BIBINPUTS = $(top_srcdir)/docs/style/
BSTINPUTS = $(top_srcdir)/docs/style/

# LaTeX commands

PDFLATEXCMD = export TEXINPUTS="$(TEXINPUTS)"; "$(PDFLATEX)"
BIBTEXCMD = export BIBINPUTS="$(BIBINPUTS)"; "$(BIBTEX)"
MAKEINDEXCMD = "$(MAKEINDEX)"

# generated index files

BIBMAIN = $(TEXMAIN:.tex=.aux)
GLOSSMAIN = $(TEXMAIN:.tex=.nlo)
INDEXMAIN = $(TEXMAIN:.tex=.nls)

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local:
	@$(PDFLATEXCMD) $(TEXMAIN)
	@$(BIBTEXCMD) $(BIBMAIN)
	@$(PDFLATEXCMD) $(TEXMAIN)
	@$(MAKEINDEXCMD) $(GLOSSMAIN) -s nomencl.ist -o $(INDEXMAIN)
	@$(PDFLATEXCMD) $(TEXMAIN)

clean-local:
	-rm -f *.out *.log *.aux *.toc *.lof *.lot \
	 *.ilg *.idx *.bbl *.blg *.ind *.sym *.nsy \
	 *.nlo *.nls

install-pdf-local:
	$(MKDIR_P) "$(DESTDIR)$(pdfdir)"
	$(INSTALL_DATA) $(PDFMAIN) $(DESTDIR)$(pdfdir)

distclean-local: clean-local
	-rm -f $(PDFMAIN)

# Generate the documentation at dist stage so that the user do not need
# to have LaTeX (and mandatory tools).

dist-hook: pdf
	$(INSTALL_DATA) $(PDFMAIN) $(distdir)

endif

# Install the documentation if available

install-data-local:
	if test -f $(srcdir)/$(PDFMAIN); then \
	  $(MKDIR_P) "$(DESTDIR)$(pdfdir)"; \
	  $(INSTALL_DATA) $(srcdir)/$(PDFMAIN) $(DESTDIR)$(pdfdir); \
	fi

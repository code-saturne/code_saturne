## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

TEXMAIN = refcard.tex

# Other LaTex documents

TEXOTHERS =

TEXINPUTFILES =

# Distributed files

EXTRA_DIST = $(TEXMAIN) $(TEXOTHERS) $(TEXINPUTFILES)

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

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local:
	@$(PDFLATEXCMD) $(TEXMAIN)
	@$(PDFLATEXCMD) $(TEXMAIN)
	@$(PDFLATEXCMD) $(TEXMAIN)

clean-local:
	-rm -f *.out *.log *.aux *.toc *.lof *.lot \
	 *.ilg *.idx *.bbl *.blg *.ind *.sym *.nsy 

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

uninstall-local:
	@$(NORMAL_UNINSTALL)
	files=$(PDFMAIN); \
	test -d "$(DESTDIR)$(pdfdir)" || exit 0; \
	test -n "$$files" || exit 0; \
	echo " ( cd '$(DESTDIR)$(pdfdir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(pdfdir)" && rm -f $$files
	rmdir $(pdfdir) || exit 0

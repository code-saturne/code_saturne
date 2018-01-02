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

# Xfig files to process

FIGFILES =

# Generated documents (pdf files)

PDF_FIGFILES = $(FIGFILES:.fig=.pdf)

# New suffixes and targets

SUFFIXES = .fig .pdf

# On Windows, fig2dev may be broken and generate a file named "-" instead
# of the target. Thus, we try to renamed such files if they exist.

.fig.pdf:
	"$(FIG2DEV)" -L pdf $< > $@
	@if test -f - ; then \mv -f - $@ ; fi

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local: $(PDF_FIGFILES)

clean-local:
	-rm -f $(PDF_FIGFILES)

distclean-local: clean-local

## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

LATEX_SUBDIRS = refcard theory user

DOXYGEN_SUBDIRS = doxygen

SUBDIRS = $(LATEX_SUBDIRS) $(DOXYGEN_SUBDIRS)

EXTRA_DIST = \
style/csbiblio.bib \
style/csdoc.cls \
style/csshortdoc.cls \
style/nomencl.cfg \
style/csdoc.sty \
style/csindex.sty \
style/csmacros.sty \
style/csvers.tex.in \
style/cs_logo.pdf \
style/cs_logo_h.pdf \
style/cs_logo_v.pdf \
style/logoedf.png

# Macro-rules so as to compile/install both PDF and HTML outputs

doc: pdf html
install-doc: install-pdf install-html

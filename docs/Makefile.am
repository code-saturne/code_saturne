## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

LATEX_SUBDIRS = refcard theory user install developer studymanager

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
style/cs_logo_black.pdf \
style/cs_logo_flux.pdf \
style/cs_logo_wire_black.pdf \
style/cs_logo_bubble.pdf \
style/cs_logo_wave.pdf \
style/logoedf.png

# Macro-rules so as to compile/install both PDF and HTML outputs

doc: pdf html
install-doc: install-pdf install-html

## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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

if HAVE_SALOME_KERNEL
if HAVE_SPHINX
SALOME = salome
endif
endif

SUBDIRS = refcard theory tutorial user install developper doxygen autovnv $(SALOME)

DIST_SUBDIRS = $(SUBDIRS) salome

EXTRA_DIST = \
style/csdoc.sty \
style/csindex.sty \
style/csmacros.sty \
style/csvers.tex.in \
style/logocs.png \
style/logoedf.png

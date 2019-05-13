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

include $(top_srcdir)/config/salome_module.mk

MODULE_NAME = CFDSTUDY

RSTFILES = \
index.rst \
intro.rst \
tutorial/coupling.rst \
tutorial/debugging.rst \
tutorial/saturne.rst \
tutorials.rst

SPHINXOPTS = -d .doctrees $(srcdir)

EXTRA_DIST = images tutorial conf.py index.rst intro.rst tutorials.rst

html-local:
	$(SPHINXBUILD) -b html $(SPHINXOPTS) html

latex:
	$(SPHINXBUILD) -b latex $(SPHINXOPTS) latex

clean-local:
	-rm -rf html latex .doctrees

install-data-local:
	$(MKDIR_P) $(DESTDIR)$(salomedocdir)
	if test -d "html"; then cp -rf html/* $(DESTDIR)$(salomedocdir); fi;

uninstall-local:
	-rm -rf $(DESTDIR)$(salomedocdir)/*

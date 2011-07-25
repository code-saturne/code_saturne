#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2011 EDF S.A., France
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

include $(top_srcdir)/config/salome_module.mk

MODULE_NAME = CFDSTUDY

RSTFILES = \
index.rst \
intro.rst \
tutorial/saturne.rst \
tutorials.rst

SPHINXBUILD = sphinx-build
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

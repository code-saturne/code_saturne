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

MODULE_NAME = FSI

edit_python = sed \
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g'

FSI_ASTER_module.py: Makefile
	rm -f $@ $@.tmp
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_python) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

FSI_ASTER_module.py: $(srcdir)/FSI_ASTER_module.py.in 

nodist_salomepython_PYTHON = FSI_ASTER_module.py
salomepython_PYTHON = FSI_ASTER_component.py

EXTRA_DIST = FSI_ASTER_module.py.in

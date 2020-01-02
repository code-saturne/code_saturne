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

# Turn off the Automake warning about the makefile portability as we need
# GNUmake rules in order to account for the PyQt user interface files target.

AUTOMAKE_OPTIONS = -Wno-portability

include $(top_srcdir)/config/salome_module.mk

MODULE_NAME = cfdstudy

# Python files
dist_salomepython_PYTHON = \
CFDSTUDYOTURNS_DistantLauncher.py \
CFDSTUDYOTURNS_LocalLauncher.py \
CFDSTUDYOTURNS_StudyInterface.py

#resources
salomeres_DATA = $(QMFILES) ${ICONS}

# Distribute icons directory, resources files and user interface files

if HAVE_QT5
  PYUIC=$(PYUIC5)
else
  PYUIC=$(PYUIC4)
endif

ui_%.py: %.ui
	$(PYUIC) $< -x -o $@

install-exec-local: $(PYUI_FILES)
	for f in $^ ; do\
		$(INSTALL) $$f "$(DESTDIR)$(salomepythondir)/$$f" ; \
	done ;

## meta object implementation files generation (moc)
#%_moc.cxx: %.h
#	$(MOC) $< -o $@

## translation (*.qm) files generation (lrelease)
%.qm: %.ts
	$(LRELEASE) $< -qm $@

mostlyclean-local:
	rm -f @builddir@/*_moc.cxx
	rm -f @builddir@/*.qm

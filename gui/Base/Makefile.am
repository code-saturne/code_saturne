## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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

# Python files

PYFILES = \
BrowserView.py \
CommandMgrDialogView.py \
Common.py \
IdView.py \
MainView.py \
QtPage.py \
Toolbox.py \
XMLengine.py \
XMLinitialize.py \
XMLmodel.py \
XMLvariables.py

# PyQt User Interface files

UIFILES = \
BrowserForm.ui \
CommandMgrLinesDisplayedDialogForm.ui \
CommandMgrDialogForm.ui \
IdForm.ui \
MainForm.ui

# PyQt Resources files

QRCFILES = \
resource_base.qrc

# Distribute icons directory, resources files and user interface files

EXTRA_DIST = \
icons \
$(UIFILES) \
$(QRCFILES)

# Generated files

FORMS = $(UIFILES:.ui=.py)
RESOURCES = $(QRCFILES:.qrc=_rc.py)

# Specific rules

PYRCC4VER = `$(PYTHON) -V 2>&1 | cut -c 8`

SUFFIXES = .ui .qrc _rc.py

.ui.py:
	$(PYUIC4) -o $@ $<

.qrc_rc.py:
	$(PYRCC4) -py$(PYRCC4VER) -o $@ $<

all-local: $(FORMS) $(RESOURCES)

clean-local:
	-rm -f $(FORMS) $(RESOURCES)

# Install Python library

basedir = $(pkgpythondir)/Base
dist_base_PYTHON = \
__init__.py \
$(PYFILES)

nodist_base_PYTHON = \
$(FORMS) \
$(RESOURCES)

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

if HAVE_QT5
  PYUIC=$(PYUIC5)
else
  PYUIC=$(PYUIC4)
endif

if HAVE_QT5
  PYRCC = $(PYRCC5)
else
  PYRCC = $(PYRCC4) -py$(PYRCC4VER)
endif

# Python files

PYFILES = \
BrowserView.py \
CommandMgrDialogView.py \
CompletionTextEditor.py \
MainView.py \
QFileEditor.py \
QtCase.py \
QtPage.py \
Toolbox.py

if HAVE_QT5
  PYFILES += qt5/QtCore.py qt5/QtGui.py qt5/QtWidgets.py
else
  PYFILES += qt4/QtCore.py qt4/QtGui.py qt4/QtWidgets.py
endif

# PyQt User Interface files

UIFILES = \
BrowserForm.ui \
CommandMgrLinesDisplayedDialogForm.ui \
CommandMgrDialogForm.ui \
MainForm.ui \
NewCaseDialogForm.ui

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

PYRCC4VER = `$(PYTHON) --version 2>&1 | cut -c 8`

SUFFIXES = .ui .qrc _rc.py

.ui.py:
	$(PYUIC) -o $@ $<

.qrc_rc.py:
	$(PYRCC) -o $@ $<

all-local: $(FORMS) $(RESOURCES)

clean-local:
	-rm -f $(FORMS) $(RESOURCES)

# Install Python library

basedir = $(pkgpythondir)/Base

dist_base_PYTHON = \
__init__.py \
$(PYFILES) \
icons/22x22/MONO-bulle-HD.png \
icons/22x22/logoneptune.png

nodist_base_PYTHON = \
$(FORMS) \
$(RESOURCES)

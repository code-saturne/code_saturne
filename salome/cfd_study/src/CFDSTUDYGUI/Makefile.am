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

# Turn off the Automake warning about the makefile portability as we need
# GNUmake rules in order to account for the PyQt user interface files target.

AUTOMAKE_OPTIONS = -Wno-portability

include $(top_srcdir)/config/salome_module.mk

MODULE_NAME = cfdstudy

# Python files
dist_salomepython_PYTHON = \
CFDSTUDYGUI.py \
CFDSTUDYGUI_ActionsHandler.py \
CFDSTUDYGUI_Agents.py \
CFDSTUDYGUI_CommandMgr.py \
CFDSTUDYGUI_Commons.py \
CFDSTUDYGUI_DataModel.py \
CFDSTUDYGUI_DesktopMgr.py \
CFDSTUDYGUI_DialogCollector.py \
CFDSTUDYGUI_SolverGUI.py \
CFDSTUDYGUI_Message.py \
CFDSTUDYGUI_Management.py

#resources
salomeres_DATA = $(QMFILES) ${ICONS}

nodist_salomeres_DATA = \
CFDSTUDY_images.qm \
CFDSTUDY_msg_en.qm \
CFDSTUDY_msg_fr.qm

# PyQt User Interface files

UIFILES = \
CFDSTUDYGUI_QProcessDialog.ui \
ECSConversionDialog.ui \
GUIActivationDialog.ui \
InfoDialog.ui \
SetTreeLocationDialog.ui

PYUI_FILES = \
ui_CFDSTUDYGUI_QProcessDialog.py \
ui_ECSConversionDialog.py \
ui_GUIActivationDialog.py \
ui_InfoDialog.py \
ui_SetTreeLocationDialog.py

# Distribute icons directory, resources files and user interface files

EXTRA_DIST = \
$(UIFILES) \
CFDSTUDY_images.ts \
CFDSTUDY_msg_en.ts \
CFDSTUDY_msg_fr.ts

if HAVE_QT5
  PYUIC=$(PYUIC5)
endif

ui_%.py: %.ui
	$(SALOMEENVCMD); \
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
	$(SALOMEENVCMD); \
	$(LRELEASE) $< -qm $@

mostlyclean-local:
	rm -f @builddir@/*_moc.cxx
	rm -f @builddir@/*.qm

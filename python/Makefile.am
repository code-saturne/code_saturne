## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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
if HAVE_GUI
ENABLE_GUI = --with-gui

# QT

if HAVE_PYSIDE6
USE_QT = --use-qt=pyside6
if HAVE_PYUIC
  CMD_UIC = --pyuic $(PYSIDE6UIC)
else
  CMD_UIC = --uic $(UIC)
endif
if HAVE_PYRCC
  CMD_CC = --pyrcc $(PYSIDE6RCC)
else
  CMD_RCC = --rcc $(RCC)
endif
endif

if HAVE_PYQT6
USE_QT = --use-qt=pyqt6
CMD_UIC = --pyuic $(PYUIC6)
CMD_RCC = --rcc $(RCC)
endif

if HAVE_PYQT5
USE_QT = --use-qt=pyqt5
CMD_UIC = --pyuic $(PYUIC5)
CMD_RCC = --pyrcc $(PYRCC5)
endif

else
ENABLE_GUI = --without-gui
endif

# Copy all files in pythondir
EXTRA_DIST = ${srcdir}/*

clean-local:
	if test -d build; then \
		echo "rm -rf build"; \
		rm -rf build; \
	fi

all-local: clean-local
	$(PYTHON) ${abs_top_srcdir}/python/setup.py build \
		--build-base ${abs_top_builddir}/python/build \
		$(ENABLE_GUI) $(USE_QT) $(CMD_UIC) $(CMD_RCC) \
		--version $(VERSION) \
		--verbose

install-exec-local:
	$(PYTHON) ${abs_top_srcdir}/python/setup.py install \
		--install-lib $(DESTDIR)/$(pythondir) \
		--single-version-externally-managed \
		--record $(DESTDIR)/$(pkgpythondir)/install_files.txt \
		--version $(VERSION) \
		--verbose

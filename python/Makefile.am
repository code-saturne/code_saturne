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
ENV_GUI = CS_ENABLE_GUI=1

# QT

if HAVE_PYSIDE6
ENV_QT = CS_USE_QT=pyside6
if HAVE_PYUIC
  ENV_UIC = CS_PYUIC="$(PYSIDE6UIC)"
else
  ENV_UIC = CS_UIC="$(UIC)"
endif
if HAVE_PYRCC
  ENV_RCC = CS_PYRCC="$(PYSIDE6RCC)"
else
  ENV_RCC = CS_RCC="$(RCC)"
endif
endif

if HAVE_PYQT6
ENV_QT = CS_USE_QT=pyqt6
ENV_UIC = CS_PYUIC="$(PYUIC6)"
ENV_RCC = CS_RCC="$(RCC)"
endif

if HAVE_PYQT5
ENV_QT = CS_USE_QT=pyqt5
ENV_UIC = CS_PYUIC="$(PYUIC5)"
ENV_RCC = CS_PYRCC="$(PYRCC5)"
endif

else
ENV_GUI = CS_ENABLE_GUI=0
endif

# Copy all files in pythondir
EXTRA_DIST = ${srcdir}/*

# Install data

imagesdir = ${pkgdatadir}/images

images_DATA = \
${abs_top_srcdir}/python/code_saturne/data/icons/splashscreen.png \
${abs_top_srcdir}/python/code_saturne/data/icons/logo_salome_cfd.png


themesdir = ${pkgdatadir}/themes

themes_DATA = \
${abs_top_srcdir}/python/code_saturne/data/themes/*.qss

clean-local:
	if test -d build; then \
		echo "rm -rf build"; \
		rm -rf build; \
	fi
	if test -d dist; then \
		echo "rm -rf dist"; \
		rm -rf dist; \
	fi
	if test "$(abs_top_builddir)" != "$(abs_top_srcdir)" && test -d $(srcdir)/build; then \
		echo "rm -rf $(srcdir)/build"; \
		rm -rf $(srcdir)/build; \
	fi

all-local: clean-local
	$(ENV_GUI) $(ENV_QT) $(ENV_UIC) $(ENV_RCC) CS_VERSION="$(VERSION)" \
	CS_PYBUILD_DIR="${abs_top_builddir}/python/build" \
	$(PYTHON) -m pip wheel --no-deps --no-build-isolation \
		--wheel-dir ${abs_top_builddir}/python/dist \
		${abs_top_srcdir}/python

install-exec-local:
	$(PYTHON) -m pip install --no-deps --no-build-isolation --no-warn-script-location --upgrade \
		--target "$(DESTDIR)$(pythondir)" \
		${abs_top_builddir}/python/dist/*.whl

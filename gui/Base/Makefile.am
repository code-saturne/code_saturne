#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009 EDF S.A., France
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

all-local:
	@$(PYTHON) $(top_srcdir)/gui/bootstrap.py --build $(srcdir)

clean-local:
	@$(PYTHON) $(top_srcdir)/gui/bootstrap.py --clean $(srcdir)

# Install Python library

basedir = $(pkgpythondir)/Base
base_PYTHON = \
__init__.py \
resource_base_rc.py \
BrowserForm.py \
BrowserView.py \
CommandLine.py \
Common.py \
IdForm.py \
IdView.py \
MainForm.py \
MainView.py \
QtPage.py \
Toolbox.py \
XMLengine.py \
XMLinitialize.py \
XMLmodel.py \
XMLvariables.py

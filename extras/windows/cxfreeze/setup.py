#!/usr/bin/env python

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

import sys, os
from cx_Freeze import setup, Executable

#-------------------------------------------------------------------------------

# Preparing environment
# ---------------------

# Module search path
path = sys.path + ["bin", "lib/python" + sys.version[:3] + "/site-packages/code_saturne", "lib/python" + sys.version[:3] + "/site-packages"]
path += [sys.prefix + "/Lib/site-packages"]

# Specific modules to be included
includes = ["sip"]

# Specific modules to be excluded
m_script = ["cs_user_scripts"]
m_studymanager = ["matplotlib", "vtk", "numpy"]
m_neptune = ["nc_package", "core.XMLinitialize", "core.MainView"]
m_syrthes = ["syrthes"]
m_salome = ["Pages.SalomeHandler"]
m_win32 = ["win32api", "win32con", "win32pipe"]
excludes = m_studymanager + m_script + m_neptune + m_syrthes + m_salome + m_win32

# Specific packages
packages = []

# Copy of some mandatory files or directories
includefiles = []
if sys.platform.startswith("linux"):
    includefiles += [(r"/usr/lib/qt4/translations", \
                       "translations")]
elif sys.platform.startswith("win"):
    includefiles += [(r"%s\Lib\site-packages\PyQt4\translations" % sys.prefix, \
                       "translations")]
else:
    pass

# Possible inclusion of additional libraries
binpathincludes = []
if sys.platform.startswith("linux"):
    binpathincludes += ["/usr/lib"]

# Build the options dictionnary
options = {"path": path,
           "includes": includes,
           "excludes": excludes,
           "packages": packages,
           "include_files": includefiles,
           "bin_path_includes": binpathincludes}

#-------------------------------------------------------------------------------

# Preparing targets
# -----------------

# Windows (win32) does not support a single executable being both a command
# line script and a graphical user interface, as code_saturne is.

# One possible trick is to generate two executables, one for the command line
# uses (ocde_saturne.com) and one to launch the graphical interface
# (code_saturne.exe)

# When calling "code_saturne" in cmd.exe or PowerShell, the "code_saturne.com"
# will first be chosen due to order rules of win32, leaving us the capability
# of launching the graphical interface through "code_saturne.exe".

# If not using this trick, "code_saturne.exe" cannot use stdout/stderr when
# run as a command line tool.

base = None
if sys.platform == "win32":
    base_gui = "Win32GUI"
    base_cli = None

target_gui = Executable(script = "bin/code_saturne",
                        targetName = "code_saturne.exe",
                        base = base_gui,
                        compress = True,
                        icon = "bin/code_saturne.ico")

target_cli = Executable(script = "bin/code_saturne",
                        targetName = "code_saturne.com",
                        base = base_cli,
                        compress = True,
                        icon = None)

#-------------------------------------------------------------------------------

# Creating the setup
# ------------------

setup(name = "Code_Saturne",
      version = "6.1",
      description = "General purpose CFD software",
      author = "EDF",
      options = {"build_exe": options},
      executables = [target_gui, target_cli])

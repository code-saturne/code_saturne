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

# This file is adapted from make_common_starter.am, generated by YACSGEN.

# Standard directory for installation
salomeincludedir   = ${includedir}/salome
salomelibdir       = ${libdir}/salome
salomebindir       = ${bindir}/salome
salomepythondir    = ${pythondir}/salome

# Directory for installing idl files
salomeidldir       = ${datarootdir}/idl/salome

# Directory for installing resource files
salomeresdir       = ${datarootdir}/salome/resources/${MODULE_NAME}

# Shared modules installation directory
sharedpkgpythondir = ${pythondir}/${MODULE_NAME}/shared_modules

# Documentation directory
salomedocdir       = ${datarootdir}/doc/salome/gui/${MODULE_NAME}

IDL_INCLUDES = $(SALOME_KERNEL_IDL) $(SALOME_GUI_IDL)
KERNEL_LIBS = $(SALOME_KERNEL_LDFLAGS) -lSalomeContainer -lOpUtil -lSalomeDSCContainer -lSalomeDSCSuperv -lSalomeDatastream -lSalomeDSCSupervBasic -lCalciumC
KERNEL_INCLUDES= $(SALOME_KERNEL_CPPFLAGS) $(OMNIORB_INCLUDES)

SALOME_LIBS= $(KERNEL_LIBS)
SALOME_IDL_LIBS= $(SALOME_KERNEL_LDFLAGS) -lSalomeIDLKernel
SALOME_INCLUDES= $(KERNEL_INCLUDES)


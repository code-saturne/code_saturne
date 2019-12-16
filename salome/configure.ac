# Process this file with autoconf to produce a configure script

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

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([cs_licence_c_comment],
[/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/])

# Note: proper determination of the version string requires the parent NEWS file
#       to be maintained properly (as it should be anyways).

m4_define([salome_cfd_version_string],
          [m4_normalize(esyscmd([../build-aux/cs_version.py]))])

AC_INIT([salome_cfd],[salome_cfd_version_string],[saturne-support@edf.fr],[],[https://code-saturne.org])
AC_CONFIG_SRCDIR([cfd_study/src/CFDSTUDY/CFDSTUDY.py])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([salome_cfd_config.h])

# Initialize automake with the following options:
# - foreign, so as to handle the absence of ChangeLog (automatically generated)
# - tar-pax, so as to handle long lines (> 99 characters) in tar archives
# Warnings can be activated at bootstrap with 'autoreconf -vi --warnings=all'
AM_INIT_AUTOMAKE([foreign tar-pax subdir-objects])

# Enable maintainer mode by default for a developer checkout
AS_IF([test -d ${srcdir}/../.git],
      [AM_MAINTAINER_MODE([enable])],
      [AM_MAINTAINER_MODE([disable])])

# Enable silent rules with "./configure --enable-silent-rules" or "make V=0"
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

#------------------------------------------------------------------------------
# Checks for Python support.
#------------------------------------------------------------------------------

#AC_ARG_VAR([PYTHON], [the Python interpreter])
AM_PATH_PYTHON(["3.4"])

# Check for compilers
#----------------------------------

AC_PROG_CC
AC_PROG_CXX

AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_CC_C_O
AC_PROG_SED

if test "x$CC" = "x" ; then
  AC_MSG_FAILURE([cannot find C compiler])
fi

# Needed to use Microsoft archiver lib.exe
# It copies the ar-lib script, similar to the compile script
#m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Debug or production compilation mode (debug by default) ?

AC_ARG_ENABLE(debug,
  [AS_HELP_STRING([--enable-debug], [enable debugging (reduces optimization)])],
  [
    case "${enableval}" in
      yes) debug=yes ;;
      no)  debug=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac
  ],
  [ debug=no ]
)
AC_SUBST(debug)

# Optionally install the code with relocatable features

AC_ARG_ENABLE(relocatable,
  [AS_HELP_STRING([--enable-relocatable], [enable relocatable installation])],
  [
    case "${enableval}" in
      yes) relocatable=yes ;;
      no)  relocatable=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-relocatable]) ;;
    esac
  ],
  [ relocatable=no ]
)
AC_SUBST(relocatable)
AM_CONDITIONAL(HAVE_RELOCATABLE, [test "${relocatable}" = yes])
if test "${relocatable}" = yes ; then
  AC_DEFINE([HAVE_RELOCATABLE], 1, [Relocatable installation])
fi

# Preprocessor flags for debugging purposes

if test "x$debug" = "xyes"; then
  CPPFLAGS="${CPPFLAGS} -DDEBUG"
else
  CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi

# We may only turn on processing for libtool now that the basic compiler
# and linker flags are set (to avoid issues with linkers with different
# modes such as 32 and 64 bit which may be modified by compiler or
# linker flags).

# Initialize libtool.

m4_pushdef([salome_cfd_lt_init_args], [disable-static, no-pic])
LT_INIT([disable-static])

#------------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------------

AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_PID_T

AC_HEADER_TIME
AC_HEADER_STDBOOL

#------------------------------------------------------------------------------
# Check for salome features.
#------------------------------------------------------------------------------

SALOME_CFD_AC_SALOME_ENV
SALOME_CFD_AC_TEST_SALOME

#------------------------------------------------------------------------------
# Checks for PyQt support.
#------------------------------------------------------------------------------

# Test for PyQt5
salome_cfd_have_gui_qt5=no

if test "x$salome_cfd_have_gui" = "xyes" ; then

  AC_MSG_CHECKING([for PyQt5 version >= 5.0])

  prog="[import sys, string
try: import PyQt5
except ImportError: sys.exit(1)
from PyQt5.QtCore import *
if list(map(int, QT_VERSION_STR.split('.'))) < [5,0,0]: sys.exit(1)
if list(map(int, PYQT_VERSION_STR.split('.'))) < [5,0,0]: sys.exit(1)
sys.exit(0)]"

  ${PYTHON} -c "${prog}"
  retval=$?

  if test $retval -ne 0 ; then
    AC_MSG_RESULT([no])
    salome_cfd_have_gui=no
  else
    AC_MSG_RESULT([yes])
    salome_cfd_have_gui=yes
    salome_cfd_have_gui_qt5=yes
  fi

fi

# Test for PyQt5 developer tools
if test "x$salome_cfd_have_gui_qt5" = "xyes" ; then

  # Try to find PyQt5 tools through Python interpreter
  prog="[import os, sys
from PyQt5.QtCore import QCoreApplication
app = QCoreApplication([])
path = app.applicationDirPath()
sys.stdout.write(path)]"

  PYPATH=`${PYTHON} -c "${prog}"`

  AC_ARG_VAR([PYUIC5], [PyQt5 user interfaces compiler])
  AC_ARG_VAR([PYRCC5], [PyQt5 resources compiler])
  AC_ARG_VAR([LRELEASE], [Qt message files generator])

  if test "x$PYUIC5" = "x" ; then
    AC_PATH_PROG([PYUIC5], [pyuic5], [$CMDPYUIC5], [${PYPATH}:$PATH])
  fi
  if test "x$PYRCC5" = "x" ; then
    AC_PATH_PROG([PYRCC5], [pyrcc5], [],           [${PYPATH}:$PATH])
  fi

  if test "x$LRELEASE" = "x" ; then
    AC_PATH_PROGS([LRELEASE], [lrelease-qt5 lrelease])
  fi

  if test "x$PYUIC5" = "x" -o "x$PYRCC5" = "x"  -o "x$LRELEASE" = "x"; then
    AC_MSG_FAILURE([cannot find PyQt5 dev tools, Graphical User Interface cannot be installed])
    salome_cfd_have_gui=no
  fi

fi

#------------------------------------------------------------------------------
# Checks for Qt tools.
#------------------------------------------------------------------------------

AM_CONDITIONAL(HAVE_GUI, [test "x$salome_cfd_have_gui" = "xyes"])
AC_SUBST(salome_cfd_have_gui)
AM_CONDITIONAL(HAVE_QT5, [test "x$salome_cfd_have_gui_qt5" = "xyes"])
AC_SUBST(salome_cfd_have_gui_qt5)

#------------------------------------------------------------------------------
# Checks for documantation tools
#------------------------------------------------------------------------------

salome_cfd_have_sphinx=yes

AC_ARG_VAR([SPHINXBUILD], [Sphinx documentation tool])

dnl where is sphinx ?
AC_PATH_PROG(SPHINXBUILD, sphinx-build)
if test "x$SPHINXBUILD" = "x"; then
  AC_MSG_WARN(sphinx-build not found)
  salome_cfd_have_sphinx=no
fi

AM_CONDITIONAL(HAVE_SPHINX, [test $salome_cfd_have_sphinx = yes])
AC_SUBST(salome_cfd_have_sphinx)

#------------------------------------------------------------------------------
# Pass compiler options to automake file
#------------------------------------------------------------------------------

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------

# Main configuration

AC_CONFIG_FILES([Makefile
                 cfd_study/Makefile
                 cfd_study/doc/Makefile
                 cfd_study/idl/Makefile
                 cfd_study/resources/Makefile
                 cfd_study/src/Makefile
                 cfd_study/src/CFDSTUDY/Makefile
                 cfd_study/src/CFDSTUDYGUI/Makefile
                 cfd_study/src/CFDSTUDYOTURNS/Makefile])

#------------------------------------------------------------------------------

AC_OUTPUT

#------------------------------------------------------------------------------

echo
echo "Configuration options:"
echo " use debugging code: "$debug""
echo " SALOME KERNEL support: "$salome_cfd_have_kernel""
echo " SALOME GUI support: "$salome_cfd_have_gui""
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the salome CFD extensions"
echo


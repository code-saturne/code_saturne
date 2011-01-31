dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2011 EDF S.A., France
dnl
dnl   The Code_Saturne Kernel is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne Kernel is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne Preprocessor; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_TEST_SPHINX
#------------------
# modifies or sets cs_have_sphinx, SPHINX depending on libraries found

AC_DEFUN([CS_AC_TEST_SPHINX],[

AC_CHECKING(for sphinx doc generator)

cs_have_sphinx=yes
dnl where is sphinx ?
AC_PATH_PROG(SPHINX, sphinx-build) 
if test "x$SPHINX" = "x"; then
  AC_MSG_WARN(sphinx not found)
  cs_have_sphinx=no
fi

AM_CONDITIONAL(HAVE_SPHINX, [test $cs_have_sphinx = yes])

])dnl

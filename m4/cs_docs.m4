dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2011 EDF S.A.
dnl
dnl This program is free software; you can redistribute it and/or modify it under
dnl the terms of the GNU General Public License as published by the Free Software
dnl Foundation; either version 2 of the License, or (at your option) any later
dnl version.
dnl
dnl This program is distributed in the hope that it will be useful, but WITHOUT
dnl ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
dnl FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
dnl details.
dnl
dnl You should have received a copy of the GNU General Public License along with
dnl this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
dnl Street, Fifth Floor, Boston, MA 02110-1301, USA.
dnl
dnl--------------------------------------------------------------------------------

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

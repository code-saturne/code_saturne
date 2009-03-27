dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2009 EDF S.A., France
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

# CS_AC_TEST_PREPRO(Minimal Release string, [Maximal Release string])
#--------------------------------------------------------------------
# Check for the Preprocessor executable ; defines ecs_prefix

AC_DEFUN([CS_AC_TEST_PREPRO], [

AC_ARG_WITH(prepro, [AS_HELP_STRING([--with-prepro=PATH], [specify prefix directory for the Preprocessor])])

AC_CHECK_PROG(ecs_prefix, cs_preprocess, $with_prepro, , $with_prepro/bin, )
AC_SUBST(ecs_prefix)

])dnl


# CS_AC_TEST_GUI(Minimal Release string, [Maximal Release string])
#-----------------------------------------------------------------
# Check for the GUI executable ; defines ics_prefix

AC_DEFUN([CS_AC_TEST_GUI], [

AC_ARG_WITH(gui, [AS_HELP_STRING([--with-gui=PATH], [specify prefix directory for the GUI])])

AC_CHECK_PROG(ics_prefix, ics, $with_gui, , $with_gui/bin, )
AC_SUBST(ics_prefix)

])dnl


# CS_AC_TEST_SYRCS(Minimal Release string, [Maximal Release string])
#-------------------------------------------------------------------
# Check for the SYR_CS library ; defines syrcs_prefix

AC_DEFUN([CS_AC_TEST_SYRCS], [

AC_ARG_WITH(syrcs, [AS_HELP_STRING([--with-syrcs=PATH], [specify prefix directory for SYR_CS library])])

if test "x$with_syrcs" != "x" ; then
  syrcs_bindir="$with_syrcs/bin"
fi

if test -f $syrcs_bindir/syr_cs_profile ; then
  syrcs_prefix=$with_syrcs
fi
AC_SUBST(syrcs_prefix)

])dnl


# CS_AC_TEST_PYTHON()
#--------------------
# Check specific Python executable (if needed); defines python_prefix

AC_DEFUN([CS_AC_TEST_PYTHON], [

AC_ARG_WITH(python, [AS_HELP_STRING([--with-python=PATH], [specify prefix directory for Python])])

AC_CHECK_PROG(python_prefix, python, $with_python, , $with_python/bin, )
AC_SUBST(python_prefix)

])dnl

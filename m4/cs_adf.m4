dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2010 EDF S.A., France
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

# CS_AC_TEST_ADF
#---------------
# modifies or sets cs_have_adf, ADF_LDFLAGS, and ADF_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_ADF], [

cs_have_adf=no

AC_ARG_WITH(adf,
            [AS_HELP_STRING([--with-adf=DIR],
                            [specify prefix directory for ADF])],
            [if test "x$withval" = "x"; then
               with_adf=yes
             fi],
            [with_adf=check])

AC_ARG_WITH(adf-lib,
            [AS_HELP_STRING([--with-adf-lib=DIR],
                            [specify directory for ADF library])],
            [if test "x$with_adf" = "xcheck"; then
               with_adf=yes
             fi
             ADF_LDFLAGS="-L$with_adf_lib"
             # Add the libdir to the runpath as ADF is not libtoolized
             ADFRUNPATH="-R$with_adf_lib"],
            [if test "x$with_adf" != "xno" -a "x$with_adf" != "xyes" \
	          -a "x$with_adf" != "xcheck"; then
               ADF_LDFLAGS="-L$with_adf/lib"
               # Add the libdir to the runpath as ADF is not libtoolized
               ADFRUNPATH="-R$with_adf/lib"
             fi])


if test "x$with_adf" != "xno" ; then

  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  ADF_LIBS="-ladf"
  LDFLAGS="${LDFLAGS} ${ADF_LDFLAGS}"
  LIBS="${LIBS} ${ADF_LIBS}"

  AC_CHECK_LIB(adf, ADF_Database_Open, 
               [ AC_DEFINE([HAVE_ADF], 1, [ADF file support])
                 cs_have_adf=yes
               ], 
               [if test "x$with_adf" != "xcheck" ; then
                  AC_MSG_FAILURE([ADF support is requested, but test for ADF failed!])
                else
                  AC_MSG_WARN([no ADF support])
                fi
               ],
               )

  if test "x$cs_have_adf" = "xno"; then
    ADF_LIBS=""
  fi

  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(cs_have_adf)
AC_SUBST(ADF_LDFLAGS)
AC_SUBST(ADF_LIBS)
AC_SUBST(ADFRUNPATH)

])dnl


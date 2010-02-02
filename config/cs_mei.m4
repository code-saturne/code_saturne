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

# CS_AC_TEST_MEI
#---------------
# modifies or sets have_mei, MEI_CPPFLAGS, MEI_LDFLAGS, and MEI_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_MEI], [

have_mei=no

AC_ARG_WITH(mei,
            [AS_HELP_STRING([--with-mei=PATH],
                            [specify prefix directory for MEI])],
            [if test "x$withval" = "x"; then
               with_mei=yes
             fi],
            [with_mei=check])

AC_ARG_WITH(mei-include,
            [AS_HELP_STRING([--with-mei-include=PATH],
                            [specify directory for MEI include files])],
            [if test "x$with_mei" = "xcheck"; then
               with_mei=yes
             fi
             MEI_CPPFLAGS="-I$with_mei_include"],
            [if test "x$with_mei" != "xno" -a "x$with_mei" != "xyes" \
	          -a "x$with_mei" != "xcheck"; then
               MEI_CPPFLAGS="-I$with_mei/include"
             fi])

AC_ARG_WITH(mei-lib,
            [AS_HELP_STRING([--with-mei-lib=PATH],
                            [specify directory for MEI library])],
            [if test "x$with_mei" = "xcheck"; then
               with_mei=yes
             fi
             MEI_LDFLAGS="-L$with_mei_lib"],
            [if test "x$with_mei" != "xno" -a "x$with_mei" != "xyes" \
	          -a "x$with_mei" != "xcheck"; then
               MEI_LDFLAGS="-L$with_mei/lib"
             fi])


if test "x$with_mei" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$with_mei" != "xyes" -a "x$with_mei" != "xcheck" ; then
    mei_prefix=$with_mei
  fi

  MEI_LIBS="-lmei"

  CPPFLAGS="${CPPFLAGS} ${MEI_CPPFLAGS} ${BFT_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${MEI_LDFLAGS} ${BFT_LDFLAGS}"
  LIBS="${LIBS} ${MEI_LIBS} -lbft"

  AC_CHECK_HEADER([mei_evaluate.h])

  AC_CHECK_LIB(mei, mei_evaluate, 
               [ AC_DEFINE([HAVE_MEI], 1, [MEI support])
                 have_mei=yes
               ], 
               [AC_MSG_FAILURE([MEI support is requested, but test for MEI failed!])
               ],
               )

  if test "x$have_mei" != "xyes"; then
    MEI_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MEI, test x$have_mei = xyes)

AC_SUBST(MEI_CPPFLAGS)
AC_SUBST(MEI_LDFLAGS)
AC_SUBST(MEI_LIBS)

AC_SUBST(mei_prefix)

])dnl


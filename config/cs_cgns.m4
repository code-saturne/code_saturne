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

# CS_AC_TEST_CGNS
#----------------
# modifies or sets have_cgns, CGNS_CPPFLAGS, CGNS_LDFLAGS, and CGNS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CGNS], [

have_cgns=no

AC_ARG_WITH(cgns,
            [AS_HELP_STRING([--with-cgns=PATH],
                            [specify prefix directory for CGNS])],
            [if test "x$withval" = "x"; then
               with_cgns=yes
             fi],
            [with_cgns=check])

AC_ARG_WITH(cgns-include,
            [AS_HELP_STRING([--with-cgns-include=PATH],
                            [specify directory for CGNS include files])],
            [if test "x$with_cgns" = "xcheck"; then
               with_cgns=yes
             fi
             CGNS_CPPFLAGS="-I$with_cgns_include"],
            [if test "x$with_cgns" != "xno" -a "x$with_cgns" != "xyes" \
	          -a "x$with_cgns" != "xcheck"; then
               CGNS_CPPFLAGS="-I$with_cgns/include"
             fi])

AC_ARG_WITH(cgns-lib,
            [AS_HELP_STRING([--with-cgns-lib=PATH],
                            [specify directory for CGNS library])],
            [if test "x$with_cgns" = "xcheck"; then
               with_cgns=yes
             fi
             CGNS_LDFLAGS="-L$with_cgns_lib"],
            [if test "x$with_cgns" != "xno" -a "x$with_cgns" != "xyes" \
	          -a "x$with_cgns" != "xcheck"; then
               CGNS_LDFLAGS="-L$with_cgns/lib"
             fi])


if test "x$with_cgns" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CGNS_LIBS="-lcgns"
  CPPFLAGS="${CPPFLAGS} ${CGNS_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${CGNS_LDFLAGS} $HDF5_LDFLAGS"
  LIBS="${LIBS} ${CGNS_LIBS} $HDF5_LIBS"

  AC_CHECK_LIB(cgns, cg_coord_partial_write, 
               [ AC_DEFINE([HAVE_CGNS], 1, [CGNS file support])
                 have_cgns=yes
               ], 
               [if test "x$with_cgns" != "xcheck" ; then
                  AC_MSG_FAILURE([CGNS support is requested (requires CGNS >= 2.4), but test for CGNS failed!])
                else
                  AC_MSG_WARN([no CGNS file support (requires CGNS >= 2.4)])
                fi
               ],
               )

  if test "x$have_cgns" != "xyes"; then
    CGNS_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_CGNS, test x$have_cgns = xyes)

AC_SUBST(CGNS_CPPFLAGS)
AC_SUBST(CGNS_LDFLAGS)
AC_SUBST(CGNS_LIBS)

])dnl


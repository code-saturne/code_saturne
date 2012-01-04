dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2012 EDF S.A.
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

# CS_AC_TEST_CGNS
#----------------
# modifies or sets cs_have_cgns, CGNS_CPPFLAGS, CGNS_LDFLAGS, and CGNS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CGNS], [

cs_have_cgns=no

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
             CGNS_LDFLAGS="-L$with_cgns_lib"
             # Add the libdir to the runpath as CGNS is not libtoolized
             CGNSRUNPATH="-R$with_cgns_lib"],
            [if test "x$with_cgns" != "xno" -a "x$with_cgns" != "xyes" \
	          -a "x$with_cgns" != "xcheck"; then
               CGNS_LDFLAGS="-L$with_cgns/lib"
               # Add the libdir to the runpath as CGNS is not libtoolized
               CGNSRUNPATH="-R$with_cgns/lib"
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
                 cs_have_cgns=yes
               ], 
               [if test "x$with_cgns" != "xcheck" ; then
                  AC_MSG_FAILURE([CGNS support is requested (requires CGNS >= 2.4), but test for CGNS failed!])
                else
                  AC_MSG_WARN([no CGNS file support (requires CGNS >= 2.4)])
                fi
               ],
               )

  if test "x$cs_have_cgns" != "xyes"; then
    CGNS_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_CGNS, test x$cs_have_cgns = xyes)

AC_SUBST(cs_have_cgns)
AC_SUBST(CGNS_CPPFLAGS)
AC_SUBST(CGNS_LDFLAGS)
AC_SUBST(CGNS_LIBS)
AC_SUBST(CGNSRUNPATH)

])dnl


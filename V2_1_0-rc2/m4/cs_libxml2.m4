dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2009 EDF S.A., France
dnl
dnl   The Code_Saturne CFD tool is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne CFD tool is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne CFD tool; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_TEST_LIBXML2
#-------------------
# modifies or sets cs_have_libxml2, LIBXML2_CPPFLAGS, LIBXML2_LDFLAGS, and LIBXML2_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_LIBXML2], [

cs_have_libxml2=no

AC_ARG_WITH(libxml2,
            [AS_HELP_STRING([--with-libxml2=PATH],
                            [specify prefix directory for LIBXML2])],
            [if test "x$withval" = "x"; then
               with_libxml2=yes
             fi],
            [with_libxml2=check])

AC_ARG_WITH(libxml2-include,
            [AS_HELP_STRING([--with-libxml2-include=PATH],
                            [specify directory for LIBXML2 include files])],
            [if test "x$with_libxml2" = "xcheck"; then
               with_libxml2=yes
             fi
             LIBXML2_CPPFLAGS="-I$with_libxml2_include"],
            [if test "x$with_libxml2" != "xno" ; then
               if test "x$with_libxml2" != "xyes" \
	               -a "x$with_libxml2" != "xcheck"; then
                 LIBXML2_CPPFLAGS="-I$with_libxml2/include/libxml2"
               else
                 LIBXML2_CPPFLAGS="-I/usr/include/libxml2"
               fi
             fi])

AC_ARG_WITH(libxml2-lib,
            [AS_HELP_STRING([--with-libxml2-lib=PATH],
                            [specify directory for LIBXML2 library])],
            [if test "x$with_libxml2" = "xcheck"; then
               with_libxml2=yes
             fi
             LIBXML2_LDFLAGS="-L$with_libxml2_lib"],
            [if test "x$with_libxml2" != "xno" -a "x$with_libxml2" != "xyes" \
	          -a "x$with_libxml2" != "xcheck"; then
               LIBXML2_LDFLAGS="-L$with_libxml2/lib"
             fi])


if test "x$with_libxml2" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  LIBXML2_LIBS="-lxml2"

  CPPFLAGS="${CPPFLAGS} ${LIBXML2_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${LIBXML2_LDFLAGS}"
  LIBS="${LIBS} ${LIBXML2_LIBS}"

  AC_CHECK_HEADER([libxml/parser.h])

  AC_CHECK_LIB(xml2, xmlInitParser, 
               [ AC_DEFINE([HAVE_LIBXML2], 1, [LIBXML2 support])
                 cs_have_libxml2=yes
               ], 
               [if test "x$with_libxml2" != "xcheck" ; then
                  AC_MSG_FAILURE([LIBXML2 support is requested, but test for LIBXML2 failed!])
                else
                  AC_MSG_WARN([no LIBXML2 support])
                fi
               ],
              )

  if test "x$cs_have_libxml2" != "xyes"; then
    LIBXML2_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_LIBXML2, test x$cs_have_libxml2 = xyes)

AC_SUBST(cs_have_libxml2)
AC_SUBST(LIBXML2_CPPFLAGS)
AC_SUBST(LIBXML2_LDFLAGS)
AC_SUBST(LIBXML2_LIBS)

])dnl


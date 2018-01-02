dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2018 EDF S.A.
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

# CS_AC_TEST_LIBXML2
#-------------------
# modifies or sets cs_have_libxml2, LIBXML2_CPPFLAGS, LIBXML2_LDFLAGS, and LIBXML2_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_LIBXML2], [

cs_have_libxml2=no
cs_have_libxml2_header=no
libxml2_prefix=""

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

  LDFLAGS="${LIBXML2_LDFLAGS} ${LDFLAGS}"
  LIBS="${LIBXML2_LIBS} ${LIBS} -lm"

  CPPFLAGS="${saved_CPPFLAGS} ${LIBXML2_CPPFLAGS}"
  AC_CHECK_HEADERS([libxml/parser.h],
                   [cs_have_libxml2_header=yes],
                   [],
                   [])

  # If header not found, try other standard configurations

  if test "x$cs_have_libxml2_header" = "xno" ; then
    unset ac_cv_header_libxml_parser_h
    LIBXML2_CPPFLAGS="-I/mingw/include/libxml2"
    CPPFLAGS="${saved_CPPFLAGS} ${LIBXML2_CPPFLAGS}"
    AC_CHECK_HEADERS([libxml/parser.h],
                     [cs_have_libxml2_header=yes],
                     [],
                     [])
  fi

  if test "x$cs_have_libxml2_header" = "xyes" ; then

    AC_CHECK_LIB(xml2, xmlInitParser,
                 [ AC_DEFINE([HAVE_LIBXML2], 1, [LIBXML2 support])
                   cs_have_libxml2=yes
                 ],
                 [AC_MSG_FAILURE([LIBXML2 support is requested, but test for LIBXML2 failed!])
                 ],
                )

  else
    AC_MSG_FAILURE([LIBXML2 support is requested, but LIBXML2 headers not found; ])
  fi

  if test "x$cs_have_libxml2" != "xyes"; then
    LIBXML2_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  case $host_os in
    mingw32)
      libxml2_prefix=`cygpath --path --windows "$with_libxml2"`;;
    *)
      ;;
  esac
fi

AM_CONDITIONAL(HAVE_LIBXML2, test x$cs_have_libxml2 = xyes)

AC_SUBST(cs_have_libxml2)
AC_SUBST(libxml2_prefix, [${libxml2_prefix}])
AC_SUBST(LIBXML2_CPPFLAGS)
AC_SUBST(LIBXML2_LDFLAGS)
AC_SUBST(LIBXML2_LIBS)

])dnl


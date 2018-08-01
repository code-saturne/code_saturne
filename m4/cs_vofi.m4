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

# CS_AC_TEST_VOFI
#----------------
# modifies or sets cs_have_vofi, VOFI_CPPFLAGS, VOFI_LDFLAGS, and VOFI_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_VOFI], [

cs_have_vofi=no
cs_have_vofi_headers=no
vofi_prefix=""

AC_ARG_WITH(vofi,
            [AS_HELP_STRING([--with-vofi=PATH],
                            [specify prefix directory for VOFI])],
            [if test "x$withval" = "x"; then
               with_vofi=yes
             fi],
            [with_vofi=check])

AC_ARG_WITH(vofi-include,
            [AS_HELP_STRING([--with-vofi-include=PATH],
                            [specify directory for VOFI include files])],
            [if test "x$with_vofi" = "xcheck"; then
               with_vofi=yes
             fi
             VOFI_CPPFLAGS="-I$with_vofi_include"],
            [if test "x$with_vofi" != "xno" -a "x$with_vofi" != "xyes" \
	          -a "x$with_vofi" != "xcheck"; then
               VOFI_CPPFLAGS="-I$with_vofi/include"
             fi])

AC_ARG_WITH(vofi-lib,
            [AS_HELP_STRING([--with-vofi-lib=PATH],
                            [specify directory for VOFI library])],
            [if test "x$with_vofi" = "xcheck"; then
               with_vofi=yes
             fi
             VOFI_LDFLAGS="-L$with_vofi_lib"
             # Add the libdir to the runpath as VOFI is not libtoolized
             VOFIRUNPATH="-R$with_vofi_lib"],
            [if test "x$with_vofi" != "xno" -a "x$with_vofi" != "xyes" \
	          -a "x$with_vofi" != "xcheck"; then
               VOFI_LDFLAGS="-L$with_vofi/lib"
               # Add the libdir to the runpath as VOFI is not libtoolized
               VOFIRUNPATH="-R$with_vofi/lib"
             fi])


if test "x$with_vofi" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  VOFI_LIBS="-lvofi"
  CPPFLAGS="${CPPFLAGS} ${VOFI_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${VOFI_LDFLAGS}"
  LIBS="${VOFI_LIBS} ${LIBS}"

  # Check for a VOFI library
  #------------------------------

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <vofi.h>
]])],
                    [AC_MSG_RESULT([Vofi headers found])
                     cs_have_vofi_headers=yes
                    ],
                    [AC_MSG_RESULT([Vofi headers not found])
                    ])

  if test "x$cs_have_vofi_headers" = "xyes"; then

    AC_CHECK_LIB(vofi, vofi_Get_fh,
                 [ AC_DEFINE([HAVE_VOFI], 1, [VOFI support])
                   cs_have_vofi=yes
                 ],
                 [])

  fi

  if test "x$cs_have_vofi" != "xyes"; then
    VOFI_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  # Report VOFI support
  #-------------------

  if test "x$cs_have_vofi" = "xno" ; then
    if test "x$with_vofi" != "xcheck" ; then
      AC_MSG_FAILURE([VOFI support is requested, but test for VOFI failed!])
    else
      AC_MSG_WARN([no VOFI support])
    fi
  fi

  case $host_os in
    mingw32)
      vofi_prefix=`cygpath --path --windows "$with_vofi"`;;
    *)
      ;;
  esac

fi

unset cs_have_vofi_headers

AM_CONDITIONAL(HAVE_VOFI, test x$cs_have_vofi = xyes)

AC_SUBST(cs_have_vofi)
AC_SUBST(vofi_prefix, [${vofi_prefix}])
AC_SUBST(VOFI_CPPFLAGS)
AC_SUBST(VOFI_LDFLAGS)
AC_SUBST(VOFI_LIBS)
AC_SUBST(VOFIRUNPATH)

])dnl


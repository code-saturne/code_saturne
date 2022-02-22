dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2022 EDF S.A.
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

# CS_AC_TEST_AMGX
#----------------
# modifies or sets cs_have_amgx, AMGX_CPPFLAGS, AMGX_LDFLAGS, and AMGX_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_AMGX], [

cs_have_amgx=no
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(amgx,
            [AS_HELP_STRING([--with-amgx=PATH],
                            [specify prefix directory for AMGX])],
            [if test "x$withval" = "x"; then
               with_amgx=no
             fi],
            [with_amgx=no])

AC_ARG_WITH(amgx-include,
            [AS_HELP_STRING([--with-amgx-include=PATH],
                            [specify directory for AMGX include files])],
            [if test "x$with_amgx" = "xcheck"; then
               with_amgx=yes
             fi
             AMGX_CPPFLAGS="-I$with_amgx_include"],
            [if test "x$with_amgx" != "xno" -a "x$with_amgx" != "xyes" \
	          -a "x$with_amgx" != "xcheck"; then
               AMGX_CPPFLAGS="-I$with_amgx/include"
             fi])

AC_ARG_WITH(amgx-lib,
            [AS_HELP_STRING([--with-amgx-lib=PATH],
                            [specify directory for AMGX library])],
            [if test "x$with_amgx" = "xcheck"; then
               with_amgx=yes
             fi
             AMGX_LDFLAGS="-L$with_amgx_lib"
             cs_amgx_libpath="$with_amgx_lib"
             # Add the libdir to the runpath as AMGX might not be libtoolized
             AMGXRUNPATH="-R$with_amgx_lib"],
            [if test "x$with_amgx" != "xno" -a "x$with_amgx" != "xyes" \
	          -a "x$with_amgx" != "xcheck"; then
               AMGX_LDFLAGS="-L$with_amgx/lib"
               # Add the libdir to the runpath as AMGX might not be libtoolized
               AMGXRUNPATH="-R$with_amgx/lib"
               cs_amgx_libpath="$with_amgx/lib"
             fi])

if test "x$with_amgx" != "xno" ; then

  # Now run tests

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${AMGX_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${AMGX_LDFLAGS} ${MPI_LDFLAGS}"

  if test "x$cs_amgx_libpath" != x ; then
    if test ! -d "$cs_amgx_libpath" ; then
      AC_MSG_FAILURE([directory specified by --with-amgx-lib=$cs_amgx_libpath does not exist!])
    fi
  fi

  unset cs_amgx_libname

  AC_MSG_CHECKING([for AMGX (dynamic library)])

  AMGX_LIBS="-lamgxsh"
  LIBS="${AMGX_LIBS} ${saved_LIBS}"

  # Now check library

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <amgx_c.h>]],
                 [[ AMGX_initialize(); ]])],
                 [AC_DEFINE([HAVE_AMGX], 1, [AMGX library support])
                  cs_have_amgx=yes],
                 [cs_have_amgx=no])

  AC_MSG_RESULT($cs_have_amgx)

  if test "x$cs_have_amgx" = "xno"; then

    AC_MSG_CHECKING([for AMGX (static library; may need flags for cuSPARSE)])

    AMGX_LIBS="-lamgx"
    LIBS="${AMGX_LIBS} ${saved_LIBS} ${MPI_LIBS} -lm"

    # Now check library

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <amgx_c.h>]],
                   [[ AMGX_initialize(); ]])],
                   [AC_DEFINE([HAVE_AMGX], 1, [AMGX library support])
                    cs_have_amgx=yes],
                   [cs_have_amgx=no])

    AC_MSG_RESULT($cs_have_amgx)

  fi

  if test "x$cs_have_amgx" = "xno"; then
    AMGX_CPPFLAGS=""
    AMGX_LDFLAGS=""
    AMGX_LIBS=""
    AMGX_RUNPATH=""
    if test "x$with_amgx" != "xcheck" ; then
      AC_MSG_FAILURE([AMGX support is requested, but test for AMGX failed!])
    else
      AC_MSG_WARN([no AMGX file support])
    fi
  fi

  unset cs_amgx_libnames
  unset cs_amgx_libpath

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_AMGX, test x$cs_have_amgx = xyes)

AC_SUBST(cs_have_amgx)
AC_SUBST(AMGX_CPPFLAGS)
AC_SUBST(AMGX_LDFLAGS)
AC_SUBST(AMGX_LIBS)
AC_SUBST(AMGXRUNPATH)

])dnl

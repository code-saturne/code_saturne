dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2025 EDF S.A.
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

# CS_AC_TEST_CUDSS
#----------------
# modifies or sets cs_have_cudss, CUDSS_CPPFLAGS, CUDSS_LDFLAGS, and CUDSS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CUDSS], [

cs_have_cudss=no
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(cudss,
            [AS_HELP_STRING([--with-cudss=PATH],
                            [specify prefix directory for cuDSS])],
            [if test "x$withval" = "x"; then
               with_cudss=no
             fi],
            [with_cudss=no])

AC_ARG_WITH(cudss-include,
            [AS_HELP_STRING([--with-cudss-include=PATH],
                            [specify directory for cuDSS include files])],
            [if test "x$with_cudss" = "xcheck"; then
               with_cudss=yes
             fi
             CUDSS_CPPFLAGS="-I$with_cudss_include"],
            [if test "x$with_cudss" != "xno" -a "x$with_cudss" != "xyes" \
	          -a "x$with_cudss" != "xcheck"; then
               CUDSS_CPPFLAGS="-I$with_cudss/include"
             fi])

AC_ARG_WITH(cudss-lib,
            [AS_HELP_STRING([--with-cudss-lib=PATH],
                            [specify directory for cuDSS library])],
            [if test "x$with_cudss" = "xcheck"; then
               with_cudss=yes
             fi
             CUDSS_LDFLAGS="-L$with_cudss_lib"
             cs_cudss_libpath="$with_cudss_lib"],
            [if test "x$with_cudss" != "xno" -a "x$with_cudss" != "xyes" \
	          -a "x$with_cudss" != "xcheck"; then
               CUDSS_LDFLAGS="-L$with_cudss/lib"
               cs_cudss_libpath="$with_cudss/lib"
             fi])

if test "x$with_cudss" != "xno" ; then

  # Now run tests

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${CUDSS_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${CUDSS_LDFLAGS} ${MPI_LDFLAGS}"

  if test "x$cs_cudss_libpath" != x ; then
    if test ! -d "$cs_cudss_libpath" ; then
      AC_MSG_FAILURE([directory specified by --with-cuDSS-lib=$cs_cudss_libpath does not exist!])
    fi
  fi

  unset cs_cudss_libname

  AC_MSG_CHECKING([for cuDSS])

  CUDSS_LIBS="-lcudss"
  LIBS="${CUDSS_LIBS} ${saved_LIBS}"

  # Now check library

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cudss.h>]],
                 [[cudssHandle_t handle;
                 cudssCreate(&handle);]])],
                 [AC_DEFINE([HAVE_CUDSS], 1, [CUDSS library support])
                  cs_have_cudss=yes],
                 [cs_have_cudss=no])

  AC_MSG_RESULT($cs_have_cudss)

  if test "x$cs_have_cudss" = "xno"; then
    CUDSS_CPPFLAGS=""
    CUDSS_LDFLAGS=""
    CUDSS_LIBS=""
    if test "x$with_cudss" != "xcheck" ; then
      AC_MSG_FAILURE([cuDSS support is requested, but test for cuDSS failed!])
    else
      AC_MSG_WARN([no cuDSS file support])
    fi
  fi

  unset cs_cudss_libnames
  unset cs_cudss_libpath

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_CUDSS, test x$cs_have_cudss = xyes)

AC_SUBST(cs_have_cudss)
AC_SUBST(CUDSS_CPPFLAGS)
AC_SUBST(CUDSS_LDFLAGS)
AC_SUBST(CUDSS_LIBS)

])dnl

dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2026 EDF S.A.
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

# CS_AC_TEST_GINKGO
#------------------
# modifies or sets cs_have_ginkgo, GINKGO_CPPFLAGS, GINKGO_LDFLAGS, and GINKGO_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_GINKGO], [

cs_have_ginkgo=no
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(ginkgo,
            [AS_HELP_STRING([--with-ginkgo=PATH],
                            [specify prefix directory for GINKGO])],
            [if test "x$withval" = "x"; then
               with_ginkgo=no
             fi],
            [with_ginkgo=no])

AC_ARG_WITH(ginkgo-include,
            [AS_HELP_STRING([--with-ginkgo-include=PATH],
                            [specify directory for GINKGO include files])],
            [if test "x$with_ginkgo" = "xcheck"; then
               with_ginkgo=yes
             fi
             GINKGO_CPPFLAGS="-I$with_ginkgo_include"],
            [if test "x$with_ginkgo" != "xno" -a "x$with_ginkgo" != "xyes" \
                  -a "x$with_ginkgo" != "xcheck"; then
               GINKGO_CPPFLAGS="-I$with_ginkgo/include"
             fi])

AC_ARG_WITH(ginkgo-lib,
            [AS_HELP_STRING([--with-ginkgo-lib=PATH],
                            [specify directory for GINKGO library])],
            [if test "x$with_ginkgo" = "xcheck"; then
               with_ginkgo=yes
             fi
             GINKGO_LDFLAGS="-L$with_ginkgo_lib"
             cs_ginkgo_libpath="$with_ginkgo_lib"],
            [if test "x$with_ginkgo" != "xno" -a "x$with_ginkgo" != "xyes" \
                  -a "x$with_ginkgo" != "xcheck"; then
               GINKGO_LDFLAGS="-L$with_ginkgo/lib"
               cs_ginkgo_libpath="$with_ginkgo/lib"
             fi])

if test "x$with_ginkgo" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_CXXFLAGS="$CXXFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # Now check library

  AC_MSG_CHECKING([for Ginkgo])

  GINKGO_LIBS="-lginkgo -lginkgo_device -lginkgo_omp -lginkgo_cuda -ldl -lginkgo_reference -lginkgo_hip -lginkgo_dpcpp"
  LIBS="${GINKGO_LIBS} ${saved_LIBS}"
  cs_ginkgo_cxxflags="`echo ${CXXFLAGS} | sed 's/-Werror=shadow/-Wshadow/'`"
  CPPFLAGS="${saved_CPPFLAGS} ${GINKGO_CPPFLAGS} ${MPI_CPPFLAGS}"
  CXXFLAGS="${cs_ginkgo_cxxflags}"
  LDFLAGS="${saved_LDFLAGS} ${GINKGO_LDFLAGS} ${MPI_LDFLAGS}"
  LIBS="${saved_LIBS} ${GINKGO_LIBS} ${MPI_LIBS} -lm"

  AC_LANG_PUSH([C++])

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <ginkgo/ginkgo.hpp>]],
                 [[const auto exec = gko::ReferenceExecutor::create();]])],
                 [AC_DEFINE([HAVE_GINKGO], 1, [GINKGO library support])
                  cs_have_ginkgo=yes],
                 [cs_have_ginkgo=no])

  AC_LANG_POP([C++])

  unset cs_ginkgo_cxxflags

  AC_MSG_RESULT($cs_have_ginkgo)

  if test "x$cs_have_ginkgo" = "xno"; then
    GINKGO_CPPFLAGS=""
    GINKGO_LDFLAGS=""
    GINKGO_LIBS=""
    if test "x$with_ginkgo" != "xcheck" ; then
      AC_MSG_FAILURE([GINKGO support is requested, but test for GINKGO failed!])
    else
      AC_MSG_WARN([no GINKGO file support])
    fi
  fi

  unset cs_ginkgo_libpath

  CPPFLAGS="$saved_CPPFLAGS"
  CXXFLAGS="$saved_CXXFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_GINKGO, test x$cs_have_ginkgo = xyes)

AC_SUBST(cs_have_ginkgo)
AC_SUBST(GINKGO_CPPFLAGS)
AC_SUBST(GINKGO_CXXFLAGS)
AC_SUBST(GINKGO_LDFLAGS)
AC_SUBST(GINKGO_LIBS)

])dnl

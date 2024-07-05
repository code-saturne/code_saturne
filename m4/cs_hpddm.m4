dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2024 EDF S.A.
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

# CS_AC_TEST_HPDDM
#----------------
# modifies or sets cs_have_hpddm, HPDDM_CPPFLAGS, HPDDM_LDFLAGS, and HPDDM_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_HPDDM], [

cs_have_hpddm=no
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(hpddm,
            [AS_HELP_STRING([--with-hpddm=PATH],
                            [specify prefix directory for HPDDM])],
            [if test "x$withval" = "x"; then
               with_hpddm=no
             fi],
            [with_hpddm=no])

AC_ARG_WITH(hpddm-include,
            [AS_HELP_STRING([--with-hpddm-include=PATH],
                            [specify directory for HPDDM include files])],
            [if test "x$with_hpddm" = "xcheck"; then
               with_hpddm=yes
             fi
             HPDDM_CPPFLAGS="-I$with_hpddm_include"],
            [if test "x$with_hpddm" != "xno" -a "x$with_hpddm" != "xyes" \
                  -a "x$with_hpddm" != "xcheck"; then
               HPDDM_CPPFLAGS="-I$with_hpddm/include"
             fi])

AC_ARG_WITH(hpddm-lib,
            [AS_HELP_STRING([--with-hpddm-lib=PATH],
                            [specify directory for HPDDM library])],
            [if test "x$with_hpddm" = "xcheck"; then
               with_hpddm=yes
             fi
             HPDDM_LDFLAGS="-L$with_hpddm_lib"
             cs_hpddm_libpath="$with_hpddm_lib"],
            [if test "x$with_hpddm" != "xno" -a "x$with_hpddm" != "xyes" \
                  -a "x$with_hpddm" != "xcheck"; then
               HPDDM_LDFLAGS="-L$with_hpddm/lib"
               cs_hpddm_libpath="$with_hpddm/lib"
             fi])

if test "x$with_hpddm" != "xno" ; then

  # If PETSc is used and HPDDM unspecified, check
  # if the associated flags include HPDDM

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$cs_have_petsc" = "xyes" -a "x$with_hpddm" = "xcheck"; then

    CPPFLAGS="${saved_CPPFLAGS} ${HPDDM_CPPFLAGS} ${PETSC_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${HPDDM_LDFLAGS} ${PETSC_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS"${saved_LIBS} ${HPDDM_LIBS} ${PETSC_LIBS} ${MPI_LIBS}"

    AC_MSG_CHECKING([for HPDDM (using PETSC flags)])

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <HPDDM.hpp>]],
                   [[ HPDDM_STR(COARSEOPERATOR); ]])],
                   [AC_DEFINE([HAVE_HPDDM], 1, [HPDDM library support])
                    cs_have_hpddm=yes],
                   [cs_have_hpddm=no])

     AC_MSG_RESULT($cs_have_hpddm)

  fi

  if test "x$cs_hpddm_libpath" != x ; then
    if test ! -d "$cs_hpddm_libpath" ; then
      AC_MSG_FAILURE([directory specified by --with-hpddm-lib=$cs_hpddm_libpath does not exist!])
    fi
  fi

  HPDDM_LIBS="-lhpddm_petsc"
  LIBS="${HPDDM_LIBS} ${saved_LIBS}"

  # Now check library

  if test "x$cs_have_hpddm" = "xno"; then

    AC_MSG_CHECKING([for HPDDM])

    HPDDM_LIBS="-lhpddm_petsc"

    CPPFLAGS="${saved_CPPFLAGS} ${HPDDM_CPPFLAGS} ${MPI_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${HPDDM_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS="${saved_LIBS} ${HPDDM_LIBS} ${MPI_LIBS} -lm"

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <HPDDM.hpp>]],
                   [[ HPDDM_STR(COARSEOPERATOR); ]])],
                   [AC_DEFINE([HAVE_HPDDM], 1, [HPDDM library support])
                    cs_have_hpddm=yes],
                   [cs_have_hpddm=no])

    AC_MSG_RESULT($cs_have_hpddm)

  fi

  if test "x$cs_have_hpddm" = "xno"; then
    HPDDM_CPPFLAGS=""
    HPDDM_LDFLAGS=""
    HPDDM_LIBS=""
    if test "x$with_hpddm" != "xcheck" ; then
      AC_MSG_FAILURE([HPDDM support is requested, but test for HPDDM failed!])
    else
      AC_MSG_WARN([no HPDDM file support])
    fi
  fi

  unset cs_hpddm_libpath

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_HPDDM, test x$cs_have_hpddm = xyes)

AC_SUBST(cs_have_hpddm)
AC_SUBST(HPDDM_CPPFLAGS)
AC_SUBST(HPDDM_LDFLAGS)
AC_SUBST(HPDDM_LIBS)

])dnl

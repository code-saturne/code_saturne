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

# CS_AC_TEST_SLEPC
#----------------
# modifies or sets cs_have_slepc, SLEPC_CPPFLAGS, SLEPC_LDFLAGS, and SLEPC_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SLEPC], [

cs_have_slepc=no
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(slepc,
            [AS_HELP_STRING([--with-slepc=PATH],
                            [specify prefix directory for SLEPC])],
            [if test "x$withval" = "x"; then
               with_slepc=no
             fi],
            [with_slepc=no])

AC_ARG_WITH(slepc-include,
            [AS_HELP_STRING([--with-slepc-include=PATH],
                            [specify directory for SLEPC include files])],
            [if test "x$with_slepc" = "xcheck"; then
               with_slepc=yes
             fi
             SLEPC_CPPFLAGS="-I$with_slepc_include"],
            [if test "x$with_slepc" != "xno" -a "x$with_slepc" != "xyes" \
                  -a "x$with_slepc" != "xcheck"; then
               SLEPC_CPPFLAGS="-I$with_slepc/include"
             fi])

AC_ARG_WITH(slepc-lib,
            [AS_HELP_STRING([--with-slepc-lib=PATH],
                            [specify directory for SLEPC library])],
            [if test "x$with_slepc" = "xcheck"; then
               with_slepc=yes
             fi
             SLEPC_LDFLAGS="-L$with_slepc_lib"
             cs_slepc_libpath="$with_slepc_lib"],
            [if test "x$with_slepc" != "xno" -a "x$with_slepc" != "xyes" \
                  -a "x$with_slepc" != "xcheck"; then
               SLEPC_LDFLAGS="-L$with_slepc/lib"
               cs_slepc_libpath="$with_slepc/lib"
             fi])

if test "x$with_slepc" != "xno" ; then

  # If PETSc is used and SLEPC unspecified, check
  # if the associated flags include SLEPC

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$cs_have_petsc" = "xyes" -a "x$with_slepc" = "xcheck"; then

    CPPFLAGS="${saved_CPPFLAGS} ${SLEPC_CPPFLAGS} ${PETSC_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${SLEPC_LDFLAGS} ${PETSC_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS"${saved_LIBS} ${SLEPC_LIBS} ${PETSC_LIBS} ${MPI_LIBS}"

    AC_MSG_CHECKING([for SLEPC (using PETSC flags)])

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <slepc.h>]],
                   [[ ]])],
                   [AC_DEFINE([HAVE_SLEPC], 1, [SLEPC library support])
                    cs_have_slepc=yes],
                   [cs_have_slepc=no])

     AC_MSG_RESULT($cs_have_slepc)

  fi

  if test "x$cs_slepc_libpath" != x ; then
    if test ! -d "$cs_slepc_libpath" ; then
      AC_MSG_FAILURE([directory specified by --with-slepc-lib=$cs_slepc_libpath does not exist!])
    fi
  fi

  SLEPC_LIBS="-lslepc"
  LIBS="${SLEPC_LIBS} ${saved_LIBS}"

  # Now check library

  if test "x$cs_have_slepc" = "xno"; then

    AC_MSG_CHECKING([for SLEPC])

    SLEPC_LIBS="-lslepc"

    CPPFLAGS="${saved_CPPFLAGS} ${SLEPC_CPPFLAGS} ${MPI_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${SLEPC_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS="${saved_LIBS} ${SLEPC_LIBS} ${MPI_LIBS} -lm"

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <slepc.h>]],
                   [[ ]])],
                   [AC_DEFINE([HAVE_SLEPC], 1, [SLEPC library support])
                    cs_have_slepc=yes],
                   [cs_have_slepc=no])

    AC_MSG_RESULT($cs_have_slepc)

  fi

  if test "x$cs_have_slepc" = "xno"; then
    SLEPC_CPPFLAGS=""
    SLEPC_LDFLAGS=""
    SLEPC_LIBS=""
    if test "x$with_slepc" != "xcheck" ; then
      AC_MSG_FAILURE([SLEPC support is requested, but test for SLEPC failed!])
    else
      AC_MSG_WARN([no SLEPC file support])
    fi
  fi

  unset cs_slepc_libpath

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_SLEPC, test x$cs_have_slepc = xyes)

AC_SUBST(cs_have_slepc)
AC_SUBST(SLEPC_CPPFLAGS)
AC_SUBST(SLEPC_LDFLAGS)
AC_SUBST(SLEPC_LIBS)

])dnl

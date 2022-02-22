dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2020 EDF S.A.
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

# CS_AC_TEST_HYPRE
#----------------
# modifies or sets cs_have_hypre, HYPRE_CPPFLAGS, HYPRE_LDFLAGS, and HYPRE_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_HYPRE], [

cs_have_hypre=no
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(hypre,
            [AS_HELP_STRING([--with-hypre=PATH],
                            [specify prefix directory for HYPRE])],
            [if test "x$withval" = "x"; then
               with_hypre=no
             fi],
            [with_hypre=no])

AC_ARG_WITH(hypre-include,
            [AS_HELP_STRING([--with-hypre-include=PATH],
                            [specify directory for HYPRE include files])],
            [if test "x$with_hypre" = "xcheck"; then
               with_hypre=yes
             fi
             HYPRE_CPPFLAGS="-I$with_hypre_include"],
            [if test "x$with_hypre" != "xno" -a "x$with_hypre" != "xyes" \
	          -a "x$with_hypre" != "xcheck"; then
               HYPRE_CPPFLAGS="-I$with_hypre/include"
             fi])

AC_ARG_WITH(hypre-lib,
            [AS_HELP_STRING([--with-hypre-lib=PATH],
                            [specify directory for HYPRE library])],
            [if test "x$with_hypre" = "xcheck"; then
               with_hypre=yes
             fi
             HYPRE_LDFLAGS="-L$with_hypre_lib"
             cs_hypre_libpath="$with_hypre_lib"
             # Add the libdir to the runpath as HYPRE might not be libtoolized
             HYPRERUNPATH="-R$with_hypre_lib"],
            [if test "x$with_hypre" != "xno" -a "x$with_hypre" != "xyes" \
	          -a "x$with_hypre" != "xcheck"; then
               HYPRE_LDFLAGS="-L$with_hypre/lib"
               # Add the libdir to the runpath as HYPRE might not be libtoolized
               HYPRERUNPATH="-R$with_hypre/lib"
               cs_hypre_libpath="$with_hypre/lib"
             fi])

if test "x$with_hypre" != "xno" ; then

  # If PETSc is used and HYPRE unspecified, check
  # if the associated flags include HYPRE

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$cs_have_petsc" = "xyes" -a "x$with_hypre" = "xcheck"; then

    CPPFLAGS="${saved_CPPFLAGS} ${HYPRE_CPPFLAGS} ${PETSC_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${HYPRE_LDFLAGS} ${PETSC_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS"${saved_LIBS} ${HYPRE_LIBS} ${PETSC_LIBS} ${MPI_LIBS}"

    AC_MSG_CHECKING([for HYPRE (using PETSC flags)])

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <HYPRE_parcsr_ls.h>]],
                   [[ HYPRE_Solver solver; HYPRE_BoomerAMGCreate(&solver); ]])],
                   [AC_DEFINE([HAVE_HYPRE], 1, [HYPRE library support])
                    cs_have_hypre=yes],
                   [cs_have_hypre=no])

     AC_MSG_RESULT($cs_have_hypre)

  fi

  if test "x$cs_hypre_libpath" != x ; then
    if test ! -d "$cs_hypre_libpath" ; then
      AC_MSG_FAILURE([directory specified by --with-hypre-lib=$cs_hypre_libpath does not exist!])
    fi
  fi

  HYPRE_LIBS="-lHYPRE"
  LIBS="${HYPRE_LIBS} ${saved_LIBS}"

  # Now check library

  if test "x$cs_have_hypre" = "xno"; then

    AC_MSG_CHECKING([for HYPRE])

    HYPRE_LIBS="-lHYPRE"

    CPPFLAGS="${saved_CPPFLAGS} ${HYPRE_CPPFLAGS} ${MPI_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${HYPRE_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS="${saved_LIBS} ${HYPRE_LIBS} ${MPI_LIBS} -lm"

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <HYPRE_parcsr_ls.h>]],
                   [[ HYPRE_Solver solver; HYPRE_BoomerAMGCreate(&solver); ]])],
                   [AC_DEFINE([HAVE_HYPRE], 1, [HYPRE library support])
                    cs_have_hypre=yes],
                   [cs_have_hypre=no])

    AC_MSG_RESULT($cs_have_hypre)

  fi

  if test "x$cs_have_hypre" = "xno"; then
    HYPRE_CPPFLAGS=""
    HYPRE_LDFLAGS=""
    HYPRE_LIBS=""
    HYPRE_RUNPATH=""
    if test "x$with_hypre" != "xcheck" ; then
      AC_MSG_FAILURE([HYPRE support is requested, but test for HYPRE failed!])
    else
      AC_MSG_WARN([no HYPRE file support])
    fi
  fi

  unset cs_hypre_libpath

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_HYPRE, test x$cs_have_hypre = xyes)

AC_SUBST(cs_have_hypre)
AC_SUBST(HYPRE_CPPFLAGS)
AC_SUBST(HYPRE_LDFLAGS)
AC_SUBST(HYPRE_LIBS)
AC_SUBST(HYPRERUNPATH)

])dnl

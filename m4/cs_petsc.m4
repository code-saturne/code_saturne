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

# CS_AC_TEST_PETSC
#----------------
# modifies or sets cs_have_petsc, PETSC_CPPFLAGS, PETSC_LDFLAGS, and PETSC_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_PETSC], [

cs_have_petsc=no
cs_have_petsc_header=no
petsc_prefix=""
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(petsc,
            [AS_HELP_STRING([--with-petsc=PATH],
                            [specify prefix directory for PETSc])],
            [if test "x$withval" = "x"; then
               with_petsc=no
             fi],
            [with_petsc=no])

AC_ARG_WITH(petsc-lib,
            [AS_HELP_STRING([--with-petsc-lib=PATH],
                            [specify directory for PETSc library])],
            [if test "x$with_petsc" != "xno"; then
               with_petsc=yes  saved_LIBS="$LIBS"
             fi
             PETSC_DIR="$with_petsc_lib/petsc"],
            [if test "x$with_petsc" != "xno" -a "x$with_petsc" != "xyes" \
	          -a "x$with_petsc" != "xcheck"; then
               PETSC_DIR="$with_petsc/lib/petsc"
             fi])

if test "x$with_petsc" != "xno" ; then
  cs_petsc_make=`which gmake`
  if test $? = 1 ; then
    cs_petsc_make=make
  fi
  if test -f ${PETSC_DIR}/conf/variables ; then
    PETSC_CPPFLAGS=$(${cs_petsc_make} -s -f "$cs_abs_srcdir/build-aux/petsc-variables.makefile" PETSC_DIR="${PETSC_DIR}" getincludedirs)
    PETSC_LDFLAGS=""
    PETSC_LIBS=$(${cs_petsc_make} -s -f "$cs_abs_srcdir/build-aux/petsc-variables.makefile"  PETSC_DIR="${PETSC_DIR}" getlinklibs)
  elif test -f ${PETSC_DIR}/conf/petscvariables ; then
    PETSC_CPPFLAGS=$(${cs_petsc_make} -s -f "$cs_abs_srcdir/build-aux/petsc-petscvariables.makefile" PETSC_DIR="${PETSC_DIR}" getincludedirs)
    PETSC_LDFLAGS=""
    PETSC_LIBS=$(${cs_petsc_make} -s -f "$cs_abs_srcdir/build-aux/petsc-petscvariables.makefile" PETSC_DIR="${PETSC_DIR}" getlinklibs)
  else
      AC_MSG_FAILURE([${PETSC_DIR}/conf/variables or ${PETSC_DIR}/conf/petscvariables not found.
Check --with-petsc or --with-petsc-lib option or PETSc directory structure
({--with-petsc}/lib/conf/variables or {--with-petsc_lib}/conf/variables or
 {--with-petsc}/lib/conf/petscvariables or {--with-petsc_lib}/conf/petscvariables
should be present).])
  fi

  unset petsc_make

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${PETSC_CPPFLAGS} ${MPI_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${MPI_LDFLAGS}"
  LIBS="${LIBS} ${PETSC_LIBS} ${MPI_LIBS}"

  AC_MSG_CHECKING([for PETSc library])

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <petscsys.h>]],
[[PetscInitializeNoArguments();]])
                   ],
                   [ AC_DEFINE([HAVE_PETSC], 1, [PETSc support])
                     cs_have_petsc=yes
                   ],
                   [ AC_MSG_WARN([no PETSc support])
                     cs_have_petsc=no
                   ],
                  )

  AC_MSG_RESULT($cs_have_petsc)

  if test "x$cs_have_petsc" = "xno"; then
    PETSC_CPPFLAGS=""
    PETSC_LDFLAGS=""
    PETSC_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_PETSC, test x$cs_have_petsc = xyes)

AC_SUBST(cs_have_petsc)
AC_SUBST(PETSC_CPPFLAGS)
AC_SUBST(PETSC_LDFLAGS)
AC_SUBST(PETSC_LIBS)
])dnl


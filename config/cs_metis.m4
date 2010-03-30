dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2010 EDF S.A., France
dnl
dnl   The Code_Saturne Kernel is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne Kernel is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne Preprocessor; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_TEST_METIS
#-----------------
# modifies or sets have_metis, METIS_CPPFLAGS, METIS_LDFLAGS, and METIS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_METIS], [

have_parmetis=no
have_metis=no

AC_ARG_WITH(metis,
            [AS_HELP_STRING([--with-metis=PATH],
                            [specify prefix directory for METIS])],
            [if test "x$withval" = "x"; then
               with_metis=yes
             fi],
            [with_metis=no])

AC_ARG_WITH(metis-include,
            [AS_HELP_STRING([--with-metis-include=PATH],
                            [specify directory for METIS include files])],
            [if test "x$with_metis" = "xcheck" -o "x$with_metis" = "xno"; then
               with_metis=yes
             fi
             METIS_CPPFLAGS="-I$with_metis_include"],
            [if test "x$with_metis" != "xno" ; then
               if test "x$with_metis" != "xyes" \
	               -a "x$with_metis" != "xcheck"; then
                 METIS_CPPFLAGS="-I$with_metis/include"
               else
                 METIS_CPPFLAGS="-I/usr/include"
               fi
             fi])

AC_ARG_WITH(metis-lib,
            [AS_HELP_STRING([--with-metis-lib=PATH],
                            [specify directory for METIS library])],
            [if test "x$with_metis" = "xcheck" -o "x$with_metis" = "xno"; then
               with_metis=yes
             fi
             METIS_LDFLAGS="-L$with_metis_lib"],
            [if test "x$with_metis" != "xno" -a "x$with_metis" != "xyes" \
	          -a "x$with_metis" != "xcheck"; then
               METIS_LDFLAGS="-L$with_metis/lib"
             fi])


if test "x$with_metis" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # Test for ParMetis first

  CPPFLAGS="${CPPFLAGS} ${METIS_CPPFLAGS} ${MPI_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${METIS_LDFLAGS} ${MPI_LDFLAGS}"
  METIS_LIBS="-lparmetis -lmetis"
  LIBS="${LIBS} ${METIS_LIBS} ${MPI_LIBS}"

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <mpi.h>
#include <parmetis.h>]],
[[  MPI_Comm comm = MPI_COMM_WORLD;
  ParMETIS_V3_PartKway((void *)0, (void *)0, (void *)0, (void *)0,     
                      (void *)0, (void *)0, (void *)0, (void *)0, (void *)0,
                      (void *)0, (void *)0, (void *)0, (void *)0, (void *)0,
                      &comm); ]])],
[have_parmetis=yes],
[have_parmetis=no])

  # Test for METIS second

  if test "x$have_parmetis" = "xno"; then

    METIS_LIBS="-lmetis -lm"
    CPPFLAGS="${saved_CPPFLAGS} ${METIS_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${METIS_LDFLAGS}"
    LIBS="${savedLIBS} ${METIS_LIBS}"

    AC_CHECK_HEADERS([metis.h],
                     [], 
                     [ AC_MSG_WARN([METIS header not found or usable])
                     ],
                     []
                    )

    AC_CHECK_LIB(metis, METIS_PartGraphKway, 
                 [ AC_DEFINE([HAVE_METIS], 1, [use METIS ])
                   have_metis=yes
                 ], 
                 [ AC_MSG_WARN([do not use METIS])
                 ],
                 )

    if test "x$have_metis" = "xno"; then
      METIS_CPPFLAGS=""
      METIS_LDFLAGS=""
      METIS_LIBS=""
    fi

  fi
fi

CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"
LIBS="$saved_LIBS"

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

if test "x$have_parmetis" = "xyes"; then
  AC_DEFINE([HAVE_PARMETIS], 1, [use ParMetis])
elif test "x$have_metis" = "xyes"; then
  AC_DEFINE([HAVE_METIS], 1, [use METIS])
else
  METIS_CPPFLAGS=""
  METIS_LDFLAGS=""
  METIS_LIBS=""
fi

AC_SUBST(METIS_CPPFLAGS)
AC_SUBST(METIS_LDFLAGS)
AC_SUBST(METIS_LIBS)

])dnl


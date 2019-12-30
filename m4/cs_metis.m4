dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2019 EDF S.A.
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

# CS_AC_TEST_METIS
#-----------------
# modifies or sets cs_have_metis, METIS_CPPFLAGS, METIS_LDFLAGS, and METIS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_METIS], [

cs_have_parmetis=no
cs_have_metis=no
metis_prefix=""

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
               fi
             fi])

AC_ARG_WITH(metis-lib,
            [AS_HELP_STRING([--with-metis-lib=PATH],
                            [specify directory for METIS library])],
            [if test "x$with_metis" = "xcheck" -o "x$with_metis" = "xno"; then
               with_metis=yes
             fi
             METIS_LDFLAGS="-L$with_metis_lib"
             # Add the libdir to the runpath as METIS is not libtoolized
             METISRUNPATH="-R$with_metis_lib"],
            [if test "x$with_metis" != "xno" -a "x$with_metis" != "xyes" \
	          -a "x$with_metis" != "xcheck"; then
               METIS_LDFLAGS="-L$with_metis/lib"
               # Add the libdir to the runpath as METIS is not libtoolized
               METISRUNPATH="-R$with_metis/lib"
             fi])


if test "x$with_metis" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${METIS_CPPFLAGS} ${MPI_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${METIS_LDFLAGS} ${MPI_LDFLAGS}"
  METIS_LIBS="-lparmetis -lmetis -lm"
  LIBS="${METIS_LIBS} ${MPI_LIBS} ${LIBS}"

  # Test for ParMetis first

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <mpi.h>
#include <parmetis.h>]],
[[#if PARMETIS_MAJOR_VERSION < 4
# error ParMETIS 4.0 or above required.
#endif
  MPI_Comm comm = MPI_COMM_WORLD;
  ParMETIS_V3_PartKway((void *)0, (void *)0, (void *)0, (void *)0,
                      (void *)0, (void *)0, (void *)0, (void *)0, (void *)0,
                      (void *)0, (void *)0, (void *)0, (void *)0, (void *)0,
                      &comm); ]])],
[cs_have_parmetis=yes
 cs_have_metis=yes],
[cs_have_parmetis=no])

  # Test for METIS second

  if test "x$cs_have_parmetis" = "xno"; then

    METIS_LIBS="-lmetis -lm"
    CPPFLAGS="${saved_CPPFLAGS} ${METIS_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${METIS_LDFLAGS}"
    LIBS="${METIS_LIBS} ${savedLIBS}"

    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <metis.h>]],
[[#if METIS_VER_MAJOR < 5
# error METIS 5.0 or above required.
#endif
  METIS_PartGraphKway((void *)0, (void *)0, (void *)0,
                      (void *)0, (void *)0, (void *)0, (void *)0,
                      (void *)0, (void *)0, (void *)0, (void *)0,
                      (void *)0, (void *)0); ]])],
[cs_have_metis=yes],
[cs_have_metis=no])

    if test "x$cs_have_metis" = "xno"; then
      METIS_CPPFLAGS=""
      METIS_LDFLAGS=""
      METIS_LIBS=""
    fi

  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  case $host_os in
    mingw64)
      metis_prefix=`cygpath --path --windows "$with_metis"`;;
    *)
      ;;
  esac
fi

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

if test "x$cs_have_parmetis" = "xyes"; then
  AC_DEFINE([HAVE_PARMETIS], 1, [use ParMetis])
elif test "x$cs_have_metis" = "xyes"; then
  AC_DEFINE([HAVE_METIS], 1, [use METIS])
else
  METIS_CPPFLAGS=""
  METIS_LDFLAGS=""
  METIS_LIBS=""
fi

AC_SUBST(cs_have_metis)
AC_SUBST(metis_prefix, [${metis_prefix}])
AC_SUBST(METIS_CPPFLAGS)
AC_SUBST(METIS_LDFLAGS)
AC_SUBST(METIS_LIBS)
AC_SUBST(METISRUNPATH)

])dnl


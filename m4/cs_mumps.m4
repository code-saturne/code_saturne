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

# CS_AC_TEST_MUMPS
#----------------
# modifies or sets cs_have_mumps, MUMPS_CPPFLAGS, MUMPS_LDFLAGS, and MUMPS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_MUMPS], [

cs_have_mumps=no
cs_have_mumps_header=no
mumps_prefix=""
cs_abs_srcdir=`cd $srcdir && pwd`

AC_ARG_WITH(mumps,
            [AS_HELP_STRING([--with-mumps=PATH],
                            [specify prefix directory for MUMPS])],
            [if test "x$withval" = "x"; then
               with_mumps=no
             fi],
            [with_mumps=no])

if test "x$with_mumps" != "xno" ; then

  if test -f ${with_mumps}/Makefile.inc ; then
    MUMPS_CPPFLAGS=$(make -s -f "$cs_abs_srcdir/build-aux/mumps.makefile" topdir="${with_mumps}" getincludedirs)
    MUMPS_LDFLAGS=$(make -s -f "$cs_abs_srcdir/build-aux/mumps.makefile" topdir="${with_mumps}" getlibdirs)
    MUMPS_LIBS=$(make -s -f "$cs_abs_srcdir/build-aux/mumps.makefile" topdir="${with_mumps}" getlinklibs)
  fi

  MUMPS="${with_mumps}"
  MUMPS_CPPFLAGS="-I${with_mumps}/include ${MUMPS_CPPFLAGS}"
  MUMPS_LDFLAGS="-L${with_mumps}/lib ${MUMPS_LDFLAGS}"
  MUMPS_LIBS="-ldmumps -lsmumps -lmumps_common ${MUMPS_LIBS}"
  if test "x$FC" = "xifort" ; then
    MUMPS_LIBS="${MUMPS_LIBS} -lifcore -lm"
  elif test "x$FC" = "xmpiifort" ; then
    MUMPS_LIBS="${MUMPS_LIBS} -lifcore -lm"
  else
    MUMPS_LIBS="${MUMPS_LIBS} -lgfortran -lm"
  fi
  MUMPSRUNPATH="-R${with_mumps}/lib"

  AC_MSG_NOTICE([MUMPS_CPP=${MUMPS_CPPFLAGS}])
  AC_MSG_NOTICE([MUMPS_LD=${MUMPS_LDFLAGS}])
  AC_MSG_NOTICE([MUMPS_LIBS=${MUMPS_LIBS}])
  AC_MSG_NOTICE([MUMPSRUNPATH=${MUMPSRUNPATH}])

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${MUMPS_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${MUMPS_LDFLAGS}"
  LIBS="${LIBS} ${MUMPS_LIBS}"

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <dmumps_c.h>]],
[[DMUMPS_STRUC_C id;id.job=-1;dmumps_c(&id);id.job=-2;dmumps_c(&id);]]
[[id.job=-1;id.par=1;id.sym=0;dmumps_c(&id);]]
[[id.job=-2;dmumps_c(&id);]])
                   ],
                   [ AC_DEFINE([HAVE_MUMPS], 1, [Mumps double-precision support])
                     cs_have_mumps=yes
                   ],
                   [ AC_MSG_WARN([no double-precision Mumps support])
                     cs_have_mumps=no
                   ],
                  )

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <smumps_c.h>]],
[[SMUMPS_STRUC_C id;id.job=-1;smumps_c(&id);id.job=-2;smumps_c(&id);]]
[[id.job=-1;id.par=1;id.sym=0;smumps_c(&id);]]
[[id.job=-2;smumps_c(&id);]])
                   ],
                   [ AC_DEFINE([HAVE_MUMPS], 1, [Mumps single-precision support])
                     cs_have_mumps=yes
                   ],
                   [ AC_MSG_WARN([no single-precision Mumps support])
                     cs_have_mumps=no
                   ],
                  )

  if test "x$cs_have_mumps" = "xno"; then
    MUMPS_CPPFLAGS=""
    MUMPS_LDFLAGS=""
    MUMPS_LIBS=""
    MUMPSRUNPATH=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MUMPS, test x$cs_have_mumps = xyes)

AC_SUBST(cs_have_mumps)
AC_SUBST(MUMPS_CPPFLAGS)
AC_SUBST(MUMPS_LDFLAGS)
AC_SUBST(MUMPS_LIBS)
AC_SUBST(MUMPSRUNPATH)

])dnl

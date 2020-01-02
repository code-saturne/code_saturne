dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
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

# CS_AC_TEST_FREESTEAM
#---------------------
# modifies or sets cs_have_freesteam, FREESTEAM_CPPFLAGS, FREESTEAM_LDFLAGS,
# and FREESTEAM_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_FREESTEAM], [

cs_have_freesteam=no
cs_have_freesteam_headers=no

AC_ARG_WITH(freesteam,
            [AS_HELP_STRING([--with-freesteam=DIR],
                            [specify prefix directory for Freesteam])],
            [if test "x$withval" = "x"; then
               with_freesteam=yes
             fi],
            [with_freesteam=check])

AC_ARG_WITH(freesteam-include,
            [AS_HELP_STRING([--with-freesteam-include=DIR],
                            [specify directory for freesteam include files])],
            [if test "x$with_freesteam" = "xcheck"; then
               with_freesteam=yes
             fi
             FREESTEAM_CPPFLAGS="-I$with_freesteam_include"],
            [if test "x$with_freesteam" != "xno" -a "x$with_freesteam" != "xyes" \
	          -a "x$with_freesteam" != "xcheck"; then
               FREESTEAM_CPPFLAGS="-I$with_freesteam/include"
             fi])

AC_ARG_WITH(freesteam-lib,
            [AS_HELP_STRING([--with-freesteam-lib=DIR],
                            [specify directory for freesteam library])],
            [if test "x$with_freesteam" = "xcheck"; then
               with_freesteam=yes
             fi
             FREESTEAM_LDFLAGS="-L$with_freesteam_lib"
             # Add the libdir to the runpath as freesteam is not libtoolized
             FREESTEAMRUNPATH="-R$with_freesteam_lib"],
            [if test "x$with_freesteam" != "xno" -a "x$with_freesteam" != "xyes" \
	          -a "x$with_freesteam" != "xcheck"; then
               FREESTEAM_LDFLAGS="-L$with_freesteam/lib"
               # Add the libdir to the runpath as freesteam is not libtoolized
               FREESTEAMRUNPATH="-R$with_freesteam/lib"
             fi])

if test "x$with_freesteam" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  FREESTEAM_LIBS="-lfreesteam"
  CPPFLAGS="${CPPFLAGS} ${FREESTEAM_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${FREESTEAM_LDFLAGS}"
  LIBS="${FREESTEAM_LIBS} ${LIBS}"


  # Check for steam.h header
  AC_CHECK_HEADERS([freesteam/steam.h],
                   [AC_DEFINE([HAVE_FREESTEAM], 1, [freesteam support])
                    cs_have_freesteam=yes],
                   [cs_have_freesteam=no],
                   [])

  if test "x$cs_have_freesteam" != "xyes"; then
    FREESTEAM_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_FREESTEAM, test x$cs_have_freesteam = xyes)

AC_SUBST(cs_have_freesteam)
AC_SUBST(freesteam_prefix, [${with_freesteam}])
AC_SUBST(FREESTEAM_CPPFLAGS)
AC_SUBST(FREESTEAM_LDFLAGS)
AC_SUBST(FREESTEAM_LIBS)
AC_SUBST(FREESTEAMRUNPATH)

])dnl


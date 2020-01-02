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

# CS_AC_TEST_EOS
#---------------
# modifies or sets cs_have_eos, EOS_CPPFLAGS, EOS_LDFLAGS, and EOS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_EOS], [

cs_have_eos=no
cs_have_eos_headers=no

AC_ARG_WITH(eos,
            [AS_HELP_STRING([--with-eos=DIR],
                            [specify prefix directory for EOS])],
            [if test "x$withval" = "x"; then
               with_eos=yes
             fi],
            [with_eos=no])

AC_ARG_WITH(eos-include,
            [AS_HELP_STRING([--with-eos-include=DIR],
                            [specify directory for EOS include files])],
            [if test "x$with_eos" = "xcheck"; then
               with_eos=yes
             fi
             EOS_CPPFLAGS="-I$with_eos_include"],
            [if test "x$with_eos" != "xno" -a "x$with_eos" != "xyes" \
	          -a "x$with_eos" != "xcheck"; then
               EOS_CPPFLAGS="-I$with_eos/include"
             fi])

AC_ARG_WITH(eos-lib,
            [AS_HELP_STRING([--with-eos-lib=DIR],
                            [specify directory for EOS library])],
            [if test "x$with_eos" = "xcheck"; then
               with_eos=yes
             fi
             EOS_LDFLAGS="-L$with_eos_lib"
             # Add the libdir to the runpath as EOS is not libtoolized
             EOSRUNPATH="-R$with_eos_lib"],
            [if test "x$with_eos" != "xno" -a "x$with_eos" != "xyes" \
	          -a "x$with_eos" != "xcheck"; then
               EOS_LDFLAGS="-L$with_eos/lib"
               # Add the libdir to the runpath as EOS is not libtoolized
               EOSRUNPATH="-R$with_eos/lib"
             fi])

if test "x$with_eos" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # Now check library

  EOS_LIBS="-leos"

  eosversion=`${with_eos}/bin/eos --version`

  if test "$eosversion" = "1.3.0"; then
    EOS_LIBS="-lCCLanguageAPI -lCCEOSAPI"
    EOS_CPPFLAGS="${EOS_CPPFLAGS} -DEOS_PRE_V1_6"
  fi
  if test "$eosversion" = "1.4.0"; then
    EOS_LIBS="-lshared_eos"
    EOS_CPPFLAGS="${EOS_CPPFLAGS} -DEOS_PRE_V1_6"
  fi
  if test "$eosversion" = "1.5.0"; then
    EOS_LIBS="-lshared_eos"
    EOS_CPPFLAGS="${EOS_CPPFLAGS} -DEOS_PRE_V1_6"
  fi

  CPPFLAGS="${CPPFLAGS} ${EOS_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${EOS_LDFLAGS}"
  LIBS="${EOS_LIBS} ${LIBS}"

  # Check that EOS files exist
  AC_LANG_PUSH([C++])

  AC_MSG_CHECKING([for EOS library)])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "EOS/API/EOS.hxx"]],
                                  [[NEPTUNE::EOS *eos]])],
                                  [ AC_DEFINE([HAVE_EOS], 1, [EOS support])
                                    cs_have_eos=yes],
                                   [cs_have_eos=no])
  AC_MSG_RESULT($cs_have_eos)

  if test "x$cs_have_eos" = "xno" ; then
    if test "x$with_eos" != "xcheck" ; then
      AC_MSG_FAILURE([EOS support is requested, but test for EOS failed!])
    fi
  fi

  if test "x$cs_have_eos" != "xyes"; then
    EOS_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_EOS, test x$cs_have_eos = xyes)

AC_SUBST(cs_have_eos)
AC_SUBST(eos_prefix, [${with_eos}])
AC_SUBST(EOS_CPPFLAGS)
AC_SUBST(EOS_LDFLAGS)
AC_SUBST(EOS_LIBS)
AC_SUBST(EOSRUNPATH)

AC_LANG_POP([C++])

])dnl


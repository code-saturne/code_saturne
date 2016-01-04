dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2016 EDF S.A.
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

# CS_AC_TEST_COOLPROP
#---------------
# modifies or sets cs_have_coolprop, COOLPROP_CPPFLAGS, COOLPROP_LDFLAGS,
# and COOLPROP_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_COOLPROP], [

cs_have_coolprop=no
cs_have_coolprop_headers=no

AC_ARG_WITH(coolprop,
            [AS_HELP_STRING([--with-coolprop=DIR],
                            [specify prefix directory for COOLPROP])],
            [if test "x$withval" = "x"; then
               with_coolprop=yes
             fi],
            [with_coolprop=check])

AC_ARG_WITH(coolprop-include,
            [AS_HELP_STRING([--with-coolprop-include=DIR],
                            [specify directory for CoolProp include files])],
            [if test "x$with_coolprop" = "xcheck"; then
               with_coolprop=yes
             fi
             COOLPROP_CPPFLAGS="-I$with_coolprop_include"],
            [if test "x$with_coolprop" != "xno" -a "x$with_coolprop" != "xyes" \
	          -a "x$with_coolprop" != "xcheck"; then
               if test -f "$with_coolprop/include/CoolPropLib.h"; then
                 COOLPROP_CPPFLAGS="-I$with_coolprop/include"
               elif test -f "$with_coolprop/CoolPropLib.h"; then
                 COOLPROP_CPPFLAGS="-I$with_coolprop"
               fi
             fi])

AC_ARG_WITH(coolprop-lib,
            [AS_HELP_STRING([--with-coolprop-lib=DIR],
                            [specify directory for CoolProp library])],
            [if test "x$with_coolprop" = "xcheck"; then
               with_coolprop=yes
             fi
             COOLPROP_LDFLAGS="-L$with_coolprop_lib"
             # Add the libdir to the runpath as CoolProp is not libtoolized
             COOLPROPRUNPATH="-R$with_coolprop_lib"],
            [if test "x$with_coolprop" != "xno" -a "x$with_coolprop" != "xyes" \
	          -a "x$with_coolprop" != "xcheck"; then
               case `uname -m` in
                  *64)        ref_name='64bit/libCoolProp' ;;
                  *32 | i*86) ref_name='32bit/libCoolProp' ;;
                  *)          ref_name='libCoolProp' ;;
               esac
               find $with_coolprop | grep ${ref_name}.* > /dev/null 2>&1;
               if test $? == "0"; then
                 cp_l_name=`find $with_coolprop | grep ${ref_name}.* | head -1`
                 cp_d_name=`dirname ${cp_l_name}`
                 unset cp_l_name
                 COOLPROP_LDFLAGS="-L$cp_d_name"
                 # Add the libdir to the runpath as CoolProp is not libtoolized
                 COOLPROPRUNPATH="-R$cp_d_name"
                 unset cp_d_name
               else
                 COOLPROP_LDFLAGS="-L$with_coolprop/lib"
                 # Add the libdir to the runpath as CoolProp is not libtoolized
                 COOLPROPRUNPATH="-R$with_coolprop/lib"
               fi
             fi])

if test "x$with_coolprop" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  COOLPROP_LIBS="-lCoolProp"
  CPPFLAGS="${CPPFLAGS} ${COOLPROP_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${COOLPROP_LDFLAGS}"
  LIBS="${COOLPROP_LIBS} ${LIBS}"

  # Check that CoolProp files exist

  # Coolprop is in C++, but provides a C wrapper, which we use for determination here.
  # AC_LANG_PUSH([C++])

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "CoolPropLib.h"]],
                                  [[double v = PropsSI("o", "n1", 1., "n2", 1., "r");]])],
                                  [ AC_DEFINE([HAVE_COOLPROP], 1, [CoolProp support])
                                    cs_have_coolprop=yes],
                                   [cs_have_coolprop=no])

  if test "x$cs_have_coolprop" != "xyes"; then
    COOLPROP_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_COOLPROP, test x$cs_have_coolprop = xyes)

AC_SUBST(cs_have_coolprop)
AC_SUBST(coolprop_prefix, [${with_coolprop}])
AC_SUBST(COOLPROP_CPPFLAGS)
AC_SUBST(COOLPROP_LDFLAGS)
AC_SUBST(COOLPROP_LIBS)
AC_SUBST(COOLPROPRUNPATH)

# AC_LANG_POP([C++])

])dnl


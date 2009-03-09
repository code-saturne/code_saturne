dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2009 EDF S.A., France
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

# CS_AC_TEST_BFT(Minimal Release string, [Maximal Release string])
#-----------------------------------------------------------------
# Check for BFT version ; defines BFT_CPPFLAGS, BFT_LDFLAGS, and BFT_LIBS
# locally (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CS_AC_TEST_BFT], [

AC_ARG_WITH(bft, [AS_HELP_STRING([--with-bft=PATH], [specify prefix directory for BFT])])
AC_ARG_WITH(bft-exec, [AS_HELP_STRING([--with-bft-exec=PATH], [specify directory for BFT executables])])
AC_ARG_WITH(bft-include, [AS_HELP_STRING([--with-bft-include=PATH], [specify directory for BFT include files])])
AC_ARG_WITH(bft-lib, [AS_HELP_STRING([--with-bft-lib=PATH], [specify directory for BFT library])])

if test "x$with_bft_exec" != "x" ; then
  bft_config="$with_bft_exec/bft-config"
elif test "x$with_bft" != "x" ; then
  bft_config="$with_bft/bin/bft-config"
else
  bft_config="bft-config"
fi

if test "x$with_bft_include" != "x" ; then
  BFT_CPPFLAGS="-I$with_bft_include"
elif test "x$with_bft" != "x" ; then
  BFT_CPPFLAGS="-I$with_bft/include"
else
  BFT_CPPFLAGS=""
fi

if test "x$with_bft_lib" != "x" ; then
  BFT_LDFLAGS="-L$with_bft_lib"
elif test "x$with_bft" != "x" ; then
  BFT_LDFLAGS="-L$with_bft/lib"
else
  BFT_LDFLAGS=""
fi
BFT_LIBS="-lbft"

type "$bft_config" > /dev/null 2>&1
if test "$?" = "0" ; then
  BFT_CPPFLAGS="$BFT_CPPFLAGS `$bft_config --cppflags`"
  BFT_LDFLAGS="$BFT_LDFLAGS `$bft_config --ldflags`"
  BFT_LIBS="$BFT_LIBS `$bft_config --libs`"
fi

bft_version_min=$1
bft_version_max=$2

if test "x$bft_version_min" != "x" ; then
  if test "x$bft_version_max" != "x" ; then
    AC_MSG_CHECKING([for bft version >= $1 and <= $2])
  else
    AC_MSG_CHECKING([for bft version >= $1])
  fi
else
  bft_version_min="0.0.0"
fi

bft_version_major_min=`echo "$bft_version_min" | cut -f1 -d.`
bft_version_minor_min=`echo "$bft_version_min" | cut -f2 -d.`
bft_version_release_min=`echo "$bft_version_min" | cut -f3 -d.`
if test    "$bft_version_major_min" = "" \
        -o "$bft_version_minor_min" = "" \
        -o "$bft_version_release_min" = ""; then
  AC_MSG_FAILURE([bad BFT version definition in configure.ac: $bft_version_min])
fi

if test "x$bft_version_max" != "x" ; then
  bft_version_major_max=`echo "$bft_version_max" | cut -f1 -d.`
  bft_version_minor_max=`echo "$bft_version_max" | cut -f2 -d.`
  bft_version_release_max=`echo "$bft_version_max" | cut -f3 -d.`
  if test    "$bft_version_major_max" = "" \
          -o "$bft_version_minor_max" = "" \
          -o "$bft_version_release_max" = ""; then
    AC_MSG_FAILURE([bad BFT version definition in configure.ac: $bft_version_max])
  fi
else
  bft_version_major_max=99999999
  bft_version_minor_max=99999999
  bft_version_release_max=99999999
fi

saved_CPPFLAGS=$CPPFLAGS
saved_LDFLAGS=$LDFLAGS
saved_LIBS=$LIBS

CPPFLAGS="${CPPFLAGS} $BFT_CPPFLAGS"
LDFLAGS="${LDFLAGS} $BFT_LDFLAGS"
LIBS="${LIBS} $BFT_LIBS"

AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <bft_config.h>
]],
[[#if BFT_MAJOR_VERSION < $bft_version_major_min
#  error BFT major version < $bft_version_major_min
#elif BFT_MAJOR_VERSION == $bft_version_major_min
#  if BFT_MINOR_VERSION < $bft_version_minor_min
#    error BFT minor version < $bft_version_minor_min
#  elif BFT_MINOR_VERSION == $bft_version_minor_min
#    if BFT_RELEASE_VERSION < $bft_version_release_min
#      error BFT release version < $bft_version_release_min
#    endif
#  endif
#endif
#if BFT_MAJOR_VERSION > $bft_version_major_max
#  error BFT major version > $bft_version_major_max
#elif BFT_MAJOR_VERSION == $bft_version_major_max
#  if BFT_MINOR_VERSION > $bft_version_minor_max
#    error BFT minor version > $bft_version_minor_max
#  elif BFT_MINOR_VERSION == $bft_version_minor_max
#    if BFT_RELEASE_VERSION > $bft_version_release_max
#      error BFT release version < $bft_version_release_max
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible bft version found])],
               [AC_MSG_FAILURE([compatible bft version not found])])

unset bft_version_major_min
unset bft_version_minor_min
unset bft_version_release_min
unset bft_version_major_max
unset bft_version_minor_max
unset bft_version_release_max

unset bft_version_min
unset bft_version_max
unset bft_config

# Restore old LIBS to add $BFT_LIBS later, as other tests
# might otherwise not run if a shared library is not found

CPPFLAGS=$saved_CPPFLAGS
LDFLAGS=$saved_LDFLAGS
LIBS=$saved_LIBS

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(BFT_CPPFLAGS)
AC_SUBST(BFT_LDFLAGS)
AC_SUBST(BFT_LIBS)

])dnl

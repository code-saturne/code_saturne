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

# CS_AC_TEST_PLE(Minimal Release string, [Maximal Release string])
#-----------------------------------------------------------------
# Check for PLE version ; defines PLE_CPPFLAGS, PLE_LDFLAGS, and PLE_LIBS
# locally (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CS_AC_TEST_PLE], [

cs_have_internal_ple=yes

AC_ARG_WITH(ple,
            [AS_HELP_STRING([--with-ple=PATH],
                            [specify prefix directory for PLE])],
            [if test "x$withval" = "x"; then
               with_ple=yes
             fi],
            [with_ple=no])

AC_ARG_WITH(ple-exec,
            [AS_HELP_STRING([--with-ple-exec=PATH],
                            [specify directory for PLE executables])],
            [ple_config="$with_ple_exec/ple-config"],
            [if test "x$with_ple" != "xyes"; then
               ple_config="$with_ple/bin/ple-config"
             else
               ple_config="ple-config"
             fi])

AC_ARG_WITH(ple-include,
            [AS_HELP_STRING([--with-ple-include=PATH],
                            [specify directory for PLE include files])],
            [PLE_CPPFLAGS="-I$with_ple_include"],
            [if test "x$with_ple" != "xno" -a "x$with_ple" != "xyes"; then
               PLE_CPPFLAGS="-I$with_ple/include"
             fi])

AC_ARG_WITH(ple-lib,
            [AS_HELP_STRING([--with-ple-lib=PATH],
                            [specify directory for PLE library])],
            [PLE_LDFLAGS="-L$with_ple_lib"],
            [if test "x$with_ple" != "xno" -a "x$with_ple" != "xyes"; then
               PLE_LDFLAGS="-L$with_ple/lib"
             fi])

AC_ARG_WITH(ple-doc,
            [AS_HELP_STRING([--with-ple-doc=PATH],
                            [specify directory for PLE documentation])],
            [ple_docdir="$with_ple_doc"],
            [if test "x$with_ple" = "xno" ; then
               ple_docdir="$datarootdir/doc/ple"
             elif test "x$with_ple" != "xyes"; then
               ple_docdir="$with_ple/share/doc/ple"
             fi])

PLE_LIBS="-lple"

if test "x$with_ple" != "xno" ; then

  type "$ple_config" > /dev/null 2>&1
  if test "$?" = "0" ; then
    PLE_CPPFLAGS="$PLE_CPPFLAGS `$ple_config --cppflags`"
    PLE_LDFLAGS="$PLE_LDFLAGS `$ple_config --ldflags`"
    PLE_LIBS="$PLE_LIBS `$ple_config --libs`"
    PLE_MPI_LDFLAGS="`$ple_config --ldflags mpi`"
    PLE_MPI_LIBS="`$ple_config --libs mpi`"
  fi

  ple_version_min=$1
  ple_version_max=$2

  if test "x$ple_version_min" != "x" ; then
    if test "x$ple_version_max" != "x" ; then
      AC_MSG_CHECKING([for ple version >= $1 and <= $2])
    else
      AC_MSG_CHECKING([for ple version >= $1])
    fi
  else
    ple_version_min="0.0.0"
  fi

  ple_version_major_min=`echo "$ple_version_min" | cut -f1 -d.`
  ple_version_minor_min=`echo "$ple_version_min" | cut -f2 -d.`
  ple_version_release_min=`echo "$ple_version_min" | cut -f3 -d.`
  if test    "$ple_version_major_min" = "" \
          -o "$ple_version_minor_min" = "" \
          -o "$ple_version_release_min" = ""; then
    AC_MSG_FAILURE([bad PLE version definition in configure.ac: $ple_version_min])
  fi

  if test "x$ple_version_max" != "x" ; then
    ple_version_major_max=`echo "$ple_version_max" | cut -f1 -d.`
    ple_version_minor_max=`echo "$ple_version_max" | cut -f2 -d.`
    ple_version_release_max=`echo "$ple_version_max" | cut -f3 -d.`
    if test    "$ple_version_major_max" = "" \
            -o "$ple_version_minor_max" = "" \
            -o "$ple_version_release_max" = ""; then
      AC_MSG_FAILURE([bad PLE version definition in configure.ac: $ple_version_max])
    fi
  else
    ple_version_major_max=99999999
    ple_version_minor_max=99999999
    ple_version_release_max=99999999
  fi

  saved_CPPFLAGS=$CPPFLAGS
  saved_LDFLAGS=$LDFLAGS
  saved_LIBS=$LIBS

  CPPFLAGS="${CPPFLAGS} $PLE_CPPFLAGS"
  LDFLAGS="$PLE_LDFLAGS $PLE_MPI_LDFLAGS ${LDFLAGS}"
  LIBS="$PLE_LIBS $PLE_MPI_LIBS ${LIBS}"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <ple_config.h>
  ]],
  [[#if PLE_MAJOR_VERSION < $ple_version_major_min
#  error PLE major version < $ple_version_major_min
#elif PLE_MAJOR_VERSION == $ple_version_major_min
#  if PLE_MINOR_VERSION < $ple_version_minor_min
#    error PLE minor version < $ple_version_minor_min
#  elif PLE_MINOR_VERSION == $ple_version_minor_min
#    if PLE_RELEASE_VERSION < $ple_version_release_min
#      error PLE release version < $ple_version_release_min
#    endif
#  endif
#endif
#if PLE_MAJOR_VERSION > $ple_version_major_max
#  error PLE major version > $ple_version_major_max
#elif PLE_MAJOR_VERSION == $ple_version_major_max
#  if PLE_MINOR_VERSION > $ple_version_minor_max
#    error PLE minor version > $ple_version_minor_max
#  elif PLE_MINOR_VERSION == $ple_version_minor_max
#    if PLE_RELEASE_VERSION > $ple_version_release_max
#      error PLE release version < $ple_version_release_max
#    endif
#  endif
#endif
  ]])],
                 [cs_have_internal_ple=no
                  AC_MSG_RESULT([using external PLE library])],
                 [AC_MSG_RESULT([external PLE library not found, using internal])])

  CPPFLAGS=$saved_CPPFLAGS
  LDFLAGS=$saved_LDFLAGS
  LIBS=$saved_LIBS

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  ple_type="external"

else

  ple_type="internal"

fi

AC_SUBST(ple_type)
AC_SUBST(ple_docdir)
AC_SUBST(PLE_CPPFLAGS)
AC_SUBST(PLE_LDFLAGS)
AC_SUBST(PLE_LIBS)

])dnl


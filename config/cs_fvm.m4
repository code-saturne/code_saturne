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

# CS_AC_TEST_FVM(Minimal Release string, [Maximal Release string])
#-----------------------------------------------------------------
# Check for FVM version ; defines FVM_CPPFLAGS, FVM_LDFLAGS, and FVM_LIBS
# locally (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CS_AC_TEST_FVM], [

AC_ARG_WITH(fvm, [AS_HELP_STRING([--with-fvm=PATH], [specify prefix directory for FVM])])
AC_ARG_WITH(fvm-exec, [AS_HELP_STRING([--with-fvm-exec=PATH], [specify directory for FVM executables])])
AC_ARG_WITH(fvm-include, [AS_HELP_STRING([--with-fvm-include=PATH], [specify directory for FVM include files])])
AC_ARG_WITH(fvm-lib, [AS_HELP_STRING([--with-fvm-lib=PATH], [specify directory for FVM library])])

if test "x$with_fvm_exec" != "x" ; then
  fvm_config="$with_fvm_exec/fvm-config"
elif test "x$with_fvm" != "x" ; then
  fvm_config="$with_fvm/bin/fvm-config"
else
  fvm_config="fvm-config"
fi

if test "x$with_fvm_include" != "x" ; then
  FVM_CPPFLAGS="-I$with_fvm_include"
elif test "x$with_fvm" != "x" ; then
  FVM_CPPFLAGS="-I$with_fvm/include"
else
  FVM_CPPFLAGS=""
fi

if test "x$with_fvm_lib" != "x" ; then
  FVM_LDFLAGS="-L$with_fvm_lib"
elif test "x$with_fvm" != "x" ; then
  FVM_LDFLAGS="-L$with_fvm/lib"
else
  FVM_LDFLAGS=""
fi
FVM_LIBS="-lfvm"

type "$fvm_config" > /dev/null 2>&1
if test "$?" = "0" ; then
  FVM_CPPFLAGS="$FVM_CPPFLAGS `$fvm_config --cppflags`"
  FVM_DEP_LDFLAGS="`$fvm_config --ldflags`"
  FVM_DEP_LDFLAGS="$FVM_DEP_LDFLAGS `$fvm_config --ldflags cgns`"
  FVM_DEP_LDFLAGS="$FVM_DEP_LDFLAGS `$fvm_config --ldflags hdf5`"
  FVM_DEP_LDFLAGS="$FVM_DEP_LDFLAGS `$fvm_config --ldflags med`"
  FVM_DEP_LDFLAGS="$FVM_DEP_LDFLAGS `$fvm_config --ldflags mpi`"
  FVM_DEP_LIBS="`$fvm_config --libs`"
  FVM_DEP_LIBS="$FVM_DEP_LIBS `$fvm_config --libs cgns`"
  FVM_DEP_LIBS="$FVM_DEP_LIBS `$fvm_config --libs hdf5`"
  FVM_DEP_LIBS="$FVM_DEP_LIBS `$fvm_config --libs med`"
  FVM_DEP_LIBS="$FVM_DEP_LIBS `$fvm_config --libs mpi`"
fi

fvm_version_min=$1
fvm_version_max=$2

if test "x$fvm_version_min" != "x" ; then
  if test "x$fvm_version_max" != "x" ; then
    AC_MSG_CHECKING([for fvm version >= $1 and <= $2])
  else
    AC_MSG_CHECKING([for fvm version >= $1])
  fi
else
  fvm_version_min="0.0.0"
fi

fvm_version_major_min=`echo "$fvm_version_min" | cut -f1 -d.`
fvm_version_minor_min=`echo "$fvm_version_min" | cut -f2 -d.`
fvm_version_release_min=`echo "$fvm_version_min" | cut -f3 -d.`
if test    "$fvm_version_major_min" = "" \
        -o "$fvm_version_minor_min" = "" \
        -o "$fvm_version_release_min" = ""; then
  AC_MSG_FAILURE([bad FVM version definition in configure.ac: $fvm_version_min])
fi

if test "x$fvm_version_max" != "x" ; then
  fvm_version_major_max=`echo "$fvm_version_max" | cut -f1 -d.`
  fvm_version_minor_max=`echo "$fvm_version_max" | cut -f2 -d.`
  fvm_version_release_max=`echo "$fvm_version_max" | cut -f3 -d.`
  if test    "$fvm_version_major_max" = "" \
          -o "$fvm_version_minor_max" = "" \
          -o "$fvm_version_release_max" = ""; then
    AC_MSG_FAILURE([bad FVM version definition in configure.ac: $fvm_version_max])
  fi
else
  fvm_version_major_max=99999999
  fvm_version_minor_max=99999999
  fvm_version_release_max=99999999
fi

saved_CPPFLAGS=$CPPFLAGS
saved_LDFLAGS=$LDFLAGS
saved_LIBS=$LIBS

CPPFLAGS="${CPPFLAGS} $FVM_CPPFLAGS"
LDFLAGS="$FVM_LDFLAGS $FVM_DEP_LDFLAGS `$fvm_config --ldflags bft` ${LDFLAGS}"
LIBS="$FVM_LIBS $FVM_DEP_LIBS `$fvm_config --libs bft` ${LIBS}"

AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fvm_config.h>
]],
[[#if FVM_MAJOR_VERSION < $fvm_version_major_min
#  error FVM major version < $fvm_version_major_min
#elif FVM_MAJOR_VERSION == $fvm_version_major_min
#  if FVM_MINOR_VERSION < $fvm_version_minor_min
#    error FVM minor version < $fvm_version_minor_min
#  elif FVM_MINOR_VERSION == $fvm_version_minor_min
#    if FVM_RELEASE_VERSION < $fvm_version_release_min
#      error FVM release version < $fvm_version_release_min
#    endif
#  endif
#endif
#if FVM_MAJOR_VERSION > $fvm_version_major_max
#  error FVM major version > $fvm_version_major_max
#elif FVM_MAJOR_VERSION == $fvm_version_major_max
#  if FVM_MINOR_VERSION > $fvm_version_minor_max
#    error FVM minor version > $fvm_version_minor_max
#  elif FVM_MINOR_VERSION == $fvm_version_minor_max
#    if FVM_RELEASE_VERSION > $fvm_version_release_max
#      error FVM release version < $fvm_version_release_max
#    endif
#  endif
#endif
]])],
               [AC_MSG_RESULT([compatible fvm version found])],
               [AC_MSG_FAILURE([compatible fvm version not found])])

# Now check if fvm_coupl library is available (using MPI)

AC_MSG_CHECKING([for fvm_coupling discovery functions])

LIBS="$FVM_LIBS -lfvm_coupl $FVM_DEP_LIBS `$fvm_config --libs bft` ${saved_LIBS}"

AC_LINK_IFELSE([AC_LANG_PROGRAM([[int fvm_coupling_mpi_world_n_apps(void *);]],
               [[fvm_coupling_mpi_world_n_apps(0); ]])],
                [fvm_have_coupl=yes],
                [fvm_have_coupl=no])

AC_MSG_RESULT($fvm_have_coupl)
if test "$fvm_have_coupl" = "yes"; then
  FVM_LIBS="$FVM_LIBS -lfvm_coupl"
fi
FVM_LDFLAGS="$FVM_LDFLAGS $FVM_DEP_LDFLAGS"
FVM_LIBS="$FVM_LIBS $FVM_DEP_LIBS"

# Unset temporary variables

unset fvm_version_major_min
unset fvm_version_minor_min
unset fvm_version_release_min
unset fvm_version_major_max
unset fvm_version_minor_max
unset fvm_version_release_max

unset fvm_version_min
unset fvm_version_max
unset fvm_config

# Restore old LIBS to add $FVM_LIBS later, as other tests
# might otherwise not run if a shared library is not found

CPPFLAGS=$saved_CPPFLAGS
LDFLAGS=$saved_LDFLAGS
LIBS=$saved_LIBS

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(FVM_CPPFLAGS)
AC_SUBST(FVM_LDFLAGS)
AC_SUBST(FVM_LIBS)

])dnl

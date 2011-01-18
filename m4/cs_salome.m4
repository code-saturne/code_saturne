dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2011 EDF S.A., France
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

# CS_AC_TEST_SALOME_KERNEL
#-------------------------
# modifies or sets cs_have_salome_kernel, SALOME_KERNEL_CPPFLAGS, SALOME_KERNEL_LDFLAGS, and SALOME_KERNEL_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SALOME_KERNEL], [

cs_have_salome=no

AC_ARG_WITH(salome-kernel,
            [AS_HELP_STRING([--with-salome-kernel=PATH],
                            [specify prefix directory for SALOME kernel])],
            [if test "x$withval" = "x"; then
               if test -z "$KERNEL_ROOT_DIR"; then
                 with_salome_kernel=yes
               else
                 with_salome_kernel=$KERNEL_ROOT_DIR
               fi
             fi],
            [if test -z "$KERNEL_ROOT_DIR"; then
               with_salome_kernel=check
             else
               with_salome_kernel=$KERNEL_ROOT_DIR
             fi])

AC_ARG_WITH(salome-kernel-include,
            [AS_HELP_STRING([--with-salome-kernel-include=PATH],
                            [specify directory for SALOME kernel include files])],
            [if test "x$with_salome_kernel" = "xcheck"; then
               with_salome_kernel=yes
             fi
             SALOME_KERNEL_CPPFLAGS="-I$with_salome_kernel_include"],
            [if test "x$with_salome_kernel" != "xno" ; then
               if test "x$with_salome_kernel" != "xyes" \
	               -a "x$with_salome_kernel" != "xcheck"; then
                 SALOME_KERNEL_CPPFLAGS="-I$with_salome_kernel/include/salome"
               else
                 SALOME_KERNEL_CPPFLAGS="-I/usr/include/salome"
               fi
             fi])

AC_ARG_WITH(salome-kernel-lib,
            [AS_HELP_STRING([--with-salome-kernel-lib=PATH],
                            [specify directory for SALOME_KERNEL library])],
            [if test "x$with_salome_kernel" = "xcheck"; then
               with_salome_kernel=yes
             fi
             SALOME_KERNEL_LDFLAGS="-L$with_salome_kernel_lib"],
            [if test "x$with_salome_kernel" != "xno" ; then
               if test "x$with_salome_kernel" != "xyes" \
	               -a "x$with_salome_kernel" != "xcheck"; then
                 SALOME_KERNEL_LDFLAGS="-L$with_salome_kernel/lib/salome"
               else
                 SALOME_KERNEL_LDFLAGS="-L/usr/lib/salome"
               fi
             fi])

if test "x$with_salome_kernel" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  SALOME_KERNEL_LIBS="-lCalciumC"
  
  CPPFLAGS="${CPPFLAGS} ${SALOME_KERNEL_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${SALOME_KERNEL_LDFLAGS}"
  LIBS="${LIBS} ${SALOME_KERNEL_LIBS}"

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <calcium.h>]],
  			             [[int iret = cp_fin(0, 0);]])],
                    [cs_have_salome_kernel=yes
                     AC_MSG_RESULT([compatible SALOME kernel found])],
                    [cs_have_salome_kernel=no
                     if test "x$with_salome_kernel" != "xcheck" ; then
                       AC_MSG_FAILURE([SALOME support is requested, but test for SALOME failed!])
                     else
                       AC_MSG_WARN([no SALOME support])
                     fi
                    ])

  if test "x$cs_have_salome_kernel" = "xno"; then
    SALOME_KERNEL_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(cs_have_salome_kernel)
AC_SUBST(SALOME_KERNEL_CPPFLAGS)
AC_SUBST(SALOME_KERNEL_IDL)
AC_SUBST(SALOME_KERNEL_LDFLAGS)
AC_SUBST(SALOME_KERNEL_LIBS)

AM_CONDITIONAL(HAVE_SALOME_KERNEL, test x$cs_have_salome_kernel = xyes)

])dnl


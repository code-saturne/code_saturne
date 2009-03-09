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

# CS_AC_TEST_BLAS
#----------------
# modifies or sets have_blas, BLAS_CPPFLAGS, BLAS_LDFLAGS, and BLAS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_BLAS], [

have_blas=no

AC_ARG_ENABLE(blas,
  [  --disable-blas          do not use BLAS when available],
  [
    case "${enableval}" in
      yes) blas=true ;;
      no)  blas=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-blas]) ;;
    esac
  ],
  [ blas=true ]
)

AC_ARG_WITH(blas, [AS_HELP_STRING([--with-blas=PATH], [specify prefix directory for BLAS])])
AC_ARG_WITH(blas-include, [AS_HELP_STRING([--with-blas-include=PATH], [specify directory for BLAS include files])])
AC_ARG_WITH(blas-lib, [AS_HELP_STRING([--with-blas-lib=PATH], [specify directory for BLAS library])])

if test "x$blas" = "xtrue" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$with_blas_include" != "x" ; then
    BLAS_CPPFLAGS="-I$with_blas_include -D_CS_HAVE_CBLAS"
    #BLAS_CPPFLAGS="-I$with_blas_include -D_CS_HAVE_ESSL"
    #BLAS_CPPFLAGS="-I$with_blas_include -D_CS_HAVE_MKL"
  elif test "x$with_blas" != "x" ; then
    BLAS_CPPFLAGS="-I$with_blas/include -D_CS_HAVE_CBLAS"
    #BLAS_CPPFLAGS="-I$with_blas/include -D_CS_HAVE_ESSL"
    #BLAS_CPPFLAGS="-I$with_blas/include -D_CS_HAVE_MKL"
  fi

  if test "x$with_blas_lib" != "x" ; then
    BLAS_LDFLAGS="-L$with_blas_lib"
  elif test "x$with_blas" != "x" ; then
    BLAS_LDFLAGS="-L$with_blas/lib"
  fi

  BLAS_LIBS="-lcblas -latlas"
  #BLAS_LIBS="-lesslbg"
  #BLAS_LIBS="-lesslbg -lesslsmpbg"
  #BLAS_LIBS="-lmkl -lmkl_blacs_intelmpi20 -lmkl_ipf -lguide -lpthread"

  CPPFLAGS="${CPPFLAGS} ${BLAS_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${BLAS_LDFLAGS}"
  LIBS="${LIBS} ${BLAS_LIBS}"

  AC_CHECK_LIB(blas, cblas_ddot, 
               [ AC_DEFINE([HAVE_BLAS], 1, [BLAS support])
                 have_blas=yes
               ], 
               [ AC_MSG_WARN([no BLAS support])
               ],
               )

  if test "x$have_blas" != "xyes"; then
    BLAS_CPPFLAGS=""
    BLAS_LDFLAGS=""
    BLAS_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_BLAS, test x$have_blas = xyes)

AC_SUBST(BLAS_CPPFLAGS)
AC_SUBST(BLAS_LDFLAGS)
AC_SUBST(BLAS_LIBS)

])dnl


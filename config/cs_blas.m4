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

# CS_AC_TEST_BLAS([use_threads])
#----------------
# modifies or sets cs_have_blas, BLAS_CPPFLAGS, BLAS_LDFLAGS, and BLAS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_BLAS], [

cs_have_blas=no

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
AC_ARG_WITH(blas-type, [AS_HELP_STRING([--with-blas-type=NAME], [force ATLAS, ESSL, MKL, ...])])
AC_ARG_WITH(blas-libs, [AS_HELP_STRING([--with-blas-libs=LIBS], [specify BLAS libraries])])

if test "x$blas" = "xtrue" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  BLAS_CPPFLAGS=""
  BLAS_LDFLAGS=""

  # Also add known paths and libraries for Blue Gene/L or P if not given

  if test "x$with_blas_include" != "x" ; then
    BLAS_CPPFLAGS="-I$with_blas_include" 
  elif test "x$with_blas" != "x" ; then
    BLAS_CPPFLAGS="-I$with_blas/include" 
  fi

  if test "x$with_blas_lib" != "x" ; then
    BLAS_LDFLAGS="-L$with_blas_lib"
  elif test "x$with_blas" != "x" ; then
    BLAS_LDFLAGS="-L$with_blas/lib" 
  fi

  if test "x$with_blas_type" = "x" ; then
    if test "x$cs_ibm_bg_type" != "x" ; then
      with_blas_type="ESSL"
    fi
  fi

  # Test for IBM ESSL BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xESSL" ; then

    # Test compilation/header separately from link, as linking may require
    # Fortran libraries, and header is for C. Test library (link) first,
    # as header is only useful if library is present.

    if test "x$BLAS_CPPFLAGS" = "x" ; then
      if test "x$cs_ibm_bg_type" = "xL" ; then
        BLAS_CPPFLAGS="-I/opt/ibmmath/essl/4.2/include"
      elif test "x$cs_ibm_bg_type" = "xP" ; then
        BLAS_CPPFLAGS="-I/opt/ibmmath/essl/4.4/include"
      fi
    fi

    if test "x$BLAS_LDFLAGS" = "x" ; then
      if test "x$cs_ibm_bg_type" = "xL" ; then
        BLAS_LDFLAGS="-L/opt/ibmmath/essl/4.2/lib"
      elif test "x$cs_ibm_bg_type" = "xP" ; then
        BLAS_LDFLAGS="-L/opt/ibmmath/essl/4.4/lib"
      fi
    fi

    AC_LANG_PUSH([Fortran])
    
    if test "$1" = "yes" -o "x$with_blas_libs" = "x"; then # Threaded version ?

      if test "x$cs_ibm_bg_type" = "xP" ; then
        BLAS_LIBS="-lesslsmpbg -lesslbg"
      else
        BLAS_LIBS="-lesslsmp"
      fi

      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

      AC_MSG_CHECKING([for smp ESSL BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                     [[      call ddot(0, 0, 0, 0, 0) ]])],
                     [ AC_DEFINE([HAVE_ESSL], 1, [ESSL BLAS support])
                       cs_have_blas=yes; with_blas_type=ESSL ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

    if test "$cs_have_blas" = "no" ; then # Test for non-threaded version
                                          # or explicitely specified libs second

      if test "x$cs_ibm_bg_type" != "x" ; then
        BLAS_LIBS="-lesslbg"
      elif test "x$with_blas_libs" != "x" -a "x$with_blas_type" = "xESSL"; then
        BLAS_LIBS="$with_blas_libs"
      else
        BLAS_LIBS="-lessl"
      fi

      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

      AC_MSG_CHECKING([for ESSL BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                     [[      call ddot(0, 0, 0, 0, 0) ]])],
                     [ AC_DEFINE([HAVE_ESSL], 1, [ESSL BLAS support])
                       cs_have_blas=yes; with_blas_type=ESSL ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

    AC_LANG_POP([Fortran])
    
    # Now check for header

    if test "$cs_have_blas" = "yes" ; then

      CPPFLAGS="${CPPFLAGS} ${BLAS_CPPFLAGS}"

      AC_MSG_CHECKING([for ESSL BLAS headers])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <essl.h>]],
                        [[ ddot(0, 0, 0, 0, 0); ]])],
                        [ AC_DEFINE([HAVE_ESSL_H], 1, [ESSL BLAS headers])
                          cs_have_essl_h=yes ],
                        [cs_have_essl_h=no])
      AC_MSG_RESULT($cs_have_essl_h)

   fi

  fi

  # Test for Intel MKL BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xMKL" ; then

    if test "x$with_blas_lib" = "x" ; then
      if test `uname -m` = ia64 ; then
        mkl_sub_lib="/64"
      elif test `uname -m` = x86_64 ; then
        mkl_sub_lib="/em64t"
      else
        mkl_sub_lib="/32"
      fi
    fi

    if test "$1" = "yes" -o "x$with_blas_libs" = "x"; then # Threaded version ?

      if test "`uname -m`" = "ia64" ; then
        BLAS_LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread"
      elif test `uname -m` = x86_64 ; then
        BLAS_LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread"
      else
        BLAS_LIBS="-lmkl_intel -lmkl_intel_thread -lmkl_core -lguide -lpthread"
      fi

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}${mkl_sub_lib}"
      LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

      AC_MSG_CHECKING([for threaded MKL BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mkl_cblas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_MKL], 1, [MKL BLAS support])
                       cs_have_blas=yes; with_blas_type=MKL ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

    if test "$cs_have_blas" = "no" ; then # Test for non-threaded version
                                          # or explicitely specified libs second

      if test "x$with_blas_libs" != "x" -a "x$with_blas_type" = "xMKL"; then
        BLAS_LIBS="$with_blas_libs"
      else
        if test "`uname -m`" = "ia64" ; then
          BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
        elif test `uname -m` = x86_64 ; then
          BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
        else
          BLAS_LIBS="-lmkl_intel -lmkl_sequential -lmkl_core"
        fi
      fi

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}${mkl_sub_lib}"
      LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

      AC_MSG_CHECKING([for MKL BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mkl_blas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_MKL], 1, [MKL BLAS support])
                       cs_have_blas=yes; with_blas_type=MKL ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

    if test "x$with_blas_type" = "xMKL" ; then
      BLAS_LDFLAGS="${BLAS_LDFLAGS}${mkl_sub_lib}"
    fi
    unset mkl_sub_lib

  fi

  # Test for ATLAS BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xATLAS" ; then

    if test "$1" = "yes" -o "x$with_blas_libs" = "x"; then # Threaded version ?

      BLAS_LIBS="-lptcblas -latlas -lpthread"

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

      AC_MSG_CHECKING([for threaded ATLAS BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_CBLAS], 1, [C BLAS support])
                       cs_have_blas=yes; with_blas_type=ATLAS ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

    if test "$cs_have_blas" = "no" ; then # Test for non-threaded version
                                          # or explicitely specified libs second
      if test "x$with_blas_libs" != "x" -a "x$with_blas_type" = "xATLAS"; then
        BLAS_LIBS="$with_blas_libs"
      else
        BLAS_LIBS="-lcblas -latlas"
      fi

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

      AC_MSG_CHECKING([for ATLAS BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_CBLAS], 1, [C BLAS support])
                       cs_have_blas=yes; with_blas_type=ATLAS ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

  fi

  # Test for generic C BLAS

  if test "x$with_blas_type" = "x" ; then

    if test "x$with_blas_libs" != "x" ; then
      BLAS_LIBS="$with_blas_libs"
    else
      BLAS_LIBS="-lblas"
    fi

    CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
    LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

    AC_MSG_CHECKING([for legacy C BLAS])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                   [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                   [ AC_DEFINE([HAVE_CBLAS], 1, [C BLAS support])
                     cs_have_blas=yes; with_blas_type=BLAS ],
                   [cs_have_blas=no])
    AC_MSG_RESULT($cs_have_blas)
  fi

  # Test for generic Fortran BLAS

  if test "x$with_blas_type" = "x" ; then

    AC_LANG_PUSH([Fortran])
    
    if test "x$with_blas_libs" != "x" ; then
      BLAS_LIBS="$with_blas_libs"
    else
      BLAS_LIBS="-lblas"
    fi

    LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
    LIBS=" ${saved_LIBS} ${BLAS_LIBS}"

    AC_MSG_CHECKING([for legacy Fortran BLAS])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                   [[      call ddot(0, 0, 0, 0, 0) ]])],
                   [ AC_DEFINE([HAVE_FBLAS], 1, [Fortran BLAS support])
                     cs_have_blas=yes; with_blas_type=BLAS ],
                   [cs_have_blas=no])
    AC_MSG_RESULT($cs_have_blas)

    AC_LANG_POP([Fortran])
    
  fi

  # Cleanup if no BLAS found

  if test "x$cs_have_blas" != "xyes"; then
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

AC_SUBST(BLAS_CPPFLAGS)
AC_SUBST(BLAS_LDFLAGS)
AC_SUBST(BLAS_LIBS)

])dnl


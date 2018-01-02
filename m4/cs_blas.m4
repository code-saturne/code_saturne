dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2018 EDF S.A.
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

# CS_AC_TEST_BLAS([use_threads])
#----------------
# modifies or sets cs_have_blas, BLAS_CPPFLAGS, BLAS_LDFLAGS, and BLAS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_BLAS], [

cs_have_blas=no

AC_ARG_WITH(blas,
            [AS_HELP_STRING([--with-blas=PATH],
                            [specify prefix directory for BLAS])],
            [if test "x$withval" = "x"; then
               with_blas=yes
             fi],
            [with_blas=no])

AC_ARG_WITH(blas-include,
            [AS_HELP_STRING([--with-blas-include=PATH],
                            [specify directory for BLAS include files])],
            [if test "x$with_blas" = "xcheck"; then
               with_blas=yes
             fi
             BLAS_CPPFLAGS="-I$with_blas_include"],
            [if test "x$with_blas" != "xno" -a "x$with_blas" != "xyes" \
	          -a "x$with_blas" != "xcheck"; then
               BLAS_CPPFLAGS="-I$with_blas/include"
             fi])

AC_ARG_WITH(blas-lib,
            [AS_HELP_STRING([--with-blas-lib=PATH],
                            [specify directory for BLAS library])],
            [if test "x$with_blas" = "xcheck"; then
               with_blas=yes
             fi
             BLAS_LDFLAGS="-L$with_blas_lib"
             # Add the libdir to the runpath as BLAS may not be libtoolized
             BLASRUNPATH="-R$with_blas_lib"],
            [if test "x$with_blas" != "xno" -a "x$with_blas" != "xyes" \
	          -a "x$with_blas" != "xcheck"; then
               BLAS_LDFLAGS="-L$with_blas/lib"
               # Add the libdir to the runpath as BLAS may not be libtoolized
               BLASRUNPATH="-R$with_blas/lib"
             fi])

AC_ARG_WITH(blas-type,
            [AS_HELP_STRING([--with-blas-type=NAME],
                            [force ATLAS, ESSL, MKL, ...])])

AC_ARG_WITH(blas-libs,
            [AS_HELP_STRING([--with-blas-libs=LIBS],
                            [specify BLAS libraries])])


if test "x$with_blas" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # Add known paths and libraries for Blue Gene/Q if not given

  # Test for IBM ESSL BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xESSL" ; then

    # Test compilation/header separately from link, as linking may require
    # Fortran libraries, and header is for C. Test library (link) first,
    # as header is only useful if library is present.

    # First, check for header

    CPPFLAGS="${CPPFLAGS} ${BLAS_CPPFLAGS}"

    AC_MSG_CHECKING([for ESSL BLAS headers])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <essl.h>]],
                      [[ ddot(0, 0, 0, 0, 0); ]])],
                      [cs_have_essl_h=yes],
                      [cs_have_essl_h=no])
    AC_MSG_RESULT($cs_have_essl_h)

    if test "$cs_have_essl_h" = "yes"; then

      AC_LANG_PUSH([Fortran])

      if test "$1" = "yes" -o "x$with_blas_libs" = "x"; then # Threaded version ?

        if test "x$cs_ibm_bg_type" = "xQ"; then
          BLAS_LIBS="-lesslbg"
          BLAS_LDFLAGS="-L$with_blas/lib64"
          BLASRUNPATH="-R$with_blas/lib64"
        else
          BLAS_LIBS="-lesslsmp"
        fi

        LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
        LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

        AC_MSG_CHECKING([for smp ESSL BLAS])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                       [[      call ddot(0, 0, 0, 0, 0) ]])],
                       [ AC_DEFINE([HAVE_ESSL], 1, [ESSL BLAS support])
                         cs_have_blas=yes; with_blas_type=ESSL ],
                       [cs_have_blas=no])
        AC_MSG_RESULT($cs_have_blas)
      fi

      if test "$cs_have_blas" = "no" ; then # Test for non-threaded version
                                            # or explicitely specified libs

        if test "x$cs_ibm_bg_type" = "xQ" ; then
          BLAS_LIBS="-lesslbg"
        elif test "x$with_blas_libs" != "x" -a "x$with_blas_type" = "xESSL"; then
          BLAS_LIBS="$with_blas_libs"
        else
          BLAS_LIBS="-lessl"
        fi

        LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
        LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

        AC_MSG_CHECKING([for ESSL BLAS])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                       [[      call ddot(0, 0, 0, 0, 0) ]])],
                       [ AC_DEFINE([HAVE_ESSL], 1, [ESSL BLAS support])
                         cs_have_blas=yes; with_blas_type=ESSL ],
                       [cs_have_blas=no])
        AC_MSG_RESULT($cs_have_blas)
      fi

      AC_LANG_POP([Fortran])

    fi

  fi

  # Test for Intel MKL libraries

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xMKL" ; then

    if test "x$with_blas_lib" = "x" ; then
      mkl_lib="$with_blas/lib"
      if test `uname -m` = ia64 ; then
        if test -d ${mkl_lib}/intel64 ; then
          mkl_sub_lib="/intel64"
        elif test -d ${mkl_lib}/64 ; then
          mkl_sub_lib="/64"
        fi
      elif test `uname -m` = x86_64 ; then
        if test -d ${mkl_lib}/intel64 ; then
          mkl_sub_lib="/intel64"
        elif test -d ${mkl_lib}/em64t ; then
          mkl_sub_lib="/em64t"
        fi
      else
        if test -d ${mkl_lib}/ia32 ; then
          mkl_sub_lib="/ia32"
        elif test -d ${mkl_lib}/32 ; then
          mkl_sub_lib="/32"
        fi
      fi
    fi

    if test "x$with_blas_libs" = "x"; then

      if test "x$with_blas_libs" != "x" -a "x$with_blas_type" = "xMKL"; then
        BLAS_LIBS="$with_blas_libs"
      elif  test "x$enable_shared" = xyes ; then
        if test "$1" = "yes" ; then # Threaded version ?
          case `uname -m` in
            *64)
              BLAS_LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm"
              ;;
            *)
              BLAS_LIBS="-lmkl_intel -lmkl_core -lmkl_intel_thread -lpthread -lm"
              ;;
          esac
        else
          case `uname -m` in
            *64)
              BLAS_LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm"
              ;;
            *)
              BLAS_LIBS="-lmkl_intel -lmkl_core -lmkl_sequential -lpthread -lm"
              ;;
          esac
        fi
      else
        if test "$1" = "yes" ; then # Threaded version ?
          case `uname -m` in
            *64)
              BLAS_LIBS="-Wl,--start-group ${mkl_lib}${mkl_sub_lib}/libmkl_intel_lp64.a ${mkl_lib}${mkl_sub_lib}/libmkl_core.a ${mkl_lib}${mkl_sub_lib}/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm"
              ;;
            *)
              BLAS_LIBS="-Wl,--start-group ${mkl_lib}${mkl_sub_lib}/libmkl_intel.a ${mkl_lib}${mkl_sub_lib}/libmkl_core.a ${mkl_lib}${mkl_sub_lib}/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm"
              ;;
          esac
        else
          case `uname -m` in
            *64)
              BLAS_LIBS="-Wl,--start-group ${mkl_lib}${mkl_sub_lib}/libmkl_intel_lp64.a ${mkl_lib}${mkl_sub_lib}/libmkl_core.a ${mkl_lib}${mkl_sub_lib}/libmkl_sequential.a -Wl,--end-group -lpthread -lm"
              ;;
            *)
              BLAS_LIBS="-Wl,--start-group ${mkl_lib}${mkl_sub_lib}/libmkl_intel.a ${mkl_lib}${mkl_sub_lib}/libmkl_core.a ${mkl_lib}${mkl_sub_lib}/libmkl_sequential.a -Wl,--end-group -lpthread -lm"
              ;;
          esac
        fi
      fi

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}${mkl_sub_lib}"
      LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

      AC_MSG_CHECKING([for MKL libraries])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mkl_cblas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_MKL], 1, [MKL libraries support])
                       cs_have_blas=yes; with_blas_type=MKL ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)

    fi

    if test "x$with_blas_type" = "xMKL" ; then
      BLAS_LDFLAGS="${BLAS_LDFLAGS}${mkl_sub_lib}"
      BLASRUNPATH="${BLASRUNPATH}${mkl_sub_lib}"
    fi
    unset mkl_sub_lib

  fi

  # Test for ATLAS BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xATLAS" ; then

    if test "$1" = "yes" -o "x$with_blas_libs" = "x"; then # Threaded version ?

      BLAS_LIBS="-lptcblas -latlas -lpthread"

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

      AC_MSG_CHECKING([for threaded ATLAS BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_ATLAS], 1, [ATLAS BLAS support])
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
      LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

      AC_MSG_CHECKING([for ATLAS BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                     [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_ATLAS], 1, [ATLAS BLAS support])
                       cs_have_blas=yes; with_blas_type=ATLAS ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

  fi

  # Test for AMD ACML BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xACML" ; then

    if test "$1" = "yes" -o "x$with_blas_libs" = "x"; then # Threaded version ?

      BLAS_LIBS="-lacml_mp"

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

      AC_MSG_CHECKING([for threaded ACML BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <acml.h>]],
                     [[ ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_ACML], 1, [ACML BLAS support])
                       cs_have_blas=yes; with_blas_type=ACML ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

    if test "$cs_have_blas" = "no" ; then # Test for non-threaded version
                                          # or explicitely specified libs second
      if test "x$with_blas_libs" != "x" -a "x$with_blas_type" = "xACML"; then
        BLAS_LIBS="$with_blas_libs"
      else
        BLAS_LIBS="-lacml"
      fi

      CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
      LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
      LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

      AC_MSG_CHECKING([for ACML BLAS])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <acml.h>]],
                     [[ ddot(0, 0, 0, 0, 0); ]])],
                     [ AC_DEFINE([HAVE_ACML], 1, [ACML BLAS support])
                       cs_have_blas=yes; with_blas_type=ACML ],
                     [cs_have_blas=no])
      AC_MSG_RESULT($cs_have_blas)
    fi

  fi

  # Test for generic C BLAS

  if test "x$with_blas_type" = "x" -o "x$with_blas_type" = "xBLAS" ; then

    if test "x$with_blas_libs" != "x" ; then
      BLAS_LIBS="$with_blas_libs"
    else
      BLAS_LIBS="-lblas"
    fi

    CPPFLAGS="${saved_CPPFLAGS} ${BLAS_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${BLAS_LDFLAGS}"
    LIBS=" ${BLAS_LIBS} ${saved_LIBS}"

    AC_MSG_CHECKING([for legacy C BLAS])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <cblas.h>]],
                   [[ cblas_ddot(0, 0, 0, 0, 0); ]])],
                   [ AC_DEFINE([HAVE_CBLAS], 1, [C BLAS support])
                     cs_tmp_have_cblas=yes; with_blas_type=BLAS ],
                   [cs_tmp_have_cblas=no])
    AC_MSG_RESULT($cs_tmp_have_cblas)

    if test "x$cs_tmp_have_cblas" = "xyes"; then
      cs_have_blas=yes
    fi

  fi

  # Cleanup if no BLAS found

  if test "x$cs_have_blas" != "xyes"; then
    if test "x$with_blas" != "xcheck" ; then
      AC_MSG_FAILURE([BLAS support is requested, but test for BLAS failed!])
    else
      AC_MSG_WARN([no BLAS support])
    fi
    BLAS_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

# ESSL requires Fortran to link, and MKL provides sparse matrix-vector
# operations (so it may be used by the Code_Saturne solver)
AM_CONDITIONAL(HAVE_ESSL, test x$with_blas_type = xESSL)
AM_CONDITIONAL(HAVE_MKL, test x$with_blas_type = xMKL)

AC_SUBST(cs_have_blas)
AC_SUBST(BLAS_CPPFLAGS)
AC_SUBST(BLAS_LDFLAGS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(BLASRUNPATH)

])dnl


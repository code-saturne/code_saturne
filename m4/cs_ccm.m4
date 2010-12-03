dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2010 EDF S.A., France
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

# CS_AC_TEST_CCM
#---------------
# modifies or sets cs_have_ccm, CCM_CPPFLAGS, CCM_LDFLAGS, and CCM_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CCM], [

cs_have_ccm=no
cs_have_ccm_headers=no

AC_ARG_WITH(ccm,
            [AS_HELP_STRING([--with-ccm=DIR],
                            [specify prefix directory for CCMIO])],
            [if test "x$withval" = "x"; then
               with_ccm=yes
             fi],
            [with_ccm=check])

AC_ARG_WITH(ccm-include,
            [AS_HELP_STRING([--with-ccm-include=DIR],
                            [specify directory for CCMIO include files])],
            [if test "x$with_ccm" = "xcheck"; then
               with_ccm=yes
             fi
             CCM_CPPFLAGS="-I$with_ccm_include"],
            [if test "x$with_ccm" != "xno" -a "x$with_ccm" != "xyes" \
	          -a "x$with_ccm" != "xcheck"; then
               CCM_CPPFLAGS="-I$with_ccm/include"
             fi])

AC_ARG_WITH(ccm-lib,
            [AS_HELP_STRING([--with-ccm-lib=DIR],
                            [specify directory for CCMIO library])],
            [if test "x$with_ccm" = "xcheck"; then
               with_ccm=yes
             fi
             CCM_LDFLAGS="-L$with_ccm_lib"
             # Add the libdir to the runpath as CCM is not libtoolized
             CCMRUNPATH="-R$with_ccm_lib"],
            [if test "x$with_ccm" != "xno" -a "x$with_ccm" != "xyes" \
	          -a "x$with_ccm" != "xcheck"; then
               CCM_LDFLAGS="-L$with_ccm/lib"
               # Add the libdir to the runpath as CCM is not libtoolized
               CCMRUNPATH="-R$with_ccm/lib"
             fi])

# ADF may be provided directly (patched ADF with libccmio)
# or through CGNS

if test "x$with_ccm" != "xno" -a "x$cs_have_adf" = "xno" -a "x$cs_have_cgns" = "xno"
then
  if test "x$with_ccm" = "xcheck"; then
    with_ccm=no
    AC_MSG_WARN([no ADF library found; will not search for CCM])
  else
    AC_MSG_ERROR([no ADF library found; required for CCM])
  fi
fi

if test "x$with_ccm" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$ADF_LIBS" != "x" ; then
    CCM_LIBS="-lccmio $ADF_LIBS"
    CPPFLAGS="${CPPFLAGS} ${CCM_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${CCM_LDFLAGS} $ADF_LDFLAGS"
  elif test "x$CGNS_LIBS" != "x" ; then
    CCM_LIBS="-lccmio"
    CPPFLAGS="${CPPFLAGS} ${CCM_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${CCM_LDFLAGS} $CGNS_LDFLAGS $HDF5_LDFLAGS"
  fi
  LIBS="${LIBS} ${CCM_LIBS} $CGNS_LIBS $HDF5_LIBS"

# Check that CCMIO header files exist

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <libccmio/ccmio.h>]],
[[int i = kCCMIONoErr;]])],
                    [AC_MSG_RESULT([CCMIO headers found])
                     cs_have_ccm_headers=yes
                    ],
                    [AC_MSG_RESULT([CCMIO headers not found])
                    ])

  if test "x$cs_have_ccm_headers" = "xyes"; then
    AC_CHECK_LIB(ccmio, CCMIOOpenFile, 
                 [ AC_DEFINE([HAVE_CCM], 1, [CCM file support])
                   cs_have_ccm=yes
                 ], 
                 [if test "x$with_ccm" != "xcheck" ; then
                    AC_MSG_FAILURE([CCM support is requested, but test for CCM failed!])
                  else
                    AC_MSG_WARN([no CCM file support])
                  fi
                 ],
                 )
  fi

  if test "x$cs_have_ccm" != "xyes"; then
    CCM_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_CCM, test x$cs_have_ccm = xyes)

AC_SUBST(cs_have_ccm)
AC_SUBST(CCM_CPPFLAGS)
AC_SUBST(CCM_LDFLAGS)
AC_SUBST(CCM_LIBS)
AC_SUBST(CCMRUNPATH)

])dnl


dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2022 EDF S.A.
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

if test "x$with_ccm" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # ADF is be provided directly (patched ADF with libccmio)
  # We must be careful not to use CGNS's adf, as this leads to nonworking
  # CCM-IO builds. CCM's LDFLAGS must thus come first...

  CCM_LIBS="-lccmio -ladf"
  CPPFLAGS="${CPPFLAGS} ${CCM_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${CCM_LDFLAGS}"
  LIBS="${CCM_LIBS} ${LIBS}"

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

    AC_MSG_CHECKING([for CCM file support])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <libccmio/ccmio.h>]],
[[CCMIOID root;
CCMIOError error = kCCMIONoErr;
CCMIOOpenFile(&error, "test.ccm", kCCMIOWrite, &root);]])
                   ],
                   [ AC_DEFINE([HAVE_CCM], 1, [CCM file support])
                     cs_have_ccm=yes
                   ],
                   [],
                   )
    AC_MSG_RESULT($cs_have_ccm)
    if test "x$cs_have_ccm" = "xno" ; then
      if test "x$with_ccm" != "xcheck" ; then
        AC_MSG_FAILURE([CCM support is requested, but test for CCM failed!])
      fi
    fi
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


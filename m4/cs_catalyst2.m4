dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
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

# CS_AC_TEST_CATALYST2
#--------------------
# modifies or sets cs_have_catalyst2, CATALYST2_CPPFLAGS, CATALYST2_LDFLAGS,
# and CATALYST2_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_CATALYST2], [

cs_have_catalyst2=no
cs_have_plugin_catalyst2=yes

# Configure options
#------------------

AC_ARG_WITH(catalyst2,
            [AS_HELP_STRING([--with-catalyst2=PATH],
                            [specify prefix directory for Catalyst2])],
            [if test "x$withval" = "x"; then
               with_catalyst2=no
             fi],
            [with_catalyst2=no])

AC_ARG_WITH(catalyst2-include,
            [AS_HELP_STRING([--with-catalyst2-include=PATH],
                            [specify directory for Catalyst2 include files])],
            [if test "x$with_catalyst2" = "xcheck"; then
               with_catalyst2=yes
             fi
             CATALYST2_CPPFLAGS="-I$with_catalyst2_include"],
            [if test "x$with_catalyst2" != "xno" -a "x$with_catalyst2" != "xyes" \
	          -a "x$with_catalyst2" != "xcheck"; then
               CATALYST2_CPPFLAGS="-I$with_catalyst2/include/catalyst-2.0"
             fi])

AC_ARG_WITH(catalyst2-lib,
            [AS_HELP_STRING([--with-catalyst2-lib=PATH],
                            [specify directory for Catalyst2 library])],
            [if test "x$with_catalyst2" = "xcheck"; then
               with_catalyst2=yes
             fi
             CATALYST2_LDFLAGS="-L$with_catalyst2_lib"
             CATALYST2RUNPATH="-R$with_catalyst2_lib"],
            [if test "x$with_catalyst2" != "xno" -a "x$with_catalyst2" != "xyes" \
	          -a "x$with_catalyst2" != "xcheck"; then
               CATALYST2_LDFLAGS="-L$with_catalyst2/lib"
               CATALYST2RUNPATH="-R$with_catalyst2/lib"
             fi])

AC_ARG_ENABLE(catalyst2-as-plugin,
  [AS_HELP_STRING([--disable-catalyst2-as-plugin], [do not use Catalyst2 as plugin])],
  [
    case "${enableval}" in
      yes) cs_have_plugin_catalyst2=yes ;;
      no)  cs_have_plugin_catalyst2=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-catalyst2-as-plugin]) ;;
    esac
  ],
  [ cs_have_plugin_catalyst2=yes ]
)

if test x$cs_have_dlloader = xno -o x$enable_shared = xno ; then
  cs_have_plugin_catalyst2=no
fi

# Check for a Catalyst2 library
#------------------------------

if test "x$with_catalyst2" != "xno" ; then

  AC_MSG_CHECKING([ParaView/Catalyst2])

  CATALYST2_LIBS=-lcatalyst

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${MPI_CPPFLAGS} ${CATALYST2_CPPFLAGS} ${CPPFLAGS}"
  LDFLAGS="${saved_LDFLAGS} ${CATALYST2_LDFLAGS}"
  LIBS=" ${CATALYST2_LIBS} ${saved_LIBS}"

  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <catalyst.hpp>
#include <catalyst_version.h>]],
[[auto node = conduit_node_create();
  catalyst_initialize(node);]])
                    ],
                    [ cs_have_catalyst2=yes ],
                    [ AC_MSG_WARN([no Catalyst2 support]) ],
                  )
  AC_LANG_POP([C++])

  if test "x$cs_have_catalyst2" = "xno"; then
    CATALYST2_CPPFLAGS=""
    CATALYST2_LDFLAGS=""
    CATALYST2_LIBS=""
  else
    CPPFLAGS="${saved_CPPFLAGS} ${CATALYST2_CPPFLAGS}"
  fi

  if test "x$CATALYST2_LIBS" != "x" ; then
    cs_have_catalyst2=yes;
  fi

  AC_MSG_RESULT($cs_have_catalyst2)

  # Report Catalyst2 support
  #-------------------------

  if test "x$cs_have_catalyst2" = "xyes" ; then
    AC_DEFINE([HAVE_CATALYST2], 1, [Catalyst2 co-processing support])
    AC_DEFINE([HAVE_CONDUIT], 1, [Conduit library support])
    if test x$cs_have_plugin_catalyst2 = xyes ; then
      AC_DEFINE([HAVE_PLUGIN_CATALYST2], 1, [Catalyst2 co-processing support as plugin])
    fi
  elif test "x$cs_have_catalyst2" = "xno" ; then
    if test "x$with_catalyst2" != "xcheck" ; then
      AC_MSG_FAILURE([Catalyst2 co-processing support requested, but test for Catalyst2 failed!])
    else
      AC_MSG_WARN([no Catalyst2 co-processing support])
    fi
  else
    AC_MSG_WARN([no Catalyst2 test as expected])
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

if test x$cs_have_catalyst2 = xno ; then
  cs_have_plugin_catalyst2=no
fi

AM_CONDITIONAL(HAVE_CATALYST2, test x$cs_have_catalyst2 = xyes)
AM_CONDITIONAL(HAVE_PLUGIN_CATALYST2, test x$cs_have_plugin_catalyst2 = xyes)

AC_SUBST(cs_have_catalyst2)
AC_SUBST(cs_have_plugin_catalyst2)
AC_SUBST(CATALYST2_CPPFLAGS)
AC_SUBST(CATALYST2_CXXFLAGS)
AC_SUBST(CATALYST2_LDFLAGS)
AC_SUBST(CATALYST2_LIBS)
AC_SUBST(CATALYST2RUNPATH)

])dnl

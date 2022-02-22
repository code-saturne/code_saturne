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

# CS_AC_TEST_SCOTCH
#-----------------
# modifies or sets cs_have_scotch, SCOTCH_CPPFLAGS, SCOTCH_LDFLAGS, and SCOTCH_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SCOTCH], [

cs_have_ptscotch_header=no
cs_have_ptscotch=no
cs_have_scotch_header=no
cs_have_scotch=no
cs_scotch_ge_6=no
scotch_prefix=""

# Some recent Linux distributions require explicit exports
# of all symbols in libraries
cs_scotch_test_ladd=''
$CC -Xlinker --help | grep "no-as-needed" > /dev/null 2>&1
if test "$?" = "0" ; then
  cs_scotch_test_ladd='-Wl,--no-as-needed '
fi

# Common library dependencies for PT-SCOTCH
cs_scotch_l0="-lm"
cs_scotch_l1="-lz -lm"
cs_scotch_l2="-lm -lpthread"
cs_scotch_l3="-lz -lm -lpthread"
cs_scotch_l4="-lm -lpthread -lrt"
cs_scotch_l5="-lz -lm -lpthread -lrt"
SCOTCH_LIBS_ADD=""

AC_ARG_WITH(scotch,
            [AS_HELP_STRING([--with-scotch=PATH],
                            [specify prefix directory for SCOTCH])],
            [if test "x$withval" = "x"; then
               with_scotch=yes
             elif test "x$withval" = "xsalome"; then
               if test -z "$SCOTCHDIR"; then
                 AC_MSG_FAILURE([no SALOME path information for SCOTCH (needed by --with-scotch=salome)!])
               else
                 with_scotch=$SCOTCHDIR
               fi
             fi],
            [with_scotch=no])

AC_ARG_WITH(scotch-include,
            [AS_HELP_STRING([--with-scotch-include=PATH],
                            [specify directory for SCOTCH include files])],
            [if test "x$with_scotch" = "xcheck" -o "x$with_scotch" = "xno"; then
               with_scotch=yes
             fi
             SCOTCH_CPPFLAGS="-I$with_scotch_include"],
            [if test "x$with_scotch" != "xno" ; then
               if test "x$with_scotch" != "xyes" \
	               -a "x$with_scotch" != "xcheck"; then
                 SCOTCH_CPPFLAGS="-I$with_scotch/include"
               fi
             fi])

AC_ARG_WITH(scotch-lib,
            [AS_HELP_STRING([--with-scotch-lib=PATH],
                            [specify directory for SCOTCH library])],
            [if test "x$with_scotch" = "xcheck" -o "x$with_scotch" = "xno"; then
               with_scotch=yes
             fi
             SCOTCH_LDFLAGS="-L$with_scotch_lib"
             # Add the libdir to the runpath as SCOTCH is not libtoolized
             SCOTCHRUNPATH="-R$with_scotch_lib"],
            [if test "x$with_scotch" != "xno" -a "x$with_scotch" != "xyes" \
	          -a "x$with_scotch" != "xcheck"; then
               SCOTCH_LDFLAGS="-L$with_scotch/lib"
               # Add the libdir to the runpath as SCOTCH is not libtoolized
               SCOTCHRUNPATH="-R$with_scotch/lib"
             fi])


if test "x$with_scotch" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # Test for PT-SCOTCH first

  # Check for ptscotch.h header
  CPPFLAGS="$saved_CPPFLAGS $SCOTCH_CPPFLAGS $MPI_CPPFLAGS"
  AC_CHECK_HEADERS([ptscotch.h],
                   [cs_have_ptscotch_header=yes],
                   [],
                   [#include <stdio.h>
                    #include <stdint.h>
                    #include <mpi.h>])

  # Second test if scotch path not specified
  if test "x$cs_have_ptscotch_header" = "xno" -a "$SCOTCH_CPPFLAGS" = "-I/usr/include" ; then
    unset ac_cv_header_ptscotch_h
    CPPFLAGS="$saved_CPPFLAGS -I/usr/include/scotch $MPI_CPPFLAGS"
    AC_CHECK_HEADERS([ptscotch.h],
                     [cs_have_ptscotch_header=yes
                      SCOTCH_CPPFLAGS=-I/usr/include/scotch],
                     [],
                     [#include <stdio.h>
                      #include <stdint.h>
                      #include <mpi.h>])
  fi

  LDFLAGS="${LDFLAGS} ${SCOTCH_LDFLAGS} ${MPI_LDFLAGS}"
  SCOTCH_LIBS="-lptscotch -lptscotcherr"

  if test "x$cs_have_ptscotch_header" = "xyes" ; then

    AC_MSG_CHECKING([for PT-SCOTCH])

    # Check if SCOTCH version is 6 or 5, as libptscotch version 5.1.x includes libscotch,
    # while version 6.0.x requires it.

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <ptscotch.h>]],
[[#if SCOTCH_VERSION < 6
# error test for SCOTCH version 6 so assume 5.1
#endif
]])],
                      [cs_scotch_ge_6=yes
                       SCOTCH_LIBS="-lptscotch -lptscotcherr -lscotch -lscotcherr"],
                      [])

    for cs_scotch_ladd in "$cs_scotch_l0" "$cs_scotch_l1" "$cs_scotch_l2" "$cs_scotch_l3" "$cs_scotch_l4" "$cs_scotch_l5"
    do
      if test "x$cs_have_ptscotch" = "xno" ; then
        LIBS="${cs_scotch_test_ladd}${SCOTCH_LIBS} ${cs_scotch_ladd} ${MPI_LIBS} ${saved_LIBS}"
        AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <stdint.h>
#include <mpi.h>
#include <ptscotch.h>]],
[[ SCOTCH_dgraphInit((void *)0, MPI_COMM_WORLD); ]])],
[cs_have_ptscotch=yes
 cs_have_scotch=yes
 SCOTCH_LIBS_ADD="${cs_scotch_ladd}"],
[cs_have_ptscotch=no])
      fi
    done

    AC_MSG_RESULT($cs_have_ptscotch)

  fi

  # Test for SCOTCH second

  if test "x$cs_have_ptscotch" = "xno"; then

    # Check for scotch.h header
    CPPFLAGS="$saved_CPPFLAGS $SCOTCH_CPPFLAGS"
    AC_CHECK_HEADERS([scotch.h],
                       [cs_have_scotch_header=yes],
                       [],
                       [])

    if test "x$cs_have_scotch_header" = "xno" ; then
      unset ac_cv_header_scotch_h
      SCOTCH_CPPFLAGS="$saved_CPPFLAGS -I/usr/include/scotch"
      CPPFLAGS="$saved_CPPFLAGS $SCOTCH_CPPFLAGS"
      AC_CHECK_HEADERS([scotch.h],
                       [cs_have_scotch_header=yes],
                       [],
                       [])
    fi

    LDFLAGS="${saved_LDFLAGS} ${SCOTCH_LDFLAGS}"
    SCOTCH_LIBS="-lscotch -lscotcherr"

    AC_MSG_CHECKING([for SCOTCH])

    for cs_scotch_ladd in "$cs_scotch_l0" "$cs_scotch_l1" "$cs_scotch_l2" "$cs_scotch_l3" "$cs_scotch_l4" "$cs_scotch_l5"
    do
      if test "x$cs_have_scotch" = "xno" ; then
        LIBS="${cs_scotch_test_ladd}${SCOTCH_LIBS} ${cs_scotch_ladd} ${saved_LIBS}"
        AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <stdio.h>
#include <stdint.h>
#include <scotch.h>]],
[[ SCOTCH_graphInit((void *)0); ]])],
[cs_have_scotch=yes
 SCOTCH_LIBS_ADD="${cs_scotch_ladd}"],
[cs_have_scotch=no])
      fi
    done

  fi

  # libptscotcherr / libscotcherr functions in cs_partition.c, so do not use these libraries

  if test "x$cs_have_ptscotch" = "xyes"; then
    AC_DEFINE([HAVE_PTSCOTCH], 1, [use SCOTCH])
    if test "x$cs_scotch_ge_6" = "xyes" ; then
      SCOTCH_LIBS="-lptscotch -lscotch ${SCOTCH_LIBS_ADD}"
    else
      SCOTCH_LIBS="-lptscotch ${SCOTCH_LIBS_ADD}"
    fi
  elif test "x$cs_have_scotch" = "xyes"; then
    AC_DEFINE([HAVE_SCOTCH], 1, [use SCOTCH])
    SCOTCH_LIBS="-lscotch ${SCOTCH_LIBS_ADD}"
  else
    SCOTCH_CPPFLAGS=""
    SCOTCH_LDFLAGS=""
    SCOTCH_LIBS=""
  fi

  # Report PT-SCOTCH/SCOTCH support
  #------------------------

  if test "x$cs_have_ptscotch" = "xno" -a "x$cs_have_scotch" = "xno" ; then
    if test "x$with_scotch" != "xcheck" ; then
      AC_MSG_FAILURE([PT-SCOTCH/SCOTCH support is requested, but test for SCOTCH failed!])
    else
      AC_MSG_WARN([no PT-SCOTCH/SCOTCH partitioner support])
    fi
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  case $host_os in
    mingw64)
      scotch_prefix=`cygpath --path --windows "$with_scotch"`;;
    *)
      ;;
  esac
fi

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS
unset cs_have_ptscotch_header
unset cs_have_scotch_header
unset cs_scotch_ge_6
unset cs_scotch_l0
unset cs_scotch_l1
unset cs_scotch_l2
unset cs_scotch_l3
unset cs_scotch_l4
unset cs_scotch_l5

AC_SUBST(cs_have_scotch)
AC_SUBST(scotch_prefix, [${scotch_prefix}])
AC_SUBST(SCOTCH_CPPFLAGS)
AC_SUBST(SCOTCH_LDFLAGS)
AC_SUBST(SCOTCH_LIBS)
AC_SUBST(SCOTCHRUNPATH)

])dnl


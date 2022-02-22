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

# CS_AC_TEST_MELISSA
#--------------------
# modifies or sets cs_have_melissa, MELISSA_CPPFLAGS, MELISSA_LDFLAGS,
# and MELISSA_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_MELISSA], [

cs_have_melissa=no
cs_have_melissa_mpi=no
cs_have_melissa_mpi_05=no
cs_have_melissa_no_mpi=no
cs_have_plugin_melissa=yes

# ZeroMQ paths
#-------------

AC_ARG_WITH(zeromq,
            [AS_HELP_STRING([--with-zeromq=PATH],
                            [specify prefix directory for ZeroMQ])],
            [if test "x$withval" = "x"; then
               with_zeromq=yes
             fi],
            [with_zeromq=check])

AC_ARG_WITH(zeromq-lib,
            [AS_HELP_STRING([--with-zeromq-lib=DIR],
                            [specify directory for ZeroMQ library])],
            [if test "x$with_zeromq" = "xcheck"; then
               with_zeromq=yes
             fi
             ZEROMQ_LDFLAGS="-L$with_zeromq_lib"
             # Add the libdir to the runpath as ZeroMQ is not libtoolized
             ZEROMQRUNPATH="-R$with_zeromq_lib"],
            [if test "x$with_zeromq" != "xno" -a "x$with_zeromq" != "xyes" \
	          -a "x$with_zeromq" != "xcheck"; then
               ZEROMQ_LDFLAGS="-L$with_zeromq/lib"
               # Add the libdir to the runpath as zeromq is not libtoolized
               ZEROMQRUNPATH="-R$with_zeromq/lib"
             fi])

# Configure options for Melissa paths
#------------------------------------

AC_ARG_WITH(melissa,
            [AS_HELP_STRING([--with-melissa=PATH],
                            [specify prefix directory for MELISSA])],
            [if test "x$withval" = "x"; then
               with_melissa=no
             fi],
            [with_melissa=no])

AC_ARG_WITH(melissa-include,
            [AS_HELP_STRING([--with-melissa-include=DIR],
                            [specify directory for melissa include files])],
            [if test "x$with_melissa" = "xcheck"; then
               with_melissa=yes
             fi
             MELISSA_CPPFLAGS="-I$with_melissa_include"],
            [if test "x$with_melissa" != "xno" -a "x$with_melissa" != "xyes" \
	          -a "x$with_melissa" != "xcheck"; then
               MELISSA_CPPFLAGS="-I$with_melissa/include"
             fi])

AC_ARG_WITH(melissa-lib,
            [AS_HELP_STRING([--with-melissa-lib=DIR],
                            [specify directory for melissa library])],
            [if test "x$with_melissa" = "xcheck"; then
               with_melissa=yes
             fi
             MELISSA_LDFLAGS="-L$with_melissa_lib"
             # Add the libdir to the runpath as melissa is not libtoolized
             MELISSARUNPATH="-R$with_melissa_lib"],
            [if test "x$with_melissa" != "xno" -a "x$with_melissa" != "xyes" \
	          -a "x$with_melissa" != "xcheck"; then
               MELISSA_LDFLAGS="-L$with_melissa/lib"
               # Add the libdir to the runpath as melissa is not libtoolized
               MELISSARUNPATH="-R$with_melissa/lib"
             fi])

AC_ARG_ENABLE(melissa-as-plugin,
  [AS_HELP_STRING([--disable-melissa-as-plugin], [do not use Melissa as plugin])],
  [
    case "${enableval}" in
      yes) cs_have_plugin_melissa=yes ;;
      no)  cs_have_plugin_melissa=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-melissa-as-plugin]) ;;
    esac
  ],
  [ cs_have_plugin_melissa=yes ]
)

if test x$cs_have_dlloader = xno -o x$enable_shared = xno ; then
  cs_have_plugin_melissa=no
fi

# Now check for libraries
#------------------------

if test "x$with_zeromq" != "xno" -a "x$with_melissa" != "xno"; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  MELISSA_LDFLAGS="${MELISSA_LDFLAGS} ${ZEROMQ_LDFLAGS}"

  if test "x$cs_have_mpi" != "xno"; then

    AC_MSG_CHECKING([for Melissa library (with MPI)])

    MELISSA_LIBS="-lmelissa -lzmq"
    CPPFLAGS="${CPPFLAGS} ${MELISSA_CPPFLAGS} ${MPI_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${MELISSA_LDFLAGS} ${MPI_LDFLAGS}"
    LIBS="${MELISSA_LIBS} ${MPI_LIBS} ${saved_LIBS}"

    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <mpi.h>
#include <melissa/api.h>]],
[[(void)melissa_init("name", 1, MPI_COMM_SELF); ]])],
                   [cs_have_melissa_mpi=yes],
                   [cs_have_melissa_mpi=no])

    # Test for version 0.5 or older, if not found
    if test "x$cs_have_melissa_mpi" = "xno"; then

      MELISSA_LIBS="-lmelissa_api -lzmq"
      LIBS="${MELISSA_LIBS} ${MPI_LIBS} ${saved_LIBS}"

      AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <mpi.h>
#include <melissa_api.h>]],
[[(void)melissa_init("name", 1, MPI_COMM_SELF); ]])],
                   [cs_have_melissa_mpi_05=yes],
                   [cs_have_melissa_mpi_05=no])

    fi

  fi

  if test "x$cs_have_melissa_mpi" = "xyes"; then
    cs_have_melissa=yes
  elif test "x$cs_have_melissa_mpi_05" = "xyes"; then
    cs_have_melissa=yes
    cs_have_melissa_no_mpi=yes
  fi

  if test "x$cs_have_melissa_no_mpi" = "xno"; then

    AC_MSG_CHECKING([for Melissa library without MPI])

    CPPFLAGS="${saved_CPPFLAGS} ${MELISSA_CPPFLAGS}"
    LDFLAGS="${saved_LDFLAGS} ${MELISSA_LDFLAGS}"

    if test "x$cs_have_melissa_mpi" = "xno"; then
      MELISSA_LIBS="-lmelissa_api -lzmq"
    fi

    LIBS="${MELISSA_LIBS} ${_saved_LIBS}"
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <melissa_api_no_mpi.h>]],
[[(void)melissa_init_no_mpi("name", 1); ]])],
                   [cs_have_melissa_no_mpi=yes],
                   [cs_have_melissa_no_mpi=no])

    if test "x$cs_have_melissa_no_mpi" = "xyes"; then
      cs_have_melissa = yes
    fi

  fi

  # Report Melissa support
  #------------------------

  if test "x$cs_have_melissa" = "xyes" ; then
    AC_DEFINE([HAVE_MELISSA], 1, [Melissa co-processing support])
    if test "x$cs_have_melissa_mpi" = "xyes" ; then
      AC_DEFINE([HAVE_MELISSA_MPI], 1, [Melissa co-processing support with MPI])
    fi
    if test "x$cs_have_melissa_mpi_05" = "xyes" ; then
      AC_DEFINE([HAVE_MELISSA_MPI_05], 1, [Melissa 0.5 or older co-processing support with MPI])
    fi
    if test "x$cs_have_melissa_no_mpi" = "xyes" ; then
      AC_DEFINE([HAVE_MELISSA_NO_MPI], 1, [Melissa co-processing support without MPI])
    fi
    if test x$cs_have_plugin_melissa = xyes ; then
      AC_DEFINE([HAVE_PLUGIN_MELISSA], 1, [Melissa co-processing support as plugin])
    fi
  elif test "x$cs_have_melissa" = "xno" ; then
    if test "x$with_melissa" != "xcheck" ; then
      AC_MSG_FAILURE([Melissa co-processing support requested, but test for Melissa failed!])
    else
      AC_MSG_WARN([no Melissa co-processing support])
    fi
  fi

  AC_MSG_RESULT($cs_have_melissa)

  if test "x$cs_have_melissa" != "xyes"; then
    MELISSA_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MELISSA, test x$cs_have_melissa = xyes)

if test x$cs_have_melissa = xno ; then
  cs_have_plugin_melissa=no
fi

AM_CONDITIONAL(HAVE_MELISSA, test x$cs_have_melissa = xyes)
AM_CONDITIONAL(HAVE_MELISSA_LEGACY, test x$cs_have_melissa_05 = xyes)
AM_CONDITIONAL(HAVE_PLUGIN_MELISSA, test x$cs_have_plugin_melissa = xyes)

cs_py_have_plugin_melissa=False
if test x$cs_have_plugin_melissa = xyes ; then
  cs_py_have_plugin_melissa=True
fi

AC_SUBST(cs_have_melissa)
AC_SUBST(cs_py_have_plugin_melissa)
AC_SUBST(MELISSA_CPPFLAGS)
AC_SUBST(MELISSA_LDFLAGS)
AC_SUBST(MELISSA_LIBS)
AC_SUBST(MELISSARUNPATH)

])dnl

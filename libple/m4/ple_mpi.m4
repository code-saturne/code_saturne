dnl Copyright (C) 2005-2010 EDF
dnl
dnl This file is part of the PLE software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl PLE source distribution.

# PLE_AC_TEST_MPI
#----------------
# optional MPI support (use CC=mpicc with configure if necessary)
# modifies or sets ple_have_mpi, MPI_CPPFLAGS, MPI_LDFLAGS, and MPI_LIBS
# depending on libraries found

AC_DEFUN([PLE_AC_TEST_MPI], [

saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"

ple_have_mpi=no
ple_have_mpi_header=no
ple_have_mpi_one_sided=no

AC_ARG_WITH(mpi,
            [AS_HELP_STRING([--with-mpi=PATH],
                            [specify prefix directory for MPI])],
            [if test "x$withval" = "x"; then
               with_mpi=yes
             fi],
            [with_mpi=check])

AC_ARG_WITH(mpi-include,
            [AS_HELP_STRING([--with-mpi-include=PATH],
                            [specify directory for MPI include files])],
            [if test "x$with_mpi" = "xcheck"; then
               with_mpi=yes
             fi
             MPI_CPPFLAGS="-I$with_mpi_include"],
            [if test "x$with_mpi" != "xno" -a "x$with_mpi" != "xyes" \
	          -a "x$with_mpi" != "xcheck"; then
               MPI_CPPFLAGS="-I$with_mpi/include"
             fi])

AC_ARG_WITH(mpi-lib,
            [AS_HELP_STRING([--with-mpi-lib=PATH],
                            [specify directory for MPI library])],
            [if test "x$with_mpi" = "xcheck"; then
               with_mpi=yes
             fi
             MPI_LDFLAGS="-L$with_mpi_lib"],
            [if test "x$with_mpi" != "xno" -a "x$with_mpi" != "xyes" \
	          -a "x$with_mpi" != "xcheck"; then
               MPI_LDFLAGS="-L$with_mpi/lib"
             fi])


# Just in case, remove excess whitespace from existing flag and libs variables.

if test "$MPI_CPPFLAGS" != "" ; then
  MPI_CPPFLAGS=`echo $MPI_CPPFLAGS | sed 's/^[ ]*//;s/[ ]*$//'`
fi
if test "$MPI_LDFLAGS" != "" ; then
  MPI_LDFLAGS=`echo $MPI_LDFLAGS | sed 's/^[ ]*//;s/[ ]*$//'`
fi
if test "$MPI_LIBS" != "" ; then
  MPI_LIBS=`echo $MPI_LIBS | sed 's/^[ ]*//;s/[ ]*$//'`
fi

# If we do not use an MPI compiler wrapper, we must add compilation
# and link flags; we try to detect the correct flags to add.

if test "x$with_mpi" != "xno" ; then

  # try several tests for MPI

  # MPI Compiler wrapper test
  AC_MSG_CHECKING([for MPI (MPI compiler wrapper test)])
  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
  LIBS="$saved_LIBS $MPI_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                 [[ MPI_Init(0, (void *)0); ]])],
                 [ple_have_mpi=yes],
                 [ple_have_mpi=no])
  AC_MSG_RESULT($ple_have_mpi)

  # If failed, test for MPICH2 first
  if test "x$ple_have_mpi" = "xno"; then
    AC_MSG_CHECKING([for MPI (MPICH2 test)])
    # First try (without MPI-IO)
    case $host_os in
      freebsd*)
        MPI_LIBS="-lmpich -lopa -lmpl -lrt $PTHREAD_LIBS";;
      *)
        MPI_LIBS="-lmpich -lopa -lmpl -lrt -lpthread";;
    esac
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                    ple_have_mpi=yes],
                   [ple_have_mpi=no])
    if test "x$ple_have_mpi" = "xno"; then
      # Second try (without ROMIO)
      case $host_os in
        freebsd*)
          MPI_LIBS="-lmpich -lopa -lmpl $PTHREAD_LIBS";;
        *)
          MPI_LIBS="-lmpich -lopa -lmpl -lpthread";;
      esac
      LIBS="$saved_LIBS $MPI_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_Init(0, (void *)0); ]])],
                     [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                      ple_have_mpi=yes],
                     [ple_have_mpi=no])
    fi
    AC_MSG_RESULT($ple_have_mpi)
  fi


  # If failed, basic test (works with OpenMPI)
  if test "x$ple_have_mpi" = "xno"; then

    CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"

    # First, check for mpi.h header
    AC_CHECK_HEADERS([mpi.h],
                     [ple_have_mpi_header=yes],
                     [], 
                     [])

    if test $ple_have_mpi_header = no ; then
      unset ac_cv_header_mpi_h
      MPI_CPPFLAGS="-I/usr/include/mpi"
      CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
      AC_CHECK_HEADERS([mpi.h],
                       [ple_have_mpi_header=yes],
                       [], 
                       [])
    fi

    # Basic test
    AC_MSG_CHECKING([for MPI (basic test)])
    if test "$MPI_LIBS" = "" ; then
      MPI_LIBS="-lmpi $PTHREAD_LIBS"
    fi
    LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [ple_have_mpi=yes],
                   [ple_have_mpi=no])
    AC_MSG_RESULT($ple_have_mpi)
  fi

  if test "x$ple_have_mpi" = "xyes"; then
    AC_MSG_CHECKING([for MPI2 one-sided communication])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Win_free((void *)0); ]])],
                   [ple_have_mpi_one_sided=yes],
                   [ple_have_mpi_one_sided=no])
    AC_MSG_RESULT($ple_have_mpi_one_sided)
  else
    if test "x$with_mpi" != "xcheck" ; then
      AC_MSG_FAILURE([MPI support is requested, but test for MPI failed!])
    else
      AC_MSG_WARN([no MPI support])
    fi
    MPI_LIBS=""
  fi

  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

fi

AM_CONDITIONAL(HAVE_MPI, test x$ple_have_mpi = xyes)
AM_CONDITIONAL(HAVE_MPI_ONE_SIDED, test x$ple_have_mpi_one_sided = xyes)

CPPFLAGS="$saved_CPPFLAGS"
LDFLAGS="$saved_LDFLAGS"
LIBS="$saved_LIBS"

unset saved_CPPFLAGS
unset saved_LDFLAGS
unset saved_LIBS

AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)

])dnl


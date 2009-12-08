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

# CS_AC_TEST_MPI
#---------------
# optional MPI support (use CC=mpicc with configure if necessary)
# modifies or sets cs_have_mpi, MPI_CPPFLAGS, MPI_LDFLAGS, and MPI_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_MPI], [

saved_CPPFLAGS="$CPPFLAGS"
saved_LDFLAGS="$LDFLAGS"
saved_LIBS="$LIBS"

cs_have_mpi=no

AC_ARG_WITH(mpi,
            [AS_HELP_STRING([--with-mpi=PATH],
                            [specify prefix directory for MPI])],
            [if test "x$withval" = "x"; then
               with_mpi=yes
             fi],
            [with_mpi=check])

AC_ARG_WITH(mpi-exec,
            [AS_HELP_STRING([--with-mpi-exec=PATH],
                            [specify prefix directory for MPI executables])],
            [if test "x$with_mpi" = "xcheck"; then
               with_mpi=yes
             fi
             mpi_bindir="$with_mpi_exec"],
            [if test "x$with_mpi" != "xno" -a "x$with_mpi" != "xyes" \
	          -a "x$with_mpi" != "xcheck"; then
               mpi_bindir="$with_mpi/bin"
             fi])

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
             MPI_LDFLAGS="-L$with_mpi_lib"
             mpi_libdir="$with_mpi_lib"],
            [if test "x$with_mpi" != "xno" -a "x$with_mpi" != "xyes" \
	          -a "x$with_mpi" != "xcheck"; then
               MPI_LDFLAGS="-L$with_mpi/lib"
               mpi_libdir="$with_mpi/lib"
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

if test "x$with_mpi" != "xno" -a "x$cs_have_mpi" = "xno" ; then

  # try several tests for MPI

  # MPI Compiler wrapper test
  AC_MSG_CHECKING([for MPI (MPI compiler wrapper test)])
  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
  LIBS="$saved_LIBS $MPI_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                 [[ MPI_Init(0, (void *)0); ]])],
                 [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                  cs_have_mpi=yes],
                 [cs_have_mpi=no])
  AC_MSG_RESULT($cs_have_mpi)

  # If failed, basic test
  if test "x$cs_have_mpi" = "xno"; then
    # Basic test
    AC_MSG_CHECKING([for MPI (basic test)])
    if test "$MPI_LIBS" = "" ; then
      MPI_LIBS="-lmpi $PTHREAD_LIBS"
    fi
    CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                    cs_have_mpi=yes],
                   [cs_have_mpi=no])
    AC_MSG_RESULT($cs_have_mpi)
  fi

  # If failed, test for mpich
  if test "x$cs_have_mpi" = "xno"; then
    AC_MSG_CHECKING([for MPI (mpich test)])
    # First try (simplest)
    MPI_LIBS="-lmpich $PTHREAD_LIBS"
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                    cs_have_mpi=yes],
                   [cs_have_mpi=no])
    if test "x$cs_have_mpi" = "xno"; then
      # Second try (with lpmpich)
      MPI_LIBS="-Wl,-lpmpich -Wl,-lmpich -Wl,-lpmpich -Wl,-lmpich"
      LIBS="$saved_LIBS $MPI_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_Init(0, (void *)0); ]])],
                     [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                      cs_have_mpi=yes],
                     [cs_have_mpi=no])
    fi
    AC_MSG_RESULT($cs_have_mpi)
  fi

  # If failed, test for lam-mpi
  if test "x$cs_have_mpi" = "xno"; then
    AC_MSG_CHECKING([for MPI (lam-mpi test)])
    # First try (without MPI-IO)
    case $host_os in
      freebsd*)
        MPI_LIBS="-lmpi -llam $PTHREAD_LIBS";;
      *)
        MPI_LIBS="-lmpi -llam -lpthread";;
    esac
    LIBS="$saved_LIBS $MPI_LIBS"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Init(0, (void *)0); ]])],
                   [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                    cs_have_mpi=yes],
                   [cs_have_mpi=no])
    if test "x$cs_have_mpi" = "xno"; then
      # Second try (with MPI-IO)
      case $host_os in
        freebsd*)
          MPI_LIBS="-lmpi -llam -lutil -ldl $PTHREAD_LIBS";;
        *)
          MPI_LIBS="-lmpi -llam -lutil -ldl -lpthread";;
      esac
      LIBS="$saved_LIBS $MPI_LIBS"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_Init(0, (void *)0); ]])],
                     [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                      cs_have_mpi=yes],
                     [cs_have_mpi=no])
    fi
    AC_MSG_RESULT($cs_have_mpi)
  fi

  if test "x$cs_have_mpi" = "xno"; then
    if test "x$with_mpi" != "xcheck" ; then
      AC_MSG_FAILURE([MPI support is requested, but test for MPI failed!])
    else
      AC_MSG_WARN([no MPI support])
    fi
    MPI_LIBS=""
  else
    # Try to detect MPI variants as this may be useful for the run scripts to
    # determine the correct mpi startup syntax (especially when multiple
    # librairies are installed on the same machine).
    CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
    mpi_type=""
    if test "x$cs_ibm_bg_type" != "x" ; then
      if test "x$cs_ibm_bg_type" = "L" ; then
        mpi_type=BGL_MPI
      elif test "x$cs_ibm_bg_type" = "P" ; then
        mpi_type=BGP_MPI
      fi
    fi
    if test "x$mpi_type" = "x"; then
      AC_EGREP_CPP([mpich2],
                   [
                    #include <mpi.h>
                    #ifdef MPICH2
                    mpich2
                    #endif
                    ],
		    [mpi_type=MPICH2])
    fi
    if test "x$mpi_type" = "x"; then
      AC_EGREP_CPP([ompi],
                   [
                    #include <mpi.h>
                    #ifdef OMPI_MAJOR_VERSION
                    ompi
                    #endif
                    ],
		    [mpi_type=OpenMPI])
    fi
    if test "x$mpi_type" = "x"; then
      AC_EGREP_CPP([mpibull2],
                   [
                    #include <mpi.h>
                    #ifdef MPIBULL2_NAME
                    mpibull2
                    #endif
                    ],
		    [mpi_type=MPIBULL2])
    fi
    if test "x$mpi_type" = "x"; then
      AC_EGREP_CPP([lam_mpi],
                   [
                    #include <mpi.h>
                    #ifdef LAM_MPI
                    lam_mpi
                    #endif
                    ],
		    [mpi_type=LAM_MPI])
    fi
    if test "x$mpi_type" = "x"; then
      AC_EGREP_CPP([hp_mpi],
                   [
                    #include <mpi.h>
                    #ifdef HP_MPI
                    hp_mpi
                    #endif
                    ],
		    [mpi_type=HP_MPI])
    fi
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MPI, test x$cs_have_mpi = xyes)

AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(mpi_type)
AC_SUBST(mpi_bindir)
AC_SUBST(mpi_libdir)

])dnl


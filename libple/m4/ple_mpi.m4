dnl Copyright (C) 2005-2013 EDF
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
# and link flags; we will try to detect the correct flags to add.
# In any case, we will try to determine the type of MPI library
# using mpi.h before running link tests.

if test "x$with_mpi" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  # MPI Compiler wrapper test

  AC_MSG_CHECKING([for MPI (MPI compiler wrapper test)])
  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
  LIBS="$MPI_LIBS $saved_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                 [[ MPI_Init(0, (void *)0); ]])],
                 [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                  ple_have_mpi=yes
                  ple_have_mpi_header=yes],
                 [ple_have_mpi=no])
  AC_MSG_RESULT($ple_have_mpi)

  # If no wrapper was used, check for mpi.h header
  if test "x$ple_have_mpi" = "xno"; then
    CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
    AC_CHECK_HEADERS([mpi.h],
                     [ple_have_mpi_header=yes],
                     [], 
                     [])
  fi

  # If no MPI options were given, guess.

  if test "x$ple_have_mpi_header" = "xno" ; then
    unset ac_cv_header_mpi_h
    if test ! -z "$MPI_INCLUDE" ; then
      MPI_CPPFLAGS="-I$MPI_INCLUDE"
      # Also assume that a similar variable is defined for libraries
      MPI_LDFLAGS="-L$MPI_LIB"
    else
      MPI_CPPFLAGS="-I/usr/include/mpi"
    fi
    CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
    AC_CHECK_HEADERS([mpi.h],
                     [ple_have_mpi_header=yes],
                     [], 
                     [])
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

fi

# If no MPI headers are available, no use in pursuing MPI tests.

if test "x$ple_have_mpi_header" = "xyes" -a  "x$ple_have_mpi" = "xno" ; then

  # Now try to determine MPI variants, to test for the correct libraries
  # (especially when multiple librairies are installed on the same machine).

  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  mpi_h_type=""
  if test "x$mpi_h_type" = "x"; then
    AC_EGREP_CPP([ple_msmpi],
                 [
                  #include <mpi.h>
                  #ifdef MSMPI_VER
                  ple_msmpi
                  #endif
                  ],
                 [mpi_h_type=MSMPI])
  fi
  if test "x$mpi_h_type" = "x"; then
    AC_EGREP_CPP([ple_mpich],
                 [
                  #include <mpi.h>
                  #ifdef MPICH_NAME
                  #if (MPICH_NAME >= 3)
                  ple_mpich
                  #endif
                  #endif
                  ],
                 [mpi_h_type=MPICH])
  fi
  if test "x$mpi_h_type" = "x"; then
    AC_EGREP_CPP([ple_mpich2],
                 [
                  #include <mpi.h>
                  #ifdef MPICH2
                  ple_mpich2
                  #endif
                  ],
                 [mpi_h_type=MPICH2])
  fi
  if test "x$mpi_h_type" = "x"; then
    AC_EGREP_CPP([ple_ompi],
                 [
                  #include <mpi.h>
                  #ifdef OMPI_MAJOR_VERSION
                  ple_ompi
                  #endif
                  ],
                  [mpi_h_type=OpenMPI])
  fi
  if test "x$mpi_h_type" = "x"; then
    AC_EGREP_CPP([ple_lam_mpi],
                 [
                  #include <mpi.h>
                  #ifdef LAM_MPI
                  ple_lam_mpi
                  #endif
                  ],
                  [mpi_h_type=LAM_MPI])
  fi
  if test "x$mpi_h_type" = "x"; then
    AC_EGREP_CPP([ple_mpich1],
                 [
                  #include <mpi.h>
                  #ifdef MPICH
                  ple_mpich1
                  #endif
                  ],
                 [mpi_h_type=MPICH1])
  fi

  # Add a specific preprocessor directive to skip the MPI C++ bindings
  case $mpi_h_type in
    OpenMPI)         MPI_CPPFLAGS="$MPI_CPPFLAGS -DOMPI_SKIP_MPICXX" ;;
    MPICH | MPICH2)  MPI_CPPFLAGS="$MPI_CPPFLAGS -DMPICH_SKIP_MPICXX" ;;
  esac

  # If only MPI headers have been detected so far (i.e. we are
  # not using an MPI compiler wrapper), try to detect MPI libraries.

  if test "x$ple_have_mpi" = "xno" ; then

    LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"

    case $mpi_h_type in

      MPICH | MPICH2)
        AC_MSG_CHECKING([for MPICH-3 or MPICH2])
        # First try (with ROMIO)
        case $host_os in
          mingw32)
            MPI_LIBS="-lmpi";;
          freebsd*)
            MPI_LIBS="-lmpich -lopa -lmpl -lrt $PTHREAD_LIBS";;
          *)
            MPI_LIBS="-lmpich -lopa -lmpl -lrt -lpthread";;
        esac
        LIBS="$MPI_LIBS $saved_LIBS"
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
          LIBS="$MPI_LIBS $saved_LIBS"
          AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        ple_have_mpi=yes],
                       [ple_have_mpi=no])
        fi
        AC_MSG_RESULT($ple_have_mpi)
        ;;

      MSMPI)
        AC_MSG_CHECKING([for MSMPI])
        case $host_os in
          mingw32)
            MPI_LIBS="-lmsmpi";;
        esac
        LIBS="$MPI_LIBS $saved_LIBS"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        ple_have_mpi=yes],
                       [ple_have_mpi=no])
        AC_MSG_RESULT($ple_have_mpi)
        ;;

      MPICH1)
        AC_MSG_CHECKING([for MPICH1)])
        # First try (simplest)
        MPI_LIBS="-lmpich $PTHREAD_LIBS"
        LIBS="$MPI_LIBS $saved_LIBS"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        ple_have_mpi=yes],
                       [ple_have_mpi=no])
        if test "x$ple_have_mpi" = "xno"; then
          # Second try (with lpmpich)
          MPI_LIBS="-Wl,-lpmpich -Wl,-lmpich -Wl,-lpmpich -Wl,-lmpich"
          LIBS="$MPI_LIBS $saved_LIBS"
          AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                         [[ MPI_Init(0, (void *)0); ]])],
                         [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                          ple_have_mpi=yes],
                         [ple_have_mpi=no])
        fi
        AC_MSG_RESULT($ple_have_mpi)
        ;;

      LAM_MPI)
        AC_MSG_CHECKING([for LAM/MPI)])
        # First try (without MPI-IO)
        case $host_os in
          freebsd*)
            MPI_LIBS="-lmpi -llam $PTHREAD_LIBS";;
          *)
            MPI_LIBS="-lmpi -llam -lpthread";;
        esac
        LIBS="$MPI_LIBS $saved_LIBS"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        ple_have_mpi=yes],
                       [ple_have_mpi=no])
        if test "x$ple_have_mpi" = "xno"; then
          # Second try (with MPI-IO)
          case $host_os in
            freebsd*)
              MPI_LIBS="-lmpi -llam -lutil -ldl $PTHREAD_LIBS";;
            *)
              MPI_LIBS="-lmpi -llam -lutil -ldl -lpthread";;
          esac
          LIBS="$MPI_LIBS $saved_LIBS"
          AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                         [[ MPI_Init(0, (void *)0); ]])],
                         [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                          ple_have_mpi=yes],
                         [ple_have_mpi=no])
        fi
        AC_MSG_RESULT($ple_have_mpi)
        ;;

      *) # General case includes OpenMPI, whose dynamic libraries
         # usually make it easy to detect.
        AC_MSG_CHECKING([for MPI (basic test)])

        if test "$MPI_LIBS" = "" ; then
          MPI_LIBS="-lmpi $PTHREAD_LIBS"
        fi
        LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
        LIBS="$MPI_LIBS $saved_LIBS"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        ple_have_mpi=yes],
                       [ple_have_mpi=no])
        AC_MSG_RESULT($ple_have_mpi)
        ;;

    esac

  fi

  # MPI libraries should now have been detected

  if test "x$ple_have_mpi" = "xno"; then
    if test "x$with_mpi" != "xcheck" ; then
      AC_MSG_FAILURE([MPI support is requested, but test for MPI failed!])
    else
      AC_MSG_WARN([no MPI support])
    fi
    MPI_LIBS=""
  else
    # Try to detect some MPI 2 features
    AC_MSG_CHECKING([for MPI2 one-sided communication])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Win_free((void *)0); ]])],
                   [AC_DEFINE([HAVE_MPI_ONE_SIDED], 1, [MPI one-sided communication])
                    ple_have_mpi_one_sided=yes],
                   [ple_have_mpi_one_sided=no])
    AC_MSG_RESULT($ple_have_mpi_one_sided)

  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MPI, test x$ple_have_mpi = xyes)
AM_CONDITIONAL(HAVE_MPI_ONE_SIDED, test x$ple_have_mpi_one_sided = xyes)

AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)

])dnl


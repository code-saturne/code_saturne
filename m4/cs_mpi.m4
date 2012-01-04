dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2012 EDF S.A.
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
cs_have_mpi_header=no
cs_have_mpi_io=no
cs_have_mpi_one_sided=no

AC_ARG_ENABLE(mpi-io,
  [AS_HELP_STRING([--disable-mpi-io], [do not use MPI I/O when available])],
  [
    case "${enableval}" in
      yes) mpi_io=true ;;
      no)  mpi_io=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-mpi-io]) ;;
    esac
  ],
  [ mpi_io=true ]
)

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
# and link flags; we will try to detect the correct flags to add.
# In any case, we will try to determine the type of MPI library
# using mpi.h before running link tests.

if test "x$with_mpi" != "xno" ; then

  # MPI Compiler wrapper test

  AC_MSG_CHECKING([for MPI (MPI compiler wrapper test)])
  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
  LIBS="$saved_LIBS $MPI_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                 [[ MPI_Init(0, (void *)0); ]])],
                 [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                  cs_have_mpi=yes
                  cs_have_mpi_header=yes],
                 [cs_have_mpi=no])
  AC_MSG_RESULT($cs_have_mpi)

  # If no wrapper was used, check for mpi.h header
  if test "x$cs_have_mpi" = "xno"; then
    CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
    AC_CHECK_HEADERS([mpi.h],
                     [cs_have_mpi_header=yes],
                     [], 
                     [])
  fi

  # If no MPI options were given, guess.

  if test "x$cs_have_mpi_header" = "xno" ; then
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
                     [cs_have_mpi_header=yes],
                     [], 
                     [])
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

fi

# If no MPI headers are available, no use in pursuing MPI tests.

if test "x$cs_have_mpi_header" = "xyes" ; then

  # Now try to determine MPI variants, as this may be useful for the
  # run scripts to determine the correct mpi startup syntax
  # (especially when multiple librairies are installed on the same machine).

  CPPFLAGS="$saved_CPPFLAGS $MPI_CPPFLAGS"
  mpi_type=""
  if test "x$cs_ibm_bg_type" != "x" ; then
    if test "x$cs_ibm_bg_type" = "xL" ; then
      mpi_type=BGL_MPI
    elif test "x$cs_ibm_bg_type" = "xP" ; then
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
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([mpich1],
                 [
                  #include <mpi.h>
                  #ifdef MPICH
                  mpich1
                  #endif
                  ],
                 [mpi_type=MPICH1])
  fi

  # Add a specific preprocessor directive to skip the MPI C++ bindings
  case $mpi_type in
    OpenMPI) MPI_CPPFLAGS="$MPI_CPPFLAGS -DOMPI_SKIP_MPICXX" ;;
    MPICH2)  MPI_CPPFLAGS="$MPI_CPPFLAGS -DMPICH_SKIP_MPICXX" ;;
  esac

  # If only MPI headers have been detected so far (i.e. we are
  # not using an MPI compiler wrapper), try to detect MPI libraries.

  if test "x$cs_have_mpi" = "xno" ; then

    LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"

    case $mpi_type in

      MPICH2)
        AC_MSG_CHECKING([for MPICH2])
        # First try (with ROMIO)
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
                      cs_have_mpi=yes],
                     [cs_have_mpi=no])
        if test "x$cs_have_mpi" = "xno"; then
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
                        cs_have_mpi=yes
                        mpi_type=MPICH2],
                       [cs_have_mpi=no])
        fi
        AC_MSG_RESULT($cs_have_mpi)
        ;;

      MPICH1)
        AC_MSG_CHECKING([for MPICH1)])
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
        ;;

      
      *) # General case include OpenMPI, whose dynamic libraries
         # make it easy to detect.
        AC_MSG_CHECKING([for MPI (basic test)])

        if test "$MPI_LIBS" = "" ; then
          MPI_LIBS="-lmpi $PTHREAD_LIBS"
        fi
        LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"
        LIBS="$saved_LIBS $MPI_LIBS"
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        cs_have_mpi=yes],
                       [cs_have_mpi=no])
        AC_MSG_RESULT($cs_have_mpi)
        ;;

    esac

  fi

  # MPI libraries should now have been detected

  if test "x$cs_have_mpi" = "xno"; then
    if test "x$with_mpi" != "xcheck" ; then
      AC_MSG_FAILURE([MPI support is requested, but test for MPI failed!])
    else
      AC_MSG_WARN([no MPI support])
    fi
    MPI_LIBS=""
  else
    # Try to detect some MPI 2 features
    if test "x$mpi_io" = "xtrue"; then
      AC_MSG_CHECKING([for MPI I/O])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                     [[ MPI_File_close((void *)0); ]])],
                     [AC_DEFINE([HAVE_MPI_IO], 1, [MPI-IO support])
                      cs_have_mpi_io=yes],
                     [cs_have_mpi_io=no])
      AC_MSG_RESULT($cs_have_mpi_io)
    fi
    AC_MSG_CHECKING([for MPI2 one-sided communication])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Win_free((void *)0); ]])],
                   [AC_DEFINE([HAVE_MPI_ONE_SIDED], 1, [MPI one-sided communication])
                    cs_have_mpi_one_sided=yes],
                   [cs_have_mpi_one_sided=no])
    AC_MSG_RESULT($cs_have_mpi_one_sided)

  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(cs_have_mpi)
AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(mpi_type)
AC_SUBST(mpi_bindir)
AC_SUBST(mpi_libdir)

])dnl


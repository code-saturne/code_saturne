dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2019 EDF S.A.
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

cs_have_mpi=no
cs_have_mpi_header=no
cs_have_mpi_io=no
cs_have_mpi_one_sided=no
cs_have_mpi_neighbor_coll=no
cs_have_mpi_ibarrier=no

mpi_prefix=""

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
             mpi_includedir="$with_mpi_include"
            [if test "x$with_mpi" != "xno" -a "x$with_mpi" != "xyes" \
                  -a "x$with_mpi" != "xcheck"; then
               MPI_CPPFLAGS="-I$with_mpi/include"
               mpi_includedir="$with_mpi/include"
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
                  cs_have_mpi=yes
                  cs_have_mpi_header=yes],
                 [cs_have_mpi=no])
  AC_MSG_RESULT($cs_have_mpi)

  # If a wrapper was used, try to determine associated install path
  # (used to test for variants)

  if test "x$cs_have_mpi" = "xyes"; then

    if test "x$mpi_includedir" = "x" -o "x$mpi_libdir" = "x" ; then
      for arg in `$CC -show`; do
        case ${arg} in
          -I*)
            if test "x$mpi_includedir" = "x";
            then
              mpi_includedir=`echo ${arg} | sed -e 's/^-I//'`
            fi
            ;;
          -L*)
            if test "x$mpi_libdir" = "x";
            then
              mpi_libdir=`echo ${arg} | sed -e 's/^-L//'`
            fi
            ;;
          *)
            ;;
        esac
      done
    fi

  # If no wrapper was used, check for mpi.h header
  elif test "x$cs_have_mpi" = "xno"; then
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
  if test "x$cs_ibm_bg_type" = "xQ" ; then
    mpi_type=BGQ_MPI
  fi
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([cs_msmpi],
                 [
                  #include <mpi.h>
                  #ifdef MSMPI_VER
                  cs_msmpi
                  #endif
                  ],
                 [mpi_type=MSMPI])
  fi
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([cs_mpich],
                 [
                  #include <mpi.h>
                  #ifdef MPICH_NAME
                  #if (MPICH_NAME >= 3)
                  cs_mpich
                  #endif
                  #endif
                  ],
                 [mpi_type=MPICH])
  fi
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([cs_mpich2],
                 [
                  #include <mpi.h>
                  #ifdef MPICH2
                  cs_mpich2
                  #endif
                  ],
                 [mpi_type=MPICH2])
  fi
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([cs_ompi],
                 [
                  #include <mpi.h>
                  #ifdef OMPI_MAJOR_VERSION
                  cs_ompi
                  #endif
                  ],
                  [mpi_type=OpenMPI])
  fi
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([cs_mpibull2],
                 [
                  #include <mpi.h>
                  #ifdef MPIBULL2_NAME
                  cs_mpibull2
                  #endif
                  ],
                  [mpi_type=MPIBULL2])
  fi
  if test "x$mpi_type" = "x"; then
    AC_EGREP_CPP([cs_platform_mpi],
                 [
                  #include <mpi.h>
                  #ifdef PLATFORM_MPI
                  cs_platform_mpi
                  #endif
                  ],
                  [mpi_type=Platform_MPI])
  fi

  # Add a specific preprocessor directive to skip the MPI C++ bindings
  case $mpi_type in
    OpenMPI)         MPI_CPPFLAGS="$MPI_CPPFLAGS -DOMPI_SKIP_MPICXX" ;;
    MPICH | MPICH2)  MPI_CPPFLAGS="$MPI_CPPFLAGS -DMPICH_SKIP_MPICXX" ;;
  esac

  # Now try to determine if we are in fact using a variant MPI,
  # which does not define its own version macros in mpi.h but still uses its
  # own numbering (very ugly, but Intel and Bull do it).

  case $mpi_type in
    OpenMPI)         if test -d "${mpi_libdir}/bullxmpi" ; then
                       mpi_type=BullxMPI
                       AC_DEFINE([MPI_VENDOR_NAME], "BullxMPI", [MPI vendor name])
                     fi
                     ;;
    MPICH | MPICH2)  cs_mpisupport=""
                     if test -f "${mpi_libdir}/../mpisupport.txt" ; then
                       # mpi_libdir may point to lib sudirectory
                       cs_mpisupport="${mpi_libdir}/../mpisupport.txt"
                     elif test -f "${mpi_libdir}/../../mpisupport.txt" ; then
                       # mpi_libdir may point to intel64/lib sudirectory
                       cs_mpisupport="${mpi_libdir}/../../mpisupport.txt"
                     fi
                     if test "x$cs_mpisupport" != "x" ; then
                       grep "Intel(R) MPI" "$cs_mpisupport" > /dev/null 2>&1
                       if test $? = 0 ; then
                         mpi_type=Intel_MPI
                         AC_DEFINE([MPI_VENDOR_NAME], "Intel MPI", [MPI vendor name])
                       fi
                     fi
                     unset cs_mpisupport
                     ;;
  esac

  # If only MPI headers have been detected so far (i.e. we are
  # not using an MPI compiler wrapper), try to detect MPI libraries.

  if test "x$cs_have_mpi" = "xno" ; then

    LDFLAGS="$saved_LDFLAGS $MPI_LDFLAGS"

    case $mpi_type in

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
          LIBS="$MPI_LIBS $saved_LIBS"
          AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                       [[ MPI_Init(0, (void *)0); ]])],
                       [AC_DEFINE([HAVE_MPI], 1, [MPI support])
                        cs_have_mpi=yes],
                       [cs_have_mpi=no])
        fi
        AC_MSG_RESULT($cs_have_mpi)
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
                        cs_have_mpi=yes],
                       [cs_have_mpi=no])
        AC_MSG_RESULT($cs_have_mpi)
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
                        cs_have_mpi=yes],
                       [cs_have_mpi=no])
        AC_MSG_RESULT($cs_have_mpi)
        ;;

    esac

  fi

  # MPI libraries should now have been detected

  if test "x$cs_have_mpi" != "xno"; then
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
    AC_MSG_CHECKING([for MPI in place])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Allreduce(MPI_IN_PLACE, (void *)0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); ]])],
                   [AC_DEFINE([HAVE_MPI_IN_PLACE], 1, [MPI_IN_PLACE support])
                   have_mpi_in_place=yes],
                   [have_mpi_in_place=no])
    AC_MSG_RESULT($have_mpi_in_place)
    # Try to detect some MPI 3 features
    AC_MSG_CHECKING([for MPI Neighborhood collectives])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ void *sbuf, *rbuf;
                      int *scounts, *sdispls, *rcounts, *rdispls;
                      MPI_Comm comm;
                      MPI_Neighbor_alltoallv(sbuf, scounts, sdispls, MPI_INT,
                                             rbuf, rcounts, rdispls, MPI_INT,
                                             comm); ]])],
                   [AC_DEFINE([HAVE_MPI_NEIGHBOR_COLL], 1, [MPI neighborhood collectives])
                    cs_have_mpi_neighbor_coll=yes],
                   [cs_have_mpi_neighbor_coll=no])
      AC_MSG_RESULT($cs_have_mpi_neighbor_coll)
    AC_MSG_CHECKING([for MPI nonblocking barrier])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
                   [[ MPI_Request r;
                      MPI_Comm comm;
                      MPI_Ibarrier(comm, &r); ]])],
                   [AC_DEFINE([HAVE_MPI_IBARRIER], 1, [MPI nonblocking barrier])
                    cs_have_mpi_ibarrier=yes],
                   [cs_have_mpi_ibarrier=no])
      AC_MSG_RESULT($cs_have_mpi_ibarrier)
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  case $host_os in
   mingw32)
      mpi_prefix=`cygpath --path --windows "$with_mpi"`;;
    *)
      ;;
  esac

fi

if test "x$cs_have_mpi" = "xno" -a "x$with_mpi" != "xno"; then
  if test "x$with_mpi" != "xcheck" ; then
    AC_MSG_FAILURE([MPI support is requested, but test for MPI failed!])
  else
    AC_MSG_WARN([no MPI support])
  fi
  MPI_LIBS=""
fi

unset mpi_includedir

AC_SUBST(cs_have_mpi)
AC_SUBST(mpi_prefix, [${mpi_prefix}])
AC_SUBST(MPI_CPPFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(mpi_type)
AC_SUBST(mpi_bindir)
AC_SUBST(mpi_libdir)

])dnl


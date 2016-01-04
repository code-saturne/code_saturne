dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2016 EDF S.A.
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

# CS_AC_TEST_HDF5
#----------------
# modifies or sets cs_have_hdf5, HDF_CPPFLAGS, HDF_LDFLAGS, and HDF_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_HDF5], [

cs_have_hdf5=no
cs_have_hdf5_header=no
hdf5_prefix=""

AC_ARG_WITH(hdf5,
            [AS_HELP_STRING([--with-hdf5=PATH],
                            [specify prefix directory for HDF5])],
            [if test "x$withval" = "x"; then
               with_hdf5=yes
             fi],
            [with_hdf5=check])

AC_ARG_WITH(hdf5-include,
            [AS_HELP_STRING([--with-hdf5-include=PATH],
                            [specify directory for HDF5 include files])],
            [if test "x$with_hdf5" = "xcheck"; then
               with_hdf5=yes
             fi
             HDF5_CPPFLAGS="-I$with_hdf5_include"],
            [if test "x$with_hdf5" != "xno" -a "x$with_hdf5" != "xyes" \
	          -a "x$with_hdf5" != "xcheck"; then
               HDF5_CPPFLAGS="-I$with_hdf5/include"
             fi])

AC_ARG_WITH(hdf5-lib,
            [AS_HELP_STRING([--with-hdf5-lib=PATH],
                            [specify directory for HDF5 library])],
            [if test "x$with_hdf5" = "xcheck"; then
               with_hdf5=yes
             fi
             HDF5_LDFLAGS="-L$with_hdf5_lib"
             # Add the libdir to the runpath as HDF5 might not be libtoolized
             HDF5RUNPATH="-R$with_hdf5_lib"],
            [if test "x$with_hdf5" != "xno" -a "x$with_hdf5" != "xyes" \
	          -a "x$with_hdf5" != "xcheck"; then
               HDF5_LDFLAGS="-L$with_hdf5/lib"
               # Add the libdir to the runpath as HDF5 might not be libtoolized
               HDF5RUNPATH="-R$with_hdf5/lib"
             fi])


if test "x$with_hdf5" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${HDF5_CPPFLAGS}"

  # First, check for hdf5.h header
  AC_CHECK_HEADERS([hdf5.h],
                   [cs_have_hdf5_header=yes],
                   [],
                   [])

  if test $cs_have_hdf5_header = no ; then
    unset ac_cv_header_hdf5_h
    HDF5_CPPFLAGS="${HDF5_CPPFLAGS} ${MPI_CPPFLAGS}"
    CPPFLAGS="${saved_CPPFLAGS} ${HDF5_CPPFLAGS}"
    AC_CHECK_HEADERS([hdf5.h],
                     [cs_have_hdf5_header=yes],
                     [],
                     [])
  fi

  HDF5_LIBS="-lhdf5 $PTHREAD_LIBS"

  LDFLAGS="${LDFLAGS} ${HDF5_LDFLAGS}"
  LIBS="${HDF5_LIBS} ${LIBS}"

  # If HDF5 was built with MPI support, it might also be needed here

  AC_EGREP_CPP([cs_hdf5_parallel],
               [
                #include <hdf5.h>
                #ifdef H5_HAVE_PARALLEL
                #if (H5_HAVE_PARALLEL > 0)
                cs_hdf5_parallel
                #endif
                #endif
                ],
                [cs_hdf5_need_mpi=yes],
                [cs_hdf5_need_mpi=no])

  if test "x$cs_hdf5_need_mpi" = "xyes" ; then
    HDF5_CPPFLAGS_MPI=$MPI_CPPFLAGS
    HDF5_LDFLAGS_MPI=$MPI_LDFLAGS
    HDF5_LIBS_MPI=$MPI_LIBS
    CPPFLAGS="${CPPFLAGS} ${HDF5_CPPFLAGS_MPI}"
    LDFLAGS="${LDFLAGS} ${HDF5_LDFLAGS_MPI}"
    LIBS="${HDF5_LIBS} ${HDF5_LIBS_MPI}"
  fi

  # Now check library

  AC_CHECK_LIB(hdf5, H5Fopen,
               [ AC_DEFINE([HAVE_HDF5], 1, [HDF5 file support])
                 cs_have_hdf5=yes
               ],
               [if test "x$with_hdf5" != "xcheck" ; then
                  AC_MSG_FAILURE([HDF5 support is requested, but test for HDF5 failed!])
                else
                  AC_MSG_WARN([no HDF5 file support])
                fi
               ],
               )

  if test "x$cs_have_hdf5" = "xno"; then
    HDF5_LIBS=""
    HDF5_CPPFLAGS_MPI=""
    HDF5_LDFLAGS_MPI=""
    HDF5_LIBS_MPI=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  case $host_os in
    mingw32)
      hdf5_prefix=`cygpath --path --windows "$with_hdf5"`;;
    *)
      ;;
  esac

fi

AC_SUBST(cs_have_hdf5)
AC_SUBST(hdf5_prefix, [${hdf5_prefix}])
AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_LIBS)
AC_SUBST(HDF5RUNPATH)

AC_SUBST(HDF5_CPPFLAGS_MPI)
AC_SUBST(HDF5_LDFLAGS_MPI)
AC_SUBST(HDF5_LIBS_MPI)

])dnl


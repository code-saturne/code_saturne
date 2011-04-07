dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2010 EDF S.A., France
dnl
dnl   The Code_Saturne CFD tool is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne CFD tool is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne CFD tool; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_TEST_HDF5
#----------------
# modifies or sets cs_have_hdf5, HDF_CPPFLAGS, HDF_LDFLAGS, and HDF_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_HDF5], [

cs_have_hdf5=no

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
             HDF5_LDFLAGS="-L$with_hdf5_lib"],
            [if test "x$with_hdf5" != "xno" -a "x$with_hdf5" != "xyes" \
	          -a "x$with_hdf5" != "xcheck"; then
               HDF5_LDFLAGS="-L$with_hdf5/lib"
             fi])


if test "x$with_hdf5" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  HDF5_LIBS="-lhdf5 $PTHREAD_LIBS"
  
  CPPFLAGS="${CPPFLAGS} ${HDF5_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${HDF5_LDFLAGS}"
  LIBS="${LIBS} ${HDF5_LIBS}"

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
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(cs_have_hdf5)
AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_LIBS)

])dnl


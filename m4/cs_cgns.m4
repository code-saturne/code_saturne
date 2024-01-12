dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2024 EDF S.A.
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

# CS_AC_TEST_CGNS
#----------------
# modifies or sets cs_have_cgns, CGNS_CPPFLAGS, CGNS_LDFLAGS, and CGNS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CGNS], [

cs_have_cgns=no
cs_have_cgns_headers=no
cgns_prefix=""

AC_ARG_VAR([CGNS_ROOT_DIR], [CGNS root directory (superseded by --with-cgns=PATH)])

AC_ARG_WITH(cgns,
            [AS_HELP_STRING([--with-cgns=PATH],
                            [specify prefix directory for CGNS])],
            [if test "x$withval" = "xyes"; then
               if test "x$CGNS_ROOT_DIR" != "x"; then
                 with_cgns=$CGNS_ROOT_DIR
               fi
             elif test "x$withval" = "xsalome"; then
               cs_salome_cgns_root_dir=`(/bin/bash -c "unset LD_LIBRARY_PATH ; $SALOMEENVCMD > /dev/null 2>&1 ; echo $CGNS_ROOT_DIR")`
               if test -z "cs_salome_cgns_root_dir"; then
                 AC_MSG_FAILURE([no SALOME path information for CGNS (needed by --with-cgns=salome)!])
               else
                 with_cgns=$cs_salome_cgns_root_dir
               fi
               unset cs_salome_cgns_root_dir
             fi],
            [if test "x$CGNS_ROOT_DIR" != "x"; then
               with_cgns=$CGNS_ROOT_DIR
             else
               with_cgns=check
             fi])

AC_ARG_WITH(cgns-include,
            [AS_HELP_STRING([--with-cgns-include=PATH],
                            [specify directory for CGNS include files])],
            [if test "x$with_cgns" = "xcheck"; then
               with_cgns=yes
             fi
             CGNS_CPPFLAGS="-I$with_cgns_include"],
            [if test "x$with_cgns" != "xno" -a "x$with_cgns" != "xyes" \
	          -a "x$with_cgns" != "xcheck"; then
               CGNS_CPPFLAGS="-I$with_cgns/include"
             fi])

AC_ARG_WITH(cgns-lib,
            [AS_HELP_STRING([--with-cgns-lib=PATH],
                            [specify directory for CGNS library])],
            [if test "x$with_cgns" = "xcheck"; then
               with_cgns=yes
             fi
             CGNS_LDFLAGS="-L$with_cgns_lib"
             # Add the libdir to the runpath for the preprocessor build
             CGNSRUNPATH="${LDRPATH}${with_cgns_lib}"],
            [if test "x$with_cgns" != "xno" -a "x$with_cgns" != "xyes" \
	          -a "x$with_cgns" != "xcheck"; then
               CGNS_LDFLAGS="-L$with_cgns/lib"
               # Add the libdir to the runpath for the preprocessor build
               CGNSRUNPATH="${LDRPATH}${with_cgns}/lib"
             fi])


if test "x$with_cgns" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CGNS_LIBS="-lcgns"
  CPPFLAGS="${CPPFLAGS} ${CGNS_CPPFLAGS} ${HDF5_CPPFLAGS_MPI}"
  LDFLAGS="${LDFLAGS} ${CGNS_LDFLAGS} ${HDF5_LDFLAGS} ${HDF5_LDFLAGS_MPI}"
  LIBS="${CGNS_LIBS} ${HDF5_LIBS} ${HDF5_LIBS_MPI} ${LIBS}"

  # Check that a header file exists and that the version is compatible
  #-------------------------------------------------------------------

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#undef HAVE_MPI
#include <stdio.h>
#include <cgnslib.h>]],
[[#if CGNS_VERSION < 3100
# error CGNS version >= 3.0 not found
#endif
]])],
                    [AC_MSG_RESULT([CGNS >= 3.1.0 headers found])
                     cs_have_cgns_headers=yes
                    ],
                    [AC_MSG_RESULT([CGNS >= 3.1.0 headers not found])
                    ])

  if test "x$cs_have_cgns_headers" = "xno"; then

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#undef HAVE_MPI
#include <stdio.h>
#include <cgnslib.h>]],
[[#if CGNS_VERSION <= 2400
# error CGNS version >= 2.4 tested here
#endif
]])],
                      [AC_MSG_FAILURE([CGNS < 3.1 headers found, but CGNS 3.1 or above is required.
If you do not need CGNS format support, you may use the --without-cgns configure option.
Otherwise, you need to provide a CGNS 3.1 library and development headers.])
                      ],
                      [])

  fi # end of test on CGNS 2 headers

  # Check for a CGNS 3.1+ library
  #------------------------------

  if test "x$cs_have_cgns_headers" = "xyes"; then

    AC_CHECK_LIB(cgns, cg_coord_partial_write,
                 [ AC_DEFINE([HAVE_CGNS], 1, [CGNS file support])
                   cs_have_cgns=yes
                 ],
                 [])

  fi

  if test "x$cs_have_cgns" != "xyes"; then
    CGNS_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  # Report CGNS support
  #-------------------

  if test "x$cs_have_cgns" = "xno" ; then
    if test "x$with_cgns" != "xcheck" ; then
      AC_MSG_FAILURE([CGNS support is requested, but test for CGNS failed!])
    else
      AC_MSG_WARN([no CGNS file support])
    fi
  fi

  case $host_os in
    mingw64)
      cgns_prefix=`cygpath --path --windows "$with_cgns"`;;
    *)
      ;;
  esac

fi

unset cs_have_cgns_headers

AM_CONDITIONAL(HAVE_CGNS, test x$cs_have_cgns = xyes)

AC_SUBST(cs_have_cgns)
AC_SUBST(cgns_prefix, [${cgns_prefix}])
AC_SUBST(CGNS_CPPFLAGS)
AC_SUBST(CGNS_LDFLAGS)
AC_SUBST(CGNS_LIBS)
AC_SUBST(CGNSRUNPATH)

])dnl


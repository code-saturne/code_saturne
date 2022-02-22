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

# CS_AC_TEST_MED
#---------------
# modifies or sets cs_have_med, MED_CPPFLAGS, MED_LDFLAGS, and MED_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_MED], [

cs_have_med=no
cs_have_med_mpi=no
cs_have_med_headers=no
cs_have_med_link_cxx=no
med_prefix=""

# Configure options
#------------------

AC_ARG_WITH(med,
            [AS_HELP_STRING([--with-med=PATH],
                            [specify prefix directory for MED])],
            [if test "x$withval" = "x"; then
               with_med=yes
             elif test "x$withval" = "xsalome"; then
               if test -z "$MEDHOME"; then
                 AC_MSG_FAILURE([no SALOME path information for MED (needed by --with-med=salome)!])
               else
                 with_med=$MEDHOME
               fi
             fi],
            [with_med=check])

AC_ARG_WITH(med-include,
            [AS_HELP_STRING([--with-med-include=PATH],
                            [specify directory for MED include files])],
            [if test "x$with_med" = "xcheck"; then
               with_med=yes
             fi
             MED_CPPFLAGS="-I$with_med_include"],
            [if test "x$with_med" != "xno" -a "x$with_med" != "xyes" \
	          -a "x$with_med" != "xcheck"; then
               MED_CPPFLAGS="-I$with_med/include"
             fi])

AC_ARG_WITH(med-lib,
            [AS_HELP_STRING([--with-med-lib=PATH],
                            [specify directory for MED library])],
            [if test "x$with_med" = "xcheck"; then
               with_med=yes
             fi
             MED_LDFLAGS="-L$with_med_lib"
             # Add the libdir to the runpath as MED libtool .la files might not be present
             MEDRUNPATH="-R$with_med_lib"],
            [if test "x$with_med" != "xno" -a "x$with_med" != "xyes" \
	          -a "x$with_med" != "xcheck"; then
               MED_LDFLAGS="-L$with_med/lib"
               # Add the libdir to the runpath as MED libtool .la files might not be present
               MEDRUNPATH="-R$with_med/lib"
             fi])

if test "x$with_med" != "xno" -a "x$cs_have_hdf5" = "xno"; then
  if test "x$with_med" = "xcheck"; then
    with_med=no
    AC_MSG_WARN([no hdf5 library found; will not search for MED])
  else
    AC_MSG_ERROR([no hdf5 library found; required for MED])
  fi
fi

if test "x$with_med" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  MED_LIBS="-lmedC"

  CPPFLAGS="${CPPFLAGS} ${MED_CPPFLAGS} ${HDF5_CPPFLAGS} ${HDF5_CPPFLAGS_MPI}"
  LDFLAGS="${MED_LDFLAGS} ${HDF5_LDFLAGS} ${HDF5_LDFLAGS_MPI} ${LDFLAGS}"
  LIBS="${MED_LIBS} ${HDF5_LIBS} ${HDF5_LIBS_MPI} ${LIBS}"

  # Check that MED header files exist and that the version is compatible
  #---------------------------------------------------------------------

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#undef HAVE_MPI
#include <med.h>]],
[[#if !defined(MED_MAJOR_NUM)
# error MED_MAJOR_NUM not defined
#endif
#if MED_MAJOR_NUM < 3
# error MED version >= 3.0.0 required
#endif
]])],
                    [AC_MSG_RESULT([MED >= 3.0 headers found])
                     cs_have_med_headers=yes
                    ],
                    [AC_MSG_RESULT([MED >= 3.0 headers not found])
                    ])

  if test "x$cs_have_med_headers" = "xno"; then

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <med.h>]],
[[#if !defined(MED_NUM_MAJEUR)
# error MED_NUM_MAJEUR not defined
#endif
#if MED_NUM_MAJEUR != 2 || MED_NUM_MINEUR != 3
# error MED version > 2.3 tested here
#endif
]])],
                      [AC_MSG_FAILURE([MED 2.3 headers found, but MED 3.0 or above is required.
If you do not need MED format support, you may use the --without-med configure option.
Otherwise, you need to provide a MED 3.0 library and development headers.])
                      ],
                      [])

  fi # end of test on cs_have_med_headers

  # Check for a MED 3.x library
  #----------------------------

  if test "x$cs_have_med_headers" = "xyes"; then

    AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#undef HAVE_MPI
#include <med.h>]],
[[(void)MEDfamilyCr(0, NULL, NULL, 0, 0, NULL);]])
                   ],
                   [ AC_DEFINE([HAVE_MED], 1, [MED file support])
                     cs_have_med=yes
                   ],
                   [ AC_MSG_WARN([no MED file support with C only link]) ],
                  )

    if test "x$cs_have_med" = "xno"; then

      # try linking with C++ in case of static MED library

      AC_LANG_PUSH(C++)
      AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#undef HAVE_MPI
#include <med.h>]],
[[(void)MEDfamilyCr(0, NULL, NULL, 0, 0, NULL);]])
                     ],
                     [ AC_DEFINE([HAVE_MED], 1, [MED file support])
                       cs_have_med=yes; cs_have_med_link_cxx=yes
                     ],
                     [ AC_MSG_WARN([no MED file support])
                     ],
                     )
      AC_LANG_POP(C++)

    fi

    # Check for parallel MED options
    #-------------------------------

    if test "x$cs_have_mpi" = "xyes" -a "x$cs_have_med" = "xyes"; then

      if test "x$cs_have_med_link_cxx" = "xno"; then
        AC_LANG_PUSH(C++)
      fi

      AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <mpi.h>
#include <med.h>]],
[[(void)MEDparFileOpen(NULL, MED_ACC_RDONLY, MPI_COMM_NULL, MPI_INFO_NULL);]])
                     ],
                     [ AC_DEFINE([HAVE_MED_MPI], 1, [MED file MPI support])
                       cs_have_med_mpi=yes
                     ],
                     [ AC_MSG_WARN([no MED file MPI support]) ],
                    )

      if test "x$cs_have_med_link_cxx" = "xno"; then
        AC_LANG_POP(C++)
      fi

    fi

  fi

  # Report MED support
  #-------------------

  if test "x$cs_have_med" = "xno" ; then
    if test "x$with_med" != "xcheck" ; then
      AC_MSG_FAILURE([MED support is requested, but test for MED failed!])
    else
      AC_MSG_WARN([no MED file support])
    fi
  fi

  if test "x$cs_have_med" = "xno"; then
    MED_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

  case $host_os in
    mingw64)
      med_prefix=`cygpath --path --windows "$with_med"`;;
    *)
      ;;
  esac
fi

unset cs_have_med_headers

AM_CONDITIONAL(HAVE_MED, test x$cs_have_med = xyes)

AC_SUBST(cs_have_med)
AC_SUBST(med_prefix, [${med_prefix}])
AC_SUBST(MED_CPPFLAGS)
AC_SUBST(MED_LDFLAGS)
AC_SUBST(MED_LIBS)
AC_SUBST(MEDRUNPATH)

])dnl


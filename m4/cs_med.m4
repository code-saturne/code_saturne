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

# CS_AC_TEST_MED
#---------------
# modifies or sets cs_have_med, MED_CPPFLAGS, MED_LDFLAGS, and MED_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_MED], [

cs_have_med=no
cs_have_med_headers=no
cs_have_med_link_cxx=no

AC_ARG_WITH(med,
            [AS_HELP_STRING([--with-med=PATH],
                            [specify prefix directory for MED])],
            [if test "x$withval" = "x"; then
               with_med=yes
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
             MED_LDFLAGS="-L$with_med_lib"],
            [if test "x$with_med" != "xno" -a "x$with_med" != "xyes" \
	          -a "x$with_med" != "xcheck"; then
               MED_LDFLAGS="-L$with_med/lib"
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

  MED_RPATH=
  MED_LIBS="-lmedC"

  if test "x$MED_LDFLAGS" != "x" ; then
    MED_RPATH="-Wl,-rpath -Wl,${MED_LDFLAGS}"
  fi

  if test "x$with_med_dep_libs" != "x" ; then
    med_dep_list="`echo $with_med_dep_libs | sed 's/,/ /g'`"
    for arg in $med_dep_list; do
      MED_LIBS="${MED_LIBS} -l${arg}"
    done
  fi

  if test "x$with_med_dep_dirs" != "x" ; then
    med_dep_list="`echo $with_med_dep_dirs | sed 's/,/ /g'`"
    for arg in $med_dep_list; do
      MED_LDFLAGS="${MED_LDFLAGS} -L${arg}"
    done
  fi

  CPPFLAGS="${CPPFLAGS} ${MED_CPPFLAGS} ${HDF5_CPPFLAGS}"
  LDFLAGS="${MED_LDFLAGS} ${HDF5_LDFLAGS} ${LDFLAGS}"
  LIBS="${MED_LIBS} ${HDF5_LIBS} ${LIBS}"

# Check that MED header files exist and that the version is compatible

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <med.h>]],
[[#if !defined(MED_NUM_MAJEUR)
# error MED_NUM_MAJEUR not defined
#endif
#if MED_NUM_MAJEUR == 2 && MED_NUM_MINEUR == 3 && MED_NUM_RELEASE < 3
# error MED version > 2.3.4 required
#endif
]])],
                    [AC_MSG_RESULT([MED >= 2.3.4 headers found])
                     cs_have_med_headers=yes
                    ],
                    [AC_MSG_RESULT([MED >= 2.3.4 headers not found])
                    ])

  if test "x$cs_have_med_headers" = "xyes"; then

    AC_CHECK_LIB(medC, MEDfamCr, 
                 [ AC_DEFINE([HAVE_MED], 1, [MED file support])
                   cs_have_med=yes
                 ], 
                 [ AC_MSG_WARN([no MED file support with C only link])
                 ],
                 )

    if test "x$cs_have_med" = "xno"; then
  
      # try linking with C++ in case of static MED library

      AC_LANG_PUSH(C++)
      AC_CHECK_LIB(med, MEDfamCr, 
                   [ AC_DEFINE([HAVE_MED], 1, [MED file support])
                     cs_have_med=yes; cs_have_med_link_cxx=yes
                   ], 
                   [ AC_MSG_WARN([no MED file support])
                   ],
                   )
      AC_LANG_POP(C++)

    fi

    if test "x$cs_have_med" = "xno" ; then
      if test "x$with_med" != "xcheck" ; then
        AC_MSG_FAILURE([MED support is requested, but test for MED failed!])
      else
        AC_MSG_WARN([no MED file support])
      fi
    fi

  fi # end of test on cs_have_med_headers

  if test "x$cs_have_med" = "xno"; then
    MED_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MED, test x$cs_have_med = xyes)

AC_SUBST(cs_have_med)
AC_SUBST(MED_CPPFLAGS)
AC_SUBST(MED_LDFLAGS)
AC_SUBST(MED_LIBS)

])dnl


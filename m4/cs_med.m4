dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2010-2011 EDF S.A., France
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
cs_have_med_mpi=no
cs_have_med2_headers=no
cs_have_med3_headers=no
cs_have_med_link_cxx=no

# Configure options
#------------------

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

  MED_LIBS="-lmedC"

  CPPFLAGS="${CPPFLAGS} ${MED_CPPFLAGS} ${HDF5_CPPFLAGS}"
  LDFLAGS="${MED_LDFLAGS} ${HDF5_LDFLAGS} ${LDFLAGS}"
  LIBS="${MED_LIBS} ${HDF5_LIBS} ${LIBS}"

  # Check that MED header files exist and that the version is compatible
  #---------------------------------------------------------------------

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#undef HAVE_MPI
#include <med.h>]],
[[#if !defined(MED_MAJOR_NUM)
# error MED_MAJOR_NUM not defined
#endif
#if MED_MAJOR_NUM == 2 && MED_MINOR_NUM < 9
# error MED version >= 2.9.0 required
#endif
]])],
                    [AC_MSG_RESULT([MED >= 2.9.0 headers found])
                     cs_have_med3_headers=yes
                    ],
                    [AC_MSG_RESULT([MED >= 2.9.0 headers not found])
                    ])

  if test "x$cs_have_med3_headers" = "xno"; then

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
                       cs_have_med2_headers=yes
                      ],
                      [AC_MSG_RESULT([MED >= 2.3.4 headers not found])
                      ])

  fi # end of test on cs_have_med3_headers

  # Check for a MED 3.x library
  #----------------------------

  if test "x$cs_have_med3_headers" = "xyes"; then

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

      CPPFLAGS="${CPPFLAGS} ${MPI_CPPFLAGS}"
      LDFLAGS="${MED_LDFLAGS} ${HDF5_LDFLAGS} ${MPI_LDFLAGS} ${saved_LDFLAGS}"
      LIBS="${MED_LIBS} ${HDF5_LIBS} ${MPI_LIBS} ${saved_LIBS}"

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

  # Check for a MED 2.3.x library
  #----------------------------

  elif test "x$cs_have_med2_headers" = "xyes"; then

    AC_CHECK_LIB(med, MEDfamCr, 
                 [ AC_DEFINE([HAVE_MED], 1, [MED file support])
                   cs_have_med=yes
                 ], 
                 [ AC_MSG_WARN([no MED file support])
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

fi

unset cs_have_med2_headers
unset cs_have_med3_headers

AM_CONDITIONAL(HAVE_MED, test x$cs_have_med = xyes)

AC_SUBST(cs_have_med)
AC_SUBST(MED_CPPFLAGS)
AC_SUBST(MED_LDFLAGS)
AC_SUBST(MED_LIBS)

])dnl


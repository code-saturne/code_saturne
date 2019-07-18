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

# CS_AC_TEST_MEDCOUPLING
#----------------------
# modifies or sets cs_have_medcoupling, MEDCOUPLING_CPPFLAGS, MEDCOUPLING_LDFLAGS,
# and MEDCOUPLING_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_MEDCOUPLING], [

cs_have_medcoupling=no
cs_have_medcoupling_loader=no
cs_have_paramedmem=no

# Configure options
#------------------

AC_ARG_WITH(medcoupling,
            [AS_HELP_STRING([--with-medcoupling=PATH],
                            [specify directory for MEDCoupling and ParaMEDMEM])],
            [if test "x$withval" = "x"; then
               with_medcoupling=yes
             fi],
            [with_medcoupling=check])

AC_ARG_ENABLE(medcoupling-as-plugin,
  [AS_HELP_STRING([--enable-medcoupling-as-plugin], [use MEDCoupling as plugin])],
  [
    case "${enableval}" in
      yes) cs_have_plugin_medcoupling=yes ;;
      no)  cs_have_plugin_medcoupling=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-medcoupling-as-plugin]) ;;
    esac
  ],
  [ cs_have_plugin_medcoupling=no ]
)

if test x$cs_have_dlloader = xno -o x$enable_shared = xno ; then
  cs_have_plugin_medcoupling=no
fi

if test "$with_medcoupling" != no ; then

  if test "$with_medcoupling" = yes -o "$with_medcoupling" = check ; then
    if test ! -z "$MEDCOUPLING_ROOT_DIR"; then
      MEDCOUPLING=$MEDCOUPLING_ROOT_DIR
    fi
  else
    if test -d "$with_medcoupling" ; then
      MEDCOUPLING="$with_medcoupling"
    else
      AC_MSG_FAILURE([directory specified by --with-medcoupling=$with_medcoupling does not exist!])
    fi
  fi

fi

if test "x$with_medcoupling" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test "x$with_medcoupling" != "x" ; then
    MEDCOUPLING_CPPFLAGS="-I$MEDCOUPLING/include"
    MEDCOUPLING_LDFLAGS="-L$MEDCOUPLING/lib"
    # Add the libdir to the runpath as libtool does not do this for modules
    MEDCOUPLINGRUNPATH="-R$MEDCOUPLING/lib"
  fi

  AC_LANG_PUSH([C++])

  # First check for MEDCoupling file support

  CPPFLAGS="${MPI_CPPFLAGS} ${MEDCOUPLING_CPPFLAGS} ${MED_CPPFLAGS} ${HDF5_CPPFLAGS} ${CPPFLAGS}"

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <MEDCouplingFieldDouble.hxx>
#include <MEDCouplingPartDefinition.hxx>
#include <MEDLoader.hxx>
#include <MEDFileField1TS.hxx>]],
[[using namespace MEDCoupling;
  std::string f_name;
  MEDCouplingField *f = ReadFieldCell("path", "name", 0, f_name, 0, 0);]])
                    ],
                    [ cs_have_medcoupling_loader=yes ],
                    [ AC_MSG_WARN([no MEDCoupling file support]) ],
                  )

  cs_medcoupling_l0=

  # Check for MEDCoupling library
  #-------------------------------

  if test "$cs_have_medcoupling" = "no"; then

    # Check for minimal MEDCoupling

    if test "$cs_have_medcoupling_loader" = "no"; then

      MEDCOUPLING_LIBS="-lmedcoupling -linterpkernel -lmedcouplingremapper"

      LDFLAGS="${MEDCOUPLING_LDFLAGS} ${LDFLAGS}"
      LIBS="${MEDCOUPLING_LIBS} ${LIBS}"

      AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <MEDCouplingUMesh.hxx>]],
[[using namespace MEDCoupling;
MEDCouplingUMesh *m = MEDCouplingUMesh::New();]])
                     ],
                     [ AC_DEFINE([HAVE_MEDCOUPLING], 1, [MEDCoupling support])
                       cs_have_medcoupling=yes
                     ],
                     [ ],
                    )

    # Check for regular MEDCoupling

    else

      MEDCOUPLING_LIBS="-lmedcoupling -linterpkernel -lmedcouplingremapper -lmedloader"

      LDFLAGS="${MEDCOUPLING_LDFLAGS} ${MED_LDFLAGS} ${HDF5_LDFLAGS} ${LDFLAGS}"
      LIBS="${MEDCOUPLING_LIBS} ${MED_LIBS} ${HDF5_LIBS} ${LIBS}"

      AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <MEDCouplingFieldDouble.hxx>
#include <MEDLoader.hxx>]],
[[using namespace MEDCoupling;
  std::string f_name;
  MEDCouplingField *f = ReadFieldCell("path", "name", 0, f_name, 0, 0);]])
                     ],
                     [ AC_DEFINE([HAVE_MEDCOUPLING], 1, [MEDCoupling support])
                       cs_have_medcoupling=yes
                     ],
                     [ ],
                    )
    fi

    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi

  if test "$cs_have_medcoupling" = "yes"; then
    if test "$cs_have_medcoupling_loader" = "yes"; then
      AC_DEFINE([HAVE_MEDCOUPLING_LOADER], 1, [MEDCoupling with loader support])
    fi
  else
    cs_have_medcoupling_loader="no"
  fi

  # Now check for MEDCoupling MPI support

  if test "$cs_have_medcoupling" = "yes" -a "$cs_have_mpi" = "yes"; then

    CPPFLAGS="${MPI_CPPFLAGS} ${MEDCOUPLING_CPPFLAGS} ${CPPFLAGS}"

    if test "$cs_have_medcoupling_loader" = "yes"; then
      cs_paramedmem_libs="-lparamedmem -lparamedloader"
    else
      cs_paramedmem_libs="-lparamedmem"
    fi

    if test "x$cs_have_paramedmem" = "xno" ; then

      if test "$cs_have_medcoupling_loader" = "no"; then
        LDFLAGS="${MEDCOUPLING_LDFLAGS} ${MPI_LDFLAGS} ${LDFLAGS}"
        LIBS="${cs_paramedmem_libs} ${MEDCOUPLING_LIBS} ${MPI_LIBS} ${LIBS}"
      else
        LDFLAGS="${MEDCOUPLING_LDFLAGS} ${MED_LDFLAGS} ${HDF5_LDFLAGS} ${MPI_LDFLAGS} ${LDFLAGS}"
        LIBS="${cs_paramedmem_libs} ${MEDCOUPLING_LIBS} ${MED_LIBS} ${HDF5_LIBS} ${MPI_LIBS} ${LIBS}"
      fi

      AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <InterpKernelDEC.hxx>
#include <set>]],
[[using namespace MEDCoupling;
int procs_source_c[1]={0};
std::set<int> procs_source(procs_source_c, procs_source_c+1);
int procs_target_c[1]={1};
std::set<int> procs_target(procs_target_c, procs_target_c+1);
InterpKernelDEC *dec = new InterpKernelDEC(procs_source, procs_target);]])
                       ],
                       [ AC_DEFINE([HAVE_PARAMEDMEM], 1, [ParaMEDMEM support])
                         cs_have_paramedmem=yes
                       ],
                       [ ],
                      )

      if test "x$cs_have_paramedmem" = "xyes"; then
        MEDCOUPLING_LIBS="${cs_paramedmem_libs} ${MEDCOUPLING_LIBS}"
      fi

      LDFLAGS="$saved_LDFLAGS"
      LIBS="$saved_LIBS"

    fi

    if test "x$cs_have_paramedmem" != "xyes"; then
      AC_MSG_WARN([no ParaMEDMEM support])
    fi

  fi

  CPPFLAGS="$saved_CPPFLAGS"

  AC_LANG_POP([C++])

  # Report MEDCOUPLING support
  #-------------------

  if test "x$cs_have_medcoupling" = "xyes" ; then
    if test x$cs_have_plugin_medcoupling = xyes ; then
      AC_DEFINE([HAVE_PLUGIN_MEDCOUPLING], 1, [MEDCoupling support as plugin])
    fi
  elif test "x$cs_have_medcoupling" = "xno" ; then
    if test "x$with_medcoupling" != "xcheck" ; then
      AC_MSG_FAILURE([MEDCoupling support is requested, but test for MEDCoupling failed!])
    fi
  fi

  if test "x$cs_have_medcoupling" = "xno"; then
    MEDCOUPLING_LIBS=""
  fi

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MEDCOUPLING, test x$cs_have_medcoupling = xyes)
AM_CONDITIONAL(HAVE_MEDCOUPLING_LOADER, test x$cs_have_medcoupling_loader = xyes)
AM_CONDITIONAL(HAVE_PARAMEDMEM, test x$cs_have_paramedmem = xyes)
AM_CONDITIONAL(HAVE_PLUGIN_MEDCOUPLING, test x$cs_have_plugin_medcoupling = xyes)

cs_py_have_plugin_medcoupling=False
if test x$cs_have_plugin_medcoupling = xyes ; then
  cs_py_have_plugin_medcoupling=True
fi

AC_SUBST(cs_have_medcoupling)
AC_SUBST(cs_have_medcoupling_loader)
AC_SUBST(cs_have_paramedmem)
AC_SUBST(cs_py_have_plugin_medcoupling)
AC_SUBST(MEDCOUPLING_CPPFLAGS)
AC_SUBST(MEDCOUPLING_LDFLAGS)
AC_SUBST(MEDCOUPLING_LIBS)
AC_SUBST(MEDCOUPLINGRUNPATH)

])dnl

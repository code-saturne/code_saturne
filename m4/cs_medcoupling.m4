dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2018 EDF S.A.
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

  if test "x$enable_shared" = "xno" ; then
    AC_MSG_WARN([MEDCoupling support plugin disabled as build is static only.])
    with_medcoupling=no
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
  MEDCOUPLING_LIBS="-lmedcoupling -linterpkernel"

  CPPFLAGS="${CPPFLAGS} ${MEDCOUPLING_CPPFLAGS}"
  LDFLAGS="${MEDCOUPLING_LDFLAGS} ${LDFLAGS}"
  LIBS="${MEDCOUPLING_LIBS} ${LIBS}"

  AC_LANG_PUSH([C++])

  # Check for MEDCoupling library
  #-------------------------------

  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#include <MEDCouplingUMesh.hxx>]],
[[using namespace MEDCoupling;
MEDCouplingUMesh *m = MEDCouplingUMesh::New();]])
                   ],
                   [ AC_DEFINE([HAVE_MEDCOUPLING], 1, [MEDCoupling support])
                     cs_have_medcoupling=yes
                   ],
                   [ AC_MSG_WARN([no MEDCoupling support]) ],
                  )

  # Now check for MEDCoupling MPI support

  if test "x$cs_have_medcoupling" = "xyes" -a "x$cs_have_mpi"; then

    CPPFLAGS="${MPI_CPPFLAGS} ${MEDCOUPLING_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${MEDCOUPLING_LDFLAGS} ${MPI_LDFLAGS} ${LDFLAGS}"
    LIBS="-lparamedmem ${MEDCOUPLING_LIBS} ${MPI_LIBS} ${LIBS}"

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
                   [ AC_MSG_WARN([no ParaMEDMEM support]) ],
                  )

    if test "x$cs_have_paramedmem" = "xyes"; then
      PARAMEDMEM_CPPFLAGS="-I$withval/include/salome"
      PARAMEDMEM_LDFLAGS="-L$withval/lib/salome"
      PARAMEDMEM_LIBS="-lparamedmem ${MEDCOUPLING_LIBS}"
    fi

  fi

  AC_LANG_POP([C++])

  # Report MEDCOUPLING support
  #-------------------

  if test "x$cs_have_medcoupling" = "xno" ; then
    if test "x$with_medcoupling" != "xcheck" ; then
      AC_MSG_FAILURE([MEDCoupling support is requested, but test for MEDCoupling failed!])
    fi
  fi

  if test "x$cs_have_medcoupling" = "xno"; then
    MEDCOUPLING_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_MEDCOUPLING, test x$cs_have_medcoupling = xyes)
AM_CONDITIONAL(HAVE_PARAMEDMEM, test x$cs_have_paramedmem = xyes)

AC_SUBST(cs_have_medcoupling)
AC_SUBST(cs_have_paramedmem)
AC_SUBST(MEDCOUPLING_CPPFLAGS)
AC_SUBST(MEDCOUPLING_LDFLAGS)
AC_SUBST(MEDCOUPLING_LIBS)
AC_SUBST(MEDCOUPLINGRUNPATH)
AC_SUBST(PARAMEDMEM_CPPFLAGS)
AC_SUBST(PARAMEDMEM_LDFLAGS)
AC_SUBST(PARAMEDMEM_LIBS)

])dnl

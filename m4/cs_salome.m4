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

# CS_AC_TEST_SALOME
#------------------
# Macro function grouping the different tests

AC_DEFUN([CS_AC_TEST_SALOME], [

AC_ARG_WITH(salome,
            [AS_HELP_STRING([--with-salome=PATH],
                            [specify prefix directory for SALOME])],
            [if test "x$withval" = "x"; then
               if test "x$ROOT_SALOME" != "x"; then
                 with_salome=$ROOT_SALOME
               elif test "x$INST_ROOT" != "x" -a "x$KERNEL_ROOT_DIR" != "x"; then
                 with_salome=$INST_ROOT
                 ROOT_SALOME=$INST_ROOT
               else
                 with_salome=no
               fi
             else
               if test "x$ROOT_SALOME" = "x"; then
                 ROOT_SALOME=$with_salome
               fi
             fi],
            [if test "x$ROOT_SALOME" = "x"; then
               with_salome=no
             else
               with_salome=$ROOT_SALOME
             fi])

if test "x$with_salome" != "xno" ; then

  # Recommended environment file for salome-platform.org installer builds
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 2 -name salome.sh | tail -1 2>/dev/null)
    if test "x$salome_env" != "x"; then
      SALOMEENVCMD=". $salome_env"
    fi
  fi

  # Environment for EDF or "universal binaries" type build for Salome7.3
  # Note that ROOT_SALOME is required but not exported in all cases.
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 1 -name salome_modules.sh 2>/dev/null)
    if test "x$salome_env" != "x"; then
      salome_pre=$(find $with_salome -maxdepth 1 -name salome_prerequisites.sh 2>/dev/null)
      if test "x$salome_pre" != "x"; then
        SALOMEENVCMD=". $salome_pre; export ROOT_SALOME=$with_salome; . $salome_env"
        SALOMEPRE="$salome_pre"
      fi
    fi
  fi

  # Environment file for EDF builds using YAMM
  # Environment for EDF or "universal binaries" type build for Salome 6.4 to 7.2
  # Note that ROOT_SALOME is required but not exported in all cases.
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 1 -name salome_modules*.sh 2>/dev/null)
    if test "x$salome_env" != "x"; then
      salome_pre=$(find $with_salome -maxdepth 1 -name salome_prerequisites_*appli.sh 2>/dev/null)
      if test "x$salome_pre" != "x"; then
        SALOMEENVCMD=". $salome_pre; export ROOT_SALOME=$with_salome; . $salome_env"
        SALOMEPRE="$salome_pre"
      fi
    fi
  fi

  # Environment file for EDF builds using some versions of YAMM
  if test "x$SALOMEENVCMD" = "x"; then
    if test -f "$with_salome/salome_modules.sh" ; then
      salome_env="$with_salome/salome_modules.sh"
    fi
    if test "x$salome_env" != "x" -a -f "$with_salome/salome_prerequis.sh" ; then
      salome_pre="$with_salome/salome_prerequis.sh"
      if test "x$salome_pre" != "x"; then
        SALOMEENVCMD=". $salome_pre;. $salome_env"
        SALOMEPRE="$salome_pre"
      fi
    fi
  fi

  # Environment file for EDF builds for Salome 6.3
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 1 -name envSalome*.sh 2>/dev/null)
    if test "x$salome_env" != "x"; then
      salome_pre=$(find $with_salome -maxdepth 1 -name prerequis*.sh 2>/dev/null)
      if test "x$salome_pre" != "x"; then
        SALOMEENVCMD=". $salome_pre;. $salome_env"
        SALOMEPRE="$salome_pre"
      fi
    fi

  fi

  if test "x$SALOMERUN" = "x"; then
    if test -f "$ROOT_SALOME/runSalome"; then
      SALOMERUN="$ROOT_SALOME/runSalome"
    fi
  fi

  unset salome_pre
  unset salome_env

  KERNEL_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $KERNEL_ROOT_DIR)
  GUI_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $GUI_ROOT_DIR)
  MEDCOUPLING_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $MEDCOUPLING_ROOT_DIR)
  YACS_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $YACS_ROOT_DIR)

  OMNIIDL=$(eval $SALOMEENVCMD ; which omniidl)

  # Make sure omniidl will work by forcing PYTHONPATH
  # Note that we must ensure we are using SALOME's python here, so use "python" with the sourced
  # environment and PATH ratther than $PYTHON here.

  if test "x$OMNIIDLPYTHONPATH" = "x"; then
   OMNIIDLPYTHONPATH=$(eval $SALOMEENVCMD ; python -B "$srcdir/build-aux/cs_config_test.py" pythonpath_filter _omniidlmodule.so _omniidlmodule.so)
  fi
  if test "x$OMNIIDLLDLIBPATH" = "x"; then
    OMNIIDLLDLIBPATH=$(eval $SALOMEENVCMD ; python -B "$srcdir/build-aux/cs_config_test.py" ld_library_path_filter libpython*)
  fi

fi

AC_SUBST(OMNIIDLPYTHONPATH)
AC_SUBST(OMNIIDLLDLIBPATH)

AC_ARG_VAR([SALOMEENVCMD], [SALOME environment setting commands])
AC_ARG_VAR([SALOMERUN], [SALOME main script (usually runSalome or runAppli)])

AC_SUBST([SALOMEPRE])

CS_AC_TEST_SALOME_KERNEL
CS_AC_TEST_SALOME_GUI
CS_AC_TEST_SALOME_YACS

AS_IF([test $cs_have_salome_kernel = yes -o $cs_have_salome_gui = yes],
      [CS_AC_TEST_OMNIORB
       CS_AC_TEST_CORBA])

AM_CONDITIONAL(HAVE_SALOME,
               test $cs_have_salome_kernel = yes \
                 -a $cs_have_salome_gui = yes \
                 -a $omniORB_ok = yes)

])dnl

# CS_AC_TEST_SALOME_KERNEL
#-------------------------
# modifies or sets cs_have_salome_kernel, SALOME_KERNEL_CPPFLAGS, SALOME_KERNEL_LDFLAGS, and SALOME_KERNEL_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SALOME_KERNEL], [

cs_have_salome_kernel=no

AC_ARG_WITH(salome-kernel,
            [AS_HELP_STRING([--with-salome-kernel=PATH],
                            [specify prefix directory for SALOME KERNEL])],
            [if test "x$withval" = "x"; then
               if test -z "$KERNEL_ROOT_DIR"; then
                 with_salome_kernel=yes
               else
                 with_salome_kernel=$KERNEL_ROOT_DIR
               fi
             fi],
            [if test -z "$KERNEL_ROOT_DIR"; then
               with_salome_kernel=check
             else
               with_salome_kernel=$KERNEL_ROOT_DIR
             fi])

if test "x$with_salome_kernel" != "xno" ; then

  cs_have_salome_kernel=yes

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  if test x"$with_salome_kernel" != xyes -a x"$with_salome_kernel" != xcheck ; then
    SALOME_KERNEL="$with_salome_kernel"
  else
    SALOME_KERNEL="/usr"
  fi

  SALOME_KERNEL_CPPFLAGS="-I$SALOME_KERNEL/include/salome"
  SALOME_KERNEL_IDL="-I$SALOME_KERNEL/idl/salome"
  SALOME_KERNEL_LDFLAGS="-L$SALOME_KERNEL/lib/salome"
  CALCIUM_LIBS="-lCalciumC"

  CPPFLAGS="${CPPFLAGS} ${SALOME_KERNEL_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${SALOME_KERNEL_LDFLAGS}"
  LIBS="${CALCIUM_LIBS} ${LIBS}"

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <calcium.h>]],
  			             [[int iret = cp_fin(0, 0);]])],
                    [cs_have_calcium=yes
                     AC_MSG_RESULT([CALCIUM support])],
                    [cs_have_calcium=no
                     AC_MSG_WARN([no CALCIUM support])
                    ])

  if test "x$cs_have_salome_kernel" = "xno"; then
    CALCIUM_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(cs_have_salome_kernel)
AC_SUBST(SALOME_KERNEL)
AC_SUBST(SALOME_KERNEL_CPPFLAGS)
AC_SUBST(SALOME_KERNEL_IDL)
AC_SUBST(SALOME_KERNEL_LDFLAGS)
AC_SUBST(SALOME_KERNEL_LIBS)

AC_SUBST(cs_have_calcium)
AC_SUBST(CALCIUM_LIBS)

AM_CONDITIONAL(HAVE_SALOME_KERNEL, test x$cs_have_salome_kernel = xyes)
AM_CONDITIONAL(HAVE_CALCIUM, test x$cs_have_calcium = xyes)

])dnl


# CS_AC_TEST_SALOME_GUI
#----------------------
# modifies or sets cs_have_salome_gui, SALOME_GUI_CPPFLAGS, SALOME_GUI_LDFLAGS, and SALOME_GUI_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SALOME_GUI], [

cs_have_salome_gui=no

AC_ARG_WITH(salome-gui,
            [AS_HELP_STRING([--with-salome-gui=PATH],
                            [specify prefix directory for SALOME GUI])],
            [if test "x$withval" = "x"; then
               if test -z "$GUI_ROOT_DIR"; then
                 with_salome_gui=yes
               else
                 with_salome_gui=$GUI_ROOT_DIR
               fi
             fi],
            [if test -z "$GUI_ROOT_DIR"; then
               with_salome_gui=check
             else
               with_salome_gui=$GUI_ROOT_DIR
             fi])

#if test "x$with_salome_gui" != "xno" ; then
# We should add a couple of tests here...
#fi

if test x"$with_salome_gui" != xno ; then

  cs_have_salome_gui=yes

  if test x"$with_salome_gui" != xyes -a x"$with_salome_gui" != xcheck ; then
    SALOME_GUI="$with_salome_gui"
  else
    SALOME_GUI="/usr"
  fi

  SALOME_GUI_CPPFLAGS="-I$SALOME_GUI/include/salome"
  SALOME_GUI_IDL="-I$SALOME_GUI/idl/salome"
  SALOME_GUI_LDFLAGS="-L$SALOME_GUI/lib/salome"
  SALOME_GUI_LIBS=

else
  cs_have_salome_gui=no
fi

AC_SUBST(cs_have_salome_gui)
AC_SUBST(SALOME_GUI)
AC_SUBST(SALOME_GUI_CPPFLAGS)
AC_SUBST(SALOME_GUI_IDL)
AC_SUBST(SALOME_GUI_LDFLAGS)
AC_SUBST(SALOME_GUI_LIBS)

AM_CONDITIONAL(HAVE_SALOME_GUI, test x$cs_have_salome_gui = xyes)

])dnl


# CS_AC_TEST_SALOME_YACS
#----------------------
# modifies or sets cs_have_salome_yacs, SALOME_YACS_CPPFLAGS, SALOME_YACS_LDFLAGS, and SALOME_YACS_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SALOME_YACS], [

cs_have_salome_yacs=no

AC_ARG_WITH(salome-yacs,
            [AS_HELP_STRING([--with-salome-yacs=PATH],
                            [specify prefix directory for SALOME YACS])],
            [if test "x$withval" = "x"; then
               if test -z "$YACS_ROOT_DIR"; then
                 with_salome_yacs=yes
               else
                 with_salome_yacs=$YACS_ROOT_DIR
               fi
             fi],
            [if test -z "$YACS_ROOT_DIR"; then
               with_salome_yacs=check
             else
               with_salome_yacs=$YACS_ROOT_DIR
             fi])

#if test "x$with_salome_yacs" != "xno" ; then
# We should add a couple of tests here...
#fi

if test x"$with_salome_yacs" != xno ; then

  cs_have_salome_yacs=yes

  if test x"$with_salome_yacs" != xyes -a x"$with_salome_yacs" != xcheck ; then
    SALOME_YACS="$with_salome_yacs"
  else
    SALOME_YACS="/usr"
  fi

  SALOME_YACS_CPPFLAGS="-I$SALOME_YACS/include/salome"
  SALOME_YACS_IDL="-I$SALOME_YACS/idl/salome"
  SALOME_YACS_LDFLAGS="-L$SALOME_YACS/lib/salome"
  SALOME_YACS_LIBS=

else
  cs_have_salome_yacs=no
fi

AC_SUBST(cs_have_salome_yacs)
AC_SUBST(SALOME_YACS)
AC_SUBST(SALOME_YACS_CPPFLAGS)
AC_SUBST(SALOME_YACS_IDL)
AC_SUBST(SALOME_YACS_LDFLAGS)
AC_SUBST(SALOME_YACS_LIBS)

AM_CONDITIONAL(HAVE_SALOME_YACS, test x$cs_have_salome_yacs = xyes)

])dnl


# CS_AC_ENABLE_PTHREAD
#---------------------
# Modify CFLAGS, CXXFLAGS and LIBS for compiling pthread-based programs.
# Use acx_pthread.m4 from GNU Autoconf Macro Archive

AC_DEFUN([CS_AC_ENABLE_PTHREADS],[
AC_REQUIRE([ACX_PTHREAD])

if test x"$cs_enable_pthreads_done" != xyes; then
  if test x"$acx_pthread_ok" = xyes; then
    CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
    LIBS="$LIBS $PTHREAD_LIBS"
    cs_threads_ok=yes
  else
    cs_threads_ok=no
  fi
  cs_enable_pthreads_done=yes
fi
])dnl
dnl


# CS_AC_TEST_CORBA
#-----------------
# Set CORBA-related variables according to omniORB detection.

AC_DEFUN([CS_AC_TEST_CORBA],[

if test x"$DEFAULT_ORB" = x"omniORB"
then

  #  Contient le nom de l'ORB
  ORB=omniorb

  AC_MSG_RESULT(default orb: omniORB)
  IDL=$OMNIORB_IDL
  AC_SUBST(IDL)

  CORBA_ROOT=$OMNIORB_ROOT
  CORBA_INCLUDES=$OMNIORB_INCLUDES
  CORBA_CXXFLAGS=$OMNIORB_CXXFLAGS
  CORBA_LIBS=$OMNIORB_LIBS
  IDLCXXFLAGS=$OMNIORB_IDLCXXFLAGS
  IDLPYFLAGS=$OMNIORB_IDLPYFLAGS

  AC_SUBST(CORBA_ROOT)
  AC_SUBST(CORBA_INCLUDES)
  AC_SUBST(CORBA_CXXFLAGS)
  AC_SUBST(CORBA_LIBS)
  AC_SUBST(IDLCXXFLAGS)
  AC_SUBST(IDLPYFLAGS)

  IDL_CLN_H=$OMNIORB_IDL_CLN_H
  IDL_CLN_CXX=$OMNIORB_IDL_CLN_CXX
  IDL_CLN_OBJ=$OMNIORB_IDL_CLN_OBJ

  AC_SUBST(IDL_CLN_H)
  AC_SUBST(IDL_CLN_CXX)
  AC_SUBST(IDL_CLN_OBJ)

  IDL_SRV_H=$OMNIORB_IDL_SRV_H
  IDL_SRV_CXX=$OMNIORB_IDL_SRV_CXX
  IDL_SRV_OBJ=$OMNIORB_IDL_SRV_OBJ

  AC_SUBST(IDL_SRV_H)
  AC_SUBST(IDL_SRV_CXX)
  AC_SUBST(IDL_SRV_OBJ)

else
  AC_MSG_RESULT($DEFAULT_ORB unknown orb)
fi

])dnl
dnl

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
               if test -z "$MEDCOUPLING_ROOT_DIR"; then
                 with_medcoupling=yes
               else
                 with_medcoupling=$MEDCOUPLING_ROOT_DIR
               fi
             fi],
            [if test -z "$MEDCOUPLING_ROOT_DIR"; then
               with_medcoupling=check
             else
               with_medcoupling=$MEDCOUPLING_ROOT_DIR
             fi])

if test "x$with_medcoupling" != "xno" -a "x$enable_shared" = "xno" ; then
  AC_MSG_WARN([MEDCoupling support plugin disabled as build is static only.])
  with_medcoupling=no
fi

if test "x$with_medcoupling" != "xno" ; then

  if test x"$with_medcoupling" != xyes -a x"$with_medcoupling" != xcheck ; then
    MEDCOUPLING="$with_medcoupling"
  else
    MEDCOUPLING="/usr"
  fi

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  MEDCOUPLING_CPPFLAGS="-I$MEDCOUPLING/include"
  MEDCOUPLING_LDFLAGS="-L$MEDCOUPLING/lib"
  # Add the libdir to the runpath as libtool does not do this for modules
  MEDCOUPLINGRUNPATH="-R$MEDCOUPLING/lib"
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
    else
      AC_MSG_WARN([no MEDCoupling support])
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

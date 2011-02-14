dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2011 EDF S.A., France
dnl
dnl   The Code_Saturne Kernel is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne Kernel is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne Preprocessor; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_TEST_SALOME
#------------------
# Macro function grouping the different tests

AC_DEFUN([CS_AC_TEST_SALOME], [

CS_AC_TEST_SALOME_KERNEL
CS_AC_TEST_SALOME_GUI


AS_IF([test $cs_have_salome_kernel = yes -o $cs_have_salome_gui = yes],
      [CS_AC_TEST_OMNIORB
       CS_AC_TEST_CORBA
       CS_AC_TEST_BOOST
       AC_SUBST(ROOT_SALOME)
       AC_SUBST(KERNEL_ROOT_DIR)
       AC_SUBST(GUI_ROOT_DIR)
       AC_SUBST(YACS_ROOT_DIR)
       AC_SUBST(SALOME_VERSION, [`basename $KERNEL_ROOT_DIR|sed -e 's/KERNEL_//' -`])])

CS_AC_TEST_SPHINX

])dnl

# CS_AC_TEST_SALOME_KERNEL
#-------------------------
# modifies or sets cs_have_salome_kernel, SALOME_KERNEL_CPPFLAGS, SALOME_KERNEL_LDFLAGS, and SALOME_KERNEL_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SALOME_KERNEL], [

cs_have_salome_kernel=no

AC_ARG_WITH(salome-kernel,
            [AS_HELP_STRING([--with-salome-kernel=PATH],
                            [specify prefix directory for SALOME kernel])],
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

AC_ARG_WITH(salome-kernel-include,
            [AS_HELP_STRING([--with-salome-kernel-include=PATH],
                            [specify directory for SALOME kernel include files])],
            [if test "x$with_salome_kernel" = "xcheck"; then
               with_salome_kernel=yes
             fi
             SALOME_KERNEL_CPPFLAGS="-I$with_salome_kernel_include"],
            [if test "x$with_salome_kernel" != "xno" ; then
               if test "x$with_salome_kernel" != "xyes" \
	               -a "x$with_salome_kernel" != "xcheck"; then
                 SALOME_KERNEL_CPPFLAGS="-I$with_salome_kernel/include/salome"
               else
                 SALOME_KERNEL_CPPFLAGS="-I/usr/include/salome"
               fi
             fi])

AC_ARG_WITH(salome-kernel-idl,
            [AS_HELP_STRING([--with-salome-kernel-idl=PATH],
                            [specify directory for SALOME_KERNEL IDL files])],
            [if test "x$with_salome_kernel" = "xcheck"; then
               with_salome_kernel=yes
             fi
             SALOME_KERNEL_IDL="-L$with_salome_kernel_idl"],
            [if test "x$with_salome_kernel" != "xno" ; then
               if test "x$with_salome_kernel" != "xyes" \
	               -a "x$with_salome_kernel" != "xcheck"; then
                 SALOME_KERNEL_IDL="-I$with_salome_kernel/idl/salome"
               else
                 SALOME_KERNEL_IDL="-I/usr/idl/salome"
               fi
             fi])

AC_ARG_WITH(salome-kernel-lib,
            [AS_HELP_STRING([--with-salome-kernel-lib=PATH],
                            [specify directory for SALOME_KERNEL library])],
            [if test "x$with_salome_kernel" = "xcheck"; then
               with_salome_kernel=yes
             fi
             SALOME_KERNEL_LDFLAGS="-L$with_salome_kernel_lib"],
            [if test "x$with_salome_kernel" != "xno" ; then
               if test "x$with_salome_kernel" != "xyes" \
	               -a "x$with_salome_kernel" != "xcheck"; then
                 SALOME_KERNEL_LDFLAGS="-L$with_salome_kernel/lib/salome"
               else
                 SALOME_KERNEL_LDFLAGS="-L/usr/lib/salome"
               fi
             fi])

if test "x$with_salome_kernel" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  SALOME_KERNEL_LIBS="-lCalciumC"
  
  CPPFLAGS="${CPPFLAGS} ${SALOME_KERNEL_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${SALOME_KERNEL_LDFLAGS}"
  LIBS="${LIBS} ${SALOME_KERNEL_LIBS}"

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <calcium.h>]],
  			             [[int iret = cp_fin(0, 0);]])],
                    [cs_have_salome_kernel=yes
                     AC_MSG_RESULT([compatible SALOME kernel found])],
                    [cs_have_salome_kernel=no
                     if test "x$with_salome_kernel" != "xcheck" ; then
                       AC_MSG_FAILURE([SALOME support is requested, but test for SALOME failed!])
                     else
                       AC_MSG_WARN([no SALOME support])
                     fi
                    ])

  if test "x$cs_have_salome_kernel" = "xno"; then
    SALOME_KERNEL_LIBS=""
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AC_SUBST(cs_have_salome_kernel)
AC_SUBST(SALOME_KERNEL_CPPFLAGS)
AC_SUBST(SALOME_KERNEL_IDL)
AC_SUBST(SALOME_KERNEL_LDFLAGS)
AC_SUBST(SALOME_KERNEL_LIBS)

AM_CONDITIONAL(HAVE_SALOME_KERNEL, test x$cs_have_salome_kernel = xyes)

])dnl


# CS_AC_TEST_SALOME_GUI
#----------------------
# modifies or sets cs_have_salome_gui, SALOME_GUI_CPPFLAGS, SALOME_GUI_LDFLAGS, and SALOME_GUI_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_SALOME_GUI], [

cs_have_salome_gui=no

AC_ARG_WITH(salome-gui,
            [AS_HELP_STRING([--with-salome-gui=PATH],
                            [specify prefix directory for SALOME gui])],
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

AC_ARG_WITH(salome-gui-include,
            [AS_HELP_STRING([--with-salome-gui-include=PATH],
                            [specify directory for SALOME gui include files])],
            [if test "x$with_salome_gui" = "xcheck"; then
               with_salome_gui=yes
             fi
             SALOME_GUI_CPPFLAGS="-I$with_salome_gui_include"],
            [if test "x$with_salome_gui" != "xno" ; then
               if test "x$with_salome_gui" != "xyes" \
	               -a "x$with_salome_gui" != "xcheck"; then
                 SALOME_GUI_CPPFLAGS="-I$with_salome_gui/include/salome"
               else
                 SALOME_GUI_CPPFLAGS="-I/usr/include/salome"
               fi
             fi])

AC_ARG_WITH(salome-gui-idl,
            [AS_HELP_STRING([--with-salome-gui-idl=PATH],
                            [specify directory for SALOME_GUI IDL files])],
            [if test "x$with_salome_gui" = "xcheck"; then
               with_salome_gui=yes
             fi
             SALOME_GUI_IDL="-L$with_salome_gui_idl"],
            [if test "x$with_salome_gui" != "xno" ; then
               if test "x$with_salome_gui" != "xyes" \
	               -a "x$with_salome_gui" != "xcheck"; then
                 SALOME_GUI_IDL="-I$with_salome_gui/idl/salome"
               else
                 SALOME_GUI_IDL="-I/usr/idl/salome"
               fi
             fi])

AC_ARG_WITH(salome-gui-lib,
            [AS_HELP_STRING([--with-salome-gui-lib=PATH],
                            [specify directory for SALOME_GUI library])],
            [if test "x$with_salome_gui" = "xcheck"; then
               with_salome_gui=yes
             fi
             SALOME_GUI_LDFLAGS="-L$with_salome_gui_lib"],
            [if test "x$with_salome_gui" != "xno" ; then
               if test "x$with_salome_gui" != "xyes" \
	               -a "x$with_salome_gui" != "xcheck"; then
                 SALOME_GUI_LDFLAGS="-L$with_salome_gui/lib/salome"
               else
                 SALOME_GUI_LDFLAGS="-L/usr/lib/salome"
               fi
             fi])

#if test "x$with_salome_gui" != "xno" ; then
# We should add a couple of tests here...
#fi

cs_have_salome_gui=yes

AC_SUBST(cs_have_salome_gui)
AC_SUBST(SALOME_GUI_CPPFLAGS)
AC_SUBST(SALOME_GUI_IDL)
AC_SUBST(SALOME_GUI_LDFLAGS)
AC_SUBST(SALOME_GUI_LIBS)

AM_CONDITIONAL(HAVE_SALOME_GUI, test x$cs_have_salome_gui = xyes)

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


# CS_AC_TEST_BOOST
#-----------------
# Boost location for a couple of headers needed by Calcium.hxx

AC_DEFUN([CS_AC_TEST_BOOST],[
AC_REQUIRE([CS_AC_ENABLE_PTHREADS])dnl

AC_LANG_SAVE
AC_LANG_CPLUSPLUS

BOOST_CPPFLAGS=""

AC_CHECKING(for BOOST headers location)

AC_ARG_WITH(boost-include,
            [AS_HELP_STRING([--with-boost-include=PATH],
                            [specify directory for BOOST include files])],
            [BOOST_CPPFLAGS="-I$with_boost_include"],
            [if ! test -z "$BOOSTDIR"; then
               BOOST_CPPFLAGS="-I${BOOSTDIR}/include"
             fi])

CPPFLAGS_old="${CPPFLAGS}"

if test "x${BOOSTDIR}" != "x" ; then
  BOOST_CPPFLAGS="-I${BOOSTDIR}/include"
fi

AC_TRY_COMPILE([#include <boost/shared_ptr.hpp>],
               [boost::shared_ptr<int>(new int)],
                boost_headers_ok=yes,
                boost_headers_ok=no)

CPPFLAGS="${CPPFLAGS_old}"

AC_MSG_RESULT(for boost headers: $boost_headers_ok)

AC_SUBST(BOOST_CPPFLAGS)

AC_LANG_RESTORE

])dnl


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

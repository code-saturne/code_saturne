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

dnl  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
dnl
dnl  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
dnl  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS

AC_DEFUN([CS_AC_TEST_OMNIORB],[
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_CXX])dnl
AC_REQUIRE([AC_PROG_CPP])dnl
AC_REQUIRE([AC_PROG_CXXCPP])dnl

AC_CHECKING(for omniORB)
omniORB_ok=yes

AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_ARG_VAR([OMNIIDL], [the omniORB IDL compiler])
AC_PATH_PROG([OMNIIDL], [omniidl])

AC_ARG_VAR([OMNIIDLPYBE], [Python backend for omniidl])

if test "x$OMNIIDL" = "x"
then
  omniORB_ok=no
  AC_MSG_RESULT(omniORB binaries not in PATH variable)
else
  omniORB_ok=yes
fi

if  test "x$omniORB_ok" = "xyes"
then
  OMNIORB_BIN=`echo ${OMNIIDL} | sed -e "s,[[^/]]*$,,;s,/$,,;s,^$,.,"`
  OMNIORB_ROOT=${OMNIORB_BIN}
  # one-level up
  OMNIORB_ROOT=`echo ${OMNIORB_ROOT}  | sed -e "s,[[^/]]*$,,;s,/$,,;s,^$,.,"`
  #
  #
  if test -d $OMNIORB_ROOT/include ; then
    # if $OMNIORB_ROOT/include exists, there are a lot of chance that
    # this is omniORB4.x installed via configure && make && make install
    OMNIORB_LIB=`echo ${OMNIORB_BIN} | sed -e "s,bin\$,lib,"`
    OMNIORB_VERSION=4
  else
    # omniORB has been installed old way
    OMNIORB_LIB=`echo ${OMNIORB_BIN} | sed -e "s,bin/,lib/,"`
    # one-level up again
    OMNIORB_ROOT=`echo ${OMNIORB_ROOT}  | sed -e "s,[[^/]]*$,,;s,/$,,;s,^$,.,"`
    if test -d $OMNIORB_ROOT/include/omniORB4 ; then
      OMNIORB_VERSION=4
    else
      OMNIORB_VERSION=3
    fi
  fi
  AC_SUBST(OMNIORB_ROOT)

  OMNIORB_INCLUDES="-I$OMNIORB_ROOT/include -I$OMNIORB_ROOT/include/omniORB${OMNIORB_VERSION} -I$OMNIORB_ROOT/include/COS"
  AC_SUBST(OMNIORB_INCLUDES)

  CS_AC_ENABLE_PTHREADS

  OMNIORB_CXXFLAGS="-DOMNIORB_VERSION=$OMNIORB_VERSION"
  case $build_cpu in
    sparc*)
      # AC_DEFINE(__sparc__)
      OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -D__sparc__"
      ;;
   *86*)
      # AC_DEFINE(__x86__)
      OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -D__x86__"
      ;;
  esac
  case $build_os in
    solaris*)
      # AC_DEFINE(__sunos__)
      __OSVERSION__=5
      AC_DEFINE([__OSVERSION__], [5], [OmniORB OS version])
      OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -D__sunos__"
      ;;
    mingw*)
      # AC_DEFINE(__WIN32__)
      __OSVERSION__=4
      AC_DEFINE([__OSVERSION__], [4], [OmniORB OS version])
      OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -D__WIN32__"
      ;;
    linux*)
      # AC_DEFINE(__linux__)
      __OSVERSION__=2
      AC_DEFINE([__OSVERSION__], [2], [OmniORB OS version])
      OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -D__linux__"
      ;;
  esac
  AC_SUBST(OMNIORB_CXXFLAGS)

  CPPFLAGS_old=$CPPFLAGS
  CPPFLAGS="$CPPFLAGS $OMNIORB_CXXFLAGS $OMNIORB_INCLUDES"

  AC_CHECK_HEADER(CORBA.h,omniORB_ok="yes",omniORB_ok="no")

  CPPFLAGS=$CPPFLAGS_old

fi

if test "x$omniORB_ok" = "xyes"
then
  if test "x$OMNIORB_LIB" = "x/usr/lib"
  then
    OMNIORB_LDFLAGS=""
    OMNIORB_RFLAGS=""
  else
    OMNIORB_LDFLAGS="-L$OMNIORB_LIB"
    OMNIORB_RFLAGS="-R$OMNIORB_LIB"
  fi

  if test $host_os = mingw32 ; then
    OMNILIBSUFFIX="_rt"
  else
    OMNILIBSUFFIX=""
  fi
  LIBS_old=$LIBS
  LIBS="$OMNIORB_LDFLAGS -lomnithread$OMNILIBSUFFIX $LIBS"

  CXXFLAGS_old=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $OMNIORB_CXXFLAGS $OMNIORB_INCLUDES"

  AC_MSG_CHECKING(whether we can link with omnithreads)
  AC_CACHE_VAL(salome_cv_lib_omnithreads,[
    AC_TRY_LINK(
#include <omnithread.h>
,   omni_mutex my_mutex,
    eval "salome_cv_lib_omnithreads=yes",eval "salome_cv_lib_omnithreads=no")
  ])

  omniORB_ok="$salome_cv_lib_omnithreads"
  if  test "x$omniORB_ok" = "xno"
  then
    AC_MSG_RESULT(omnithreads not found)
  else
    AC_MSG_RESULT(yes)
  fi

  LIBS=$LIBS_old
  CXXFLAGS=$CXXFLAGS_old
fi


dnl omniORB_ok=yes
if test "x$omniORB_ok" = "xyes"
then

  AC_CHECK_LIB(socket,socket, LIBS="-lsocket $LIBS",,)
  AC_CHECK_LIB(nsl,gethostbyname, LIBS="-lnsl $LIBS",,)

  LIBS_old=$LIBS
  OMNIORB_LIBS="$OMNIORB_LDFLAGS"
  OMNIORB_LIBS="$OMNIORB_LIBS -lomniORB${OMNIORB_VERSION}$OMNILIBSUFFIX"
  OMNIORB_LIBS="$OMNIORB_LIBS -lomniDynamic${OMNIORB_VERSION}$OMNILIBSUFFIX"
  OMNIORB_LIBS="$OMNIORB_LIBS -lCOS${OMNIORB_VERSION}$OMNILIBSUFFIX"
  OMNIORB_LIBS="$OMNIORB_LIBS -lCOSDynamic${OMNIORB_VERSION}$OMNILIBSUFFIX"
  OMNIORB_LIBS="$OMNIORB_LIBS -lomnithread$OMNILIBSUFFIX"
  OMNIORB_LIBS="$OMNIORB_LIBS ${OMNIORB_RFLAGS}"
  if test $OMNIORB_VERSION = 3 ; then
    OMNIORB_LIBS="$OMNIORB_LIBS -ltcpwrapGK"
  fi
  AC_SUBST(OMNIORB_LIBS)

  LIBS_old=$LIBS
  LIBS="$OMNIORB_LIBS $LIBS"

  CXXFLAGS_old=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $OMNIORB_CXXFLAGS $OMNIORB_INCLUDES"

  AC_MSG_CHECKING(whether we can link with omniORB)
  AC_CACHE_VAL(salome_cv_lib_omniorb,[
    AC_TRY_LINK(
#include <CORBA.h>
,   CORBA::ORB_var orb,
    eval "salome_cv_lib_omniorb3=yes",eval "salome_cv_lib_omniorb3=no")
  ])
  omniORB_ok="$salome_cv_lib_omniorb3"

  omniORB_ok=yes
  if test "x$omniORB_ok" = "xno"
  then
    AC_MSG_RESULT(omniORB library linking failed)
    omniORB_ok=no
  else
    AC_MSG_RESULT(yes)
  fi
  LIBS="$LIBS_old"
  CXXFLAGS=$CXXFLAGS_old
fi


if test "x$omniORB_ok" = "xyes"
then

  OMNIORB_IDLCXXFLAGS="-Wba -nf -I${OMNIORB_ROOT}/idl"
  OMNIORB_IDLPYFLAGS_1='-bpython -nf '
  OMNIORB_IDLPYFLAGS_2=" -I${OMNIORB_ROOT}/idl"
  OMNIORB_IDLPYFLAGS=${OMNIORB_IDLPYFLAGS_1}${OMNIORB_IDLPYFLAGS_2}

  AC_SUBST(OMNIORB_IDLCXXFLAGS)
  AC_SUBST(OMNIORB_IDLPYFLAGS)

  OMNIORB_IDL_CLN_H=.hh
  OMNIORB_IDL_CLN_CXX=SK.cc
  OMNIORB_IDL_CLN_OBJ=SK.o
  AC_SUBST(OMNIORB_IDL_CLN_H)
  AC_SUBST(OMNIORB_IDL_CLN_CXX)
  AC_SUBST(OMNIORB_IDL_CLN_OBJ)

  OMNIORB_IDL_SRV_H=.hh
  OMNIORB_IDL_SRV_CXX=SK.cc
  OMNIORB_IDL_SRV_OBJ=SK.o
  AC_SUBST(OMNIORB_IDL_SRV_H)
  AC_SUBST(OMNIORB_IDL_SRV_CXX)
  AC_SUBST(OMNIORB_IDL_SRV_OBJ)

  OMNIORB_IDL_TIE_H=
  OMNIORB_IDL_TIE_CXX=
  AC_SUBST(OMNIORB_IDL_TIE_H)
  AC_SUBST(OMNIORB_IDL_TIE_CXX)

  AC_DEFINE([OMNIORB], [], [Description])

  CORBA_HAVE_POA=1
  AC_DEFINE([CORBA_HAVE_POA], [], [Description])

  CORBA_ORB_INIT_HAVE_3_ARGS=1
  AC_DEFINE([CORBA_ORB_INIT_HAVE_3_ARGS], [], [Description])
  CORBA_ORB_INIT_THIRD_ARG='"omniORB"'
  AC_DEFINE([CORBA_ORB_INIT_THIRD_ARG], [], [omniORB])

fi

omniORBpy_ok=no
if  test "x$omniORB_ok" = "xyes"
then
  AC_MSG_CHECKING(omniORBpy)
  $PYTHON -c "import omniORB" &> /dev/null
  if test $? = 0 ; then
    AC_MSG_RESULT(yes)
    omniORBpy_ok=yes
  else
    AC_MSG_RESULT(no, check your installation of omniORBpy)
    omniORBpy_ok=no
  fi
fi

dnl AC_LANG_RESTORE

AC_MSG_RESULT(for omniORBpy: $omniORBpy_ok)
AC_MSG_RESULT(for omniORB: $omniORB_ok)

# Save cache
AC_CACHE_SAVE

if test "x$omniORB_ok" = "xyes"
then

  CXXFLAGS_old=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $OMNIORB_CXXFLAGS $OMNIORB_INCLUDES"
  LIBS_old=$LIBS
  LIBS="$OMNIORB_LDFLAGS $OMNIORB_LIBS $LIBS"
  AC_MSG_CHECKING(whether we have double and CORBA::Double compatibility)

  # Note: in cross-compilation mode, run-time checks are replaced
  # by compile-time checks; the latter are based on the fact that
  # the case statement does not allow twice the same label
  # (autoconf uses a similar test for some of its macros)

  AC_RUN_IFELSE([AC_LANG_PROGRAM([[#include <stdlib.h>
#include <CORBA.h>]],
                [[CORBA::Double *a=new CORBA::Double(2.5);
  double c=2.5;
  double *b;
  b=(double *)a;

  if( (c==*b) && (sizeof(double)==sizeof(CORBA::Double)) ){
    delete a;
    exit(0);
  }
  else{
    delete a;
    exit(1);
  }]]) ],
  DOUBLECOMP="yes",
  DOUBLECOMP="no",
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdlib.h>
#include <CORBA.h>]],
                     [[switch (0) case 0: case (sizeof(double)==sizeof(CORBA::Double)):;]])],
                     [ DOUBLECOMP="yes" ],
                     [ DOUBLECOMP="no"])]
  )
  if test "$DOUBLECOMP" = yes; then
    OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -DCOMP_CORBA_DOUBLE"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi
  AC_MSG_CHECKING(whether we have int and CORBA::Long compatibility)
  AC_RUN_IFELSE([AC_LANG_PROGRAM([[#include <stdlib.h>
#include <CORBA.h>]],
                [[CORBA::Long *a=new CORBA::Long(2);
  int c=2;
  int *b;
  b=(int *)a;

  if( (c==*b) && (sizeof(int)==sizeof(CORBA::Long)) )
    exit(0);
  else
    exit(1);]]) ],
  LONGCOMP="yes",
  LONGCOMP="no",
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdlib.h>
#include <CORBA.h>]],
  [[switch (0) case 0: case (sizeof(int)==sizeof(CORBA::Long)):;]])],
                     [ LONGCOMP="yes" ],
                     [ LONGCOMP="no"])]
  )
  if test "$LONGCOMP" = yes; then
    OMNIORB_CXXFLAGS="$OMNIORB_CXXFLAGS -DCOMP_CORBA_LONG"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
  fi

  CXXFLAGS=$CXXFLAGS_old
  LIBS=$LIBS_old
fi

AC_LANG_RESTORE

AC_SUBST(OMNIORB_CXXFLAGS)

])dnl
dnl

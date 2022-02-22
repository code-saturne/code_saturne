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

# CS_AC_TEST_CATALYST
#--------------------
# modifies or sets cs_have_catalyst, CATALYST_CPPFLAGS, CATALYST_LDFLAGS,
# and CATALYST_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_CATALYST], [

cs_have_catalyst=no
cs_have_plugin_catalyst=yes

cs_catalyst_version=""

# Configure options
#------------------

AC_ARG_WITH(catalyst,
            [AS_HELP_STRING([--with-catalyst=PATH],
                            [specify prefix directory for Catalyst])],
            [if test "x$withval" = "x"; then
               with_catalyst=no
             elif test "x$withval" = "xsalome"; then
               if test -z "$CATALYST_ROOT_DIR"; then
                 AC_MSG_FAILURE([no SALOME  path information for Catalyst
(CATALYST_ROOT_DIR environment variable needed by --with-catalyst=salome)!])
               else
                 with_catalyst=$CATALYST_ROOT_DIR
               fi
             fi],
            [with_catalyst=no])

AC_ARG_WITH(catalyst-version,
            [AS_HELP_STRING([--with-catalyst-version=x.y],
                            [Specify Catalyst version (if <= 5.6)])],
            [if test "x$withval" != "x"; then
               cs_catalyst_version="$withval"
             fi],
            [cs_catalyst_version=""])

AC_ARG_ENABLE(catalyst-as-plugin,
  [AS_HELP_STRING([--disable-catalyst-as-plugin], [do not use Catalyst as plugin])],
  [
    case "${enableval}" in
      yes) cs_have_plugin_catalyst=yes ;;
      no)  cs_have_plugin_catalyst=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-catalyst-as-plugin]) ;;
    esac
  ],
  [ cs_have_plugin_catalyst=yes ]
)

if test x$cs_have_dlloader = xno -o x$enable_shared = xno ; then
  cs_have_plugin_catalyst=no
fi

# Check for CMake first
#----------------------

if test "x$with_catalyst" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CATALYST_LIBS=""

  CPPFLAGS="${CPPFLAGS} ${CATALYST_CPPFLAGS}"
  LDFLAGS="${CATALYST_LDFLAGS} ${LDFLAGS}"
  LIBS="${CATALYST_LIBS} ${LIBS}"

  # Test for CMake

  AC_ARG_VAR([CMAKE], [CMake cross platform make])

  if test "x$CMAKE" = "x" ; then
    AC_PATH_PROGS([CMAKE], [cmake])
  fi

  if test "$CMAKE" = : ; then
    AC_MSG_FAILURE([cannot find CMake, Catalyst plugin cannot be installed])
    with_catalyst=no
  fi

fi

if test "x$with_catalyst" != "xno" ; then

  # Check for a Catalyst library (using different CMakeLists variants depending on version)
  #-----------------------------

  detection_variant=catalyst

  if test "x$cs_catalyst_version" != "x"; then
    case "$cs_catalyst_version" in
    5.6|5.5|5.4)
      detection_variant=catalyst-5.6
      ;;
    esac
  fi

  AC_MSG_CHECKING([ParaView/Catalyst])

  cs_prv_dir=`pwd`
  rm -rf "$cs_prv_dir"/catalyst_test
  rm -rf "$cs_prv_dir"/CMakeFiles
  rm -rf "$cs_prv_dir"/CMakeCache.txt

  cs_abs_srcdir=`cd $srcdir && pwd`

  # Work around some detection issues on some systems
  catalyst_cmake_options=""
  if test "x$TBB_INCLUDE_DIR" != "x" ; then
    catalyst_cmake_options="$catalyst_cmake_options -DTBB_INCLUDE_DIR=${TBB_INCLUDE_DIR}"
  fi

  mkdir catalyst_test && cd catalyst_test
  "$CMAKE" -DCMAKE_PREFIX_PATH="$with_catalyst" "$cs_abs_srcdir/build-aux/$detection_variant" $catalyst_cmake_options >&5

  if test $? = 0 ; then

    CATALYST_CPPFLAGS=`grep CXX_DEFINES CMakeFiles/CoProcessingTest.dir/flags.make | sed -e 's/^@<:@^=@:>@*=@<:@\ @:>@*//'`
    CATALYST_INCLUDES=`grep CXX_INCLUDES CMakeFiles/CoProcessingTest.dir/flags.make | sed -e 's/^@<:@^=@:>@*=@<:@\ @:>@*//'`
    CATALYST_CPPFLAGS="${CATALYST_CPPFLAGS} ${CATALYST_INCLUDES}"
    CATALYST_CXXFLAGS=`grep CXX_FLAGS CMakeFiles/CoProcessingTest.dir/flags.make | sed -e 's/^@<:@^=@:>@*=@<:@\ @:>@*//'`
    cs_link_opts=`sed -e 's/^.*CoProcessingTest//' CMakeFiles/CoProcessingTest.dir/link.txt`
    for a in $cs_link_opts; do
      case "{$a}" in
      -l*) CATALYST_LIBS="${CATALYST_LIBS} ${a}"
           ;;
      *)   if test -f "$a" ; then
             CATALYST_LIBS="${CATALYST_LIBS} -Wl,${a}"
           else
             CATALYST_LDFLAGS="${CATALYST_LDFLAGS} ${a}"
           fi
           ;;
      esac
    done
    unset cs_link_opts

  fi

  cd "$cs_prv_dir"

  if test "x$CATALYST_LIBS" != "x" ; then
    cs_have_catalyst=yes;
  fi

  AC_MSG_RESULT($cs_have_catalyst)

  # Report Catalyst support
  #------------------------

  if test "x$cs_have_catalyst" = "xyes" ; then
    AC_DEFINE([HAVE_CATALYST], 1, [Catalyst co-processing support])
    if test x$cs_have_plugin_catalyst = xyes ; then
      AC_DEFINE([HAVE_PLUGIN_CATALYST], 1, [Catalyst co-processing support as plugin])
    fi
  elif test "x$cs_have_catalyst" = "xno" ; then
    if test "x$with_catalyst" != "xcheck" ; then
      AC_MSG_FAILURE([Catalyst co-processing support requested, but test for Catalyst failed!])
    else
      AC_MSG_WARN([no Catalyst co-processing support])
    fi
  else
    AC_MSG_WARN([no Catalyst test as expected])
  fi

  rm -rf "$cs_prv_dir"/catalyst_test
  rm -rf "$cs_prv_dir"/CMakeFiles
  rm -rf "$cs_prv_dir"/CMakeCache.txt

  if test "x$cs_have_catalyst" = "xno"; then
    CATALYST_LIBS=""
  else
    AC_LANG_PUSH([C++])
    CPPFLAGS="${saved_CPPFLAGS} ${CATALYST_CPPFLAGS}"
    AC_CHECK_HEADERS([vtkCPPythonScriptV2Pipeline.h])
    AC_LANG_POP([C++])
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

unset cs_catalyst_version

if test x$cs_have_catalyst = xno ; then
  cs_have_plugin_catalyst=no
fi

AM_CONDITIONAL(HAVE_CATALYST, test x$cs_have_catalyst = xyes)
AM_CONDITIONAL(HAVE_PLUGIN_CATALYST, test x$cs_have_plugin_catalyst = xyes)

cs_py_have_plugin_catalyst=False
if test x$cs_have_plugin_catalyst = xyes ; then
  cs_py_have_plugin_catalyst=True
fi

AC_SUBST(cs_have_catalyst)
AC_SUBST(cs_py_have_plugin_catalyst)
AC_SUBST(CATALYST_CPPFLAGS)
AC_SUBST(CATALYST_CXXFLAGS)
AC_SUBST(CATALYST_LDFLAGS)
AC_SUBST(CATALYST_LIBS)
AC_SUBST(CATALYSTRUNPATH)

])dnl

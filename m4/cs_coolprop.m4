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

# CS_AC_TEST_COOLPROP
#---------------
# modifies or sets cs_have_coolprop, COOLPROP_CPPFLAGS, COOLPROP_LDFLAGS,
# and COOLPROP_LIBS depending on libraries found

AC_DEFUN([CS_AC_TEST_COOLPROP], [

cs_gui_coolprop=false

cs_have_coolprop=no
cs_have_coolprop_headers=no

AC_ARG_ENABLE(gui-coolprop,
  [AS_HELP_STRING([--enable-gui-coolprop], [Allow CoolProp fluids in GUI even when not available])],
  [
    case "${enableval}" in
      yes) cs_gui_coolprop=true ;;
      no)  cs_gui_coolprop=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-gui-coolprop]) ;;
    esac
  ],
  [ cs_gui_coolprop=false ]
)

AC_ARG_WITH(coolprop,
            [AS_HELP_STRING([--with-coolprop=DIR],
                            [specify prefix directory for COOLPROP])],
            [if test "x$withval" = "x"; then
               with_coolprop=yes
             elif test "x$withval" = "xsalome"; then
               if test -z "$COOLPROPHOME"; then
                 AC_MSG_FAILURE([no SALOME path information for COOLPROP (needed by --with-coolprop=salome)!])
               else
                 with_coolprop=$COOLPROPHOME
               fi
             fi],
            [with_coolprop=check])

AC_ARG_WITH(coolprop-include,
            [AS_HELP_STRING([--with-coolprop-include=DIR],
                            [specify directory for CoolProp include files])],
            [if test "x$with_coolprop" = "xcheck"; then
               with_coolprop=yes
             fi
             COOLPROP_CPPFLAGS="-I$with_coolprop_include"],
            [if test "x$with_coolprop" != "xno" -a "x$with_coolprop" != "xyes" \
	          -a "x$with_coolprop" != "xcheck"; then
               if test -f "$with_coolprop/include/CoolProp.h"; then
                 COOLPROP_CPPFLAGS="-I$with_coolprop/include"
               elif test -f "$with_coolprop/CoolProp.h"; then
                 COOLPROP_CPPFLAGS="-I$with_coolprop"
               else
                 find $with_coolprop | grep CoolProp.h > /dev/null 2>&1;
                 if test $? == "0"; then
                   cs_coolprop_tmp=`find $with_coolprop -name CoolProp.h | head -1`
                   COOLPROP_CPPFLAGS="-I`dirname $cs_coolprop_tmp`"
                   unset cs_coolprop_tmp
                 fi
               fi
             fi])

AC_ARG_WITH(coolprop-lib,
            [AS_HELP_STRING([--with-coolprop-lib=DIR],
                            [specify directory for CoolProp library])],
            [if test "x$with_coolprop" = "xcheck"; then
               with_coolprop=yes
             fi
             COOLPROP_LDFLAGS="-L$with_coolprop_lib"
             # Add the libdir to the runpath as CoolProp is not libtoolized
             COOLPROPRUNPATH="-R$with_coolprop_lib"],
            [if test "x$with_coolprop" != "xno" -a "x$with_coolprop" != "xyes" \
	          -a "x$with_coolprop" != "xcheck"; then
               case `uname -m` in
                  ppc64)      ref_name='libCoolProp' ;;
                  *64)        ref_name='64bit/libCoolProp' ;;
                  *32 | i*86) ref_name='32bit/libCoolProp' ;;
                  *)          ref_name='libCoolProp' ;;
               esac
               find $with_coolprop | grep ${ref_name}.* > /dev/null 2>&1;
               if test $? == "0"; then
                 cp_l_name=`find $with_coolprop | grep ${ref_name}.* | head -1`
                 cp_d_name=`dirname ${cp_l_name}`
                 unset cp_l_name
                 COOLPROP_LDFLAGS="-L$cp_d_name"
                 # Add the libdir to the runpath as CoolProp is not libtoolized
                 COOLPROPRUNPATH="-R$cp_d_name"
                 unset cp_d_name
               else
                 COOLPROP_LDFLAGS="-L$with_coolprop/lib"
                 # Add the libdir to the runpath as CoolProp is not libtoolized
                 COOLPROPRUNPATH="-R$with_coolprop/lib"
               fi
             fi])

AC_ARG_WITH(coolprop-pythonpath,
            [AS_HELP_STRING([--with-coolprop-pythonpath=DIR],
                            [specify directory for CoolProp Python bindings])],
            [if test "x$with_coolprop" = "xcheck"; then
               with_coolprop=yes
             fi
             COOLPROPPYTHONPATH="$with_coolprop_pythonpath"],
            [if test "x$with_coolprop" != "xno" -a "x$with_coolprop" != "xyes" \
	          -a "x$with_coolprop" != "xcheck"; then
               find $with_coolprop -name CoolProp.py > /dev/null 2>&1;
               if test $? == "0"; then
                 cs_coolprop_tmp=`find $with_coolprop -name CoolProp.py | head -1`
                 if test "x$cs_coolprop_tmp" != x ; then
                   COOLPROPPYTHONPATH="`dirname $cs_coolprop_tmp`"
                 fi
                 unset cs_coolprop_tmp
               fi
             fi])

if test "x$with_coolprop" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${COOLPROP_CPPFLAGS}"

  COOLPROP_LIBS="-lCoolProp"
  LIBS="${COOLPROP_LIBS} ${LIBS}"

  # Check that CoolProp files exist

  # CoolProp is in C++, though it also provides a C wrapper.
  AC_LANG_PUSH([C++])

  AC_MSG_CHECKING([for CoolProp library)])

  for coolprop_ldadd in "" "-ldl"
  do
    if test "x$cs_have_coolprop" != "xyes"; then
      LDFLAGS="${saved_LDFLAGS} ${COOLPROP_LDFLAGS} ${coolprop_ldadd}"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "CoolProp.h"]],
                                      [[std::vector<double> t;
                                        std::vector<double> p;
                                        std::vector<double> z;
                                        std::vector<std::string> f;
                                        std::vector<std::string> o;
                                        CoolProp::PropsSImulti(o, "T", t, "P", p, "", f, z);]])],
                                      [ AC_DEFINE([HAVE_COOLPROP], 1, [CoolProp support])
                                        cs_have_coolprop=yes],
                                       [cs_have_coolprop=no])
    fi
  done

  AC_MSG_RESULT($cs_have_coolprop)

  AC_LANG_POP

  if test "x$cs_have_coolprop" != "xyes"; then
    if test "x$with_coolprop" != "xcheck"; then
      AC_MSG_FAILURE([CoolProp support requested, but test for CoolProp failed])
    fi
    COOLPROP_CPPFLAGS=""
    COOLPROP_LDFLAGS=""
    COOLPROP_LIBS=""
    COOLPROPRUNPATH=""
    COOLPROPPYTHONPATH=""
    if test "x$cs_gui_coolprop" != "xfalse"; then
      cs_have_coolprop=gui_only
    fi
  fi

  CPPFLAGS="$saved_CPPFLAGS"
  LDFLAGS="$saved_LDFLAGS"
  LIBS="$saved_LIBS"

  unset saved_CPPFLAGS
  unset saved_LDFLAGS
  unset saved_LIBS

fi

AM_CONDITIONAL(HAVE_COOLPROP, test x$cs_have_coolprop = xyes)

AC_SUBST(cs_have_coolprop)
AC_SUBST(coolprop_prefix, [${with_coolprop}])
AC_SUBST(COOLPROP_CPPFLAGS)
AC_SUBST(COOLPROP_LDFLAGS)
AC_SUBST(COOLPROP_LIBS)
AC_SUBST(COOLPROPRUNPATH)
AC_SUBST(COOLPROPPYTHONPATH)

# AC_LANG_POP([C++])

])dnl


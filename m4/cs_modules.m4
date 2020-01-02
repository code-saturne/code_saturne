dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2020 EDF S.A.
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

# CS_AC_TEST_ENV_MODULES
#-----------------------
# checks for environement modules

AC_DEFUN([CS_AC_TEST_ENV_MODULES], [

AC_ARG_WITH(modules,
            [AS_HELP_STRING([--with-modules=LIST],
                            [colon-separated list of environment modules])],
            [with_modules=$withval],
            [with_modules=check])

# Attempt at auto-detection

cs_env_modules="no"

if test "x$with_modules" = "xcheck" ; then

  # Test for environment modules

  if test "x$MODULESHOME" != "x" ; then

    AC_MSG_CHECKING([for environment modules])

    cs_env_modules=""
    try_modules=""
    try_modules_p=""

    outfile=cs_ac_config_modules-tmp

    cs_abs_srcdir=`cd $srcdir && pwd`
    "$SHELL" "$cs_abs_srcdir/build-aux/list_modules.sh" "$outfile" >&5

    cs_env_modules=`cat $outfile`
    rm -fr $outfile

  fi

elif test "x$with_modules" != "xno" ; then

  cs_env_modules=""
  oldIFS=$IFS; IFS=:
  for m in $with_modules; do cs_env_modules="$cs_env_modules $m"; done
  IFS=$oldIFS

fi

# Find the modulecmd executable; handle both modules and lmod

if test "x$LMOD_CMD" != "x" ; then
  MODULECMD=$LMOD_CMD
  AC_SUBST(MODULECMD)
  "$MODULECMD" python avail > /dev/null 2>&1
  if test $? != 0 ; then
    AC_MSG_WARN([LMOD_CMD defined but module commmand does not seem usable])
    cs_env_modules="no"
  fi
elif test "x$MODULESHOME" != "x" ; then
  AC_PATH_PROG([MODULECMD], [modulecmd], [], [${MODULESHOME}/bin:$PATH])
else
  AC_MSG_WARN([no supported environment module commmand detected])
fi

AC_SUBST(cs_env_modules)

])dnl


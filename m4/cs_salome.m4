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

# CS_AC_SALOME_ENV
#-----------------
# Macro determining the presence of a Salome environment

AC_DEFUN([CS_AC_SALOME_ENV], [

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

  if test ! -d "$with_salome" ; then
    AC_MSG_FAILURE([directory specified by --with-salome=$with_salome does not exist!])
  fi

  # Recommended environment file for salome-platform.org installer builds
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 2 -name salome.sh | tail -1 2>/dev/null)
    if test "x$salome_env" != "x"; then
      SALOMEENVCMD="source $salome_env"
    fi
  fi

  # Environment for EDF or "universal binaries" type build for Salome
  # Note that ROOT_SALOME is required but not exported in all cases.
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 1 -name salome_modules.sh 2>/dev/null)
    if test "x$salome_env" != "x"; then
      salome_pre=$(find $with_salome -maxdepth 1 -name salome_prerequisites.sh 2>/dev/null)
      if test "x$salome_pre" != "x"; then
        SALOMEENVCMD="source $salome_pre; export ROOT_SALOME=$with_salome; source $salome_env"
      fi
    fi
  fi

  # Environment for CAS (salome-platform.org) builds for Salome
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env="${with_salome}/env_launch.sh"
    if test -f "$salome_env" ; then
      SALOMEENVCMD="source $salome_env"
    fi
  fi

  unset salome_pre
  unset salome_env

  # Paths for libraries provided by SALOME distibution, for automatic checks

  if test "x$SALOMEENVCMD" != "x" ; then

    (/bin/bash -c "$SALOMEENVCMD ; env > conftest.salome_env")

    if test -z "$MEDCOUPLING_ROOT_DIR" ; then
      MEDCOUPLING_ROOT_DIR=$(grep MEDCOUPLING_ROOT_DIR conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$HDF5HOME" ; then
      HDF5HOME=$(grep HDF5HOME conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$MEDHOME" ; then
      MEDHOME=$(grep MEDHOME conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$CGNSHOME" ; then
      CGNSHOME=$(grep CGNSHOME conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$CATALYST_ROOT_DIR" ; then
      CATALYST_ROOT_DIR=$(grep CATALYST_ROOT_DIR conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$COOLPROPHOME" ; then
      COOLPROPHOME=$(grep COOLPROPHOME conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$METISDIR" ; then
      METISDIR=$(grep METISDIR conftest.salome_env | cut -f2 -d'=')
    fi

    if test -z "$SCOTCHDIR" ; then
      SCOTCHDIR=$(grep SCOTCHDIR conftest.salome_env | cut -f2 -d'=')
    fi

    \rm -rf conftest.salome_env

  fi

  AC_ARG_VAR([SALOMEENVCMD], [SALOME environment setting commands])

fi

])dnl

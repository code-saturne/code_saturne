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

  unset salome_pre
  unset salome_env

  # Paths for libraries provided by SALOME distibution, for automatic checks

  if test -z "$MEDCOUPLING_ROOT_DIR" ; then
    if test "x$SALOMEENVCMD" != "x" ; then
      MEDCOUPLING_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $MEDCOUPLING_ROOT_DIR)
    else
      MEDCOUPLING_ROOT_DIR=$(echo $MEDCOUPLING_ROOT_DIR)
    fi
  fi

  if test -z "$HDF5HOME" ; then
    if test "x$SALOMEENVCMD" != "x" ; then
      HDF5HOME=$(eval $SALOMEENVCMD ; echo $HDF5HOME)
    else
      HDF5HOME=$(echo $HDF5HOME)
    fi
  fi

  if test -z "$MED3HOME" ; then
    if test "x$SALOMEENVCMD" != "x" ; then
      MED3HOME=$(eval $SALOMEENVCMD ; echo $MED3HOME)
    else
      MED3HOME=$(echo $MED3HOME)
    fi
  fi

  if test -z "$CGNSHOME" ; then
    if test "x$SALOMEENVCMD" != "x" ; then
      CGNSHOME=$(eval $SALOMEENVCMD ; echo $CGNSHOME)
    else
      CGNSHOME=$(echo $CGNSHOME)
    fi
  fi

  if test -z "$CATALYST_ROOT_DIR" ; then
    if test "x$SALOMEENVCMD" != "x" ; then
      CATALYST_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $CATALYST_ROOT_DIR)
    else
      CATALYST_ROOT_DIR=$(echo $CATALYST_ROOT_DIR)
    fi
  fi

  AC_ARG_VAR([SALOMEENVCMD], [SALOME environment setting commands])

fi

])dnl

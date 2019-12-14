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

# CS_AC_TEST_SALOME
#------------------
# Macro function grouping the different tests
# Requires CS_AC_TEST_SALOME_ENV be called first

AC_DEFUN([CS_AC_TEST_SALOME], [

AC_REQUIRE([CS_AC_SALOME_ENV])

if test x$with_salome != xno ; then

  if test "x$SALOMEENVCMD" != "x" ; then
    KERNEL_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $KERNEL_ROOT_DIR)
    GUI_ROOT_DIR=$(eval $SALOMEENVCMD ; echo $GUI_ROOT_DIR)
    OMNIIDL=$(eval $SALOMEENVCMD ; which omniidl)
  fi

  # Make sure omniidl will work by forcing PYTHONPATH
  # Note that we must ensure we are using SALOME's python here, so use "python"
  # with the sourced environment and PATH rather than $PYTHON here.

  if test "x$SALOMEENVCMD" != "x" ; then
    if test "x$OMNIIDLPYTHONPATH" = "x"; then
     OMNIIDLPYTHONPATH=$(eval $SALOMEENVCMD ; python -B "$srcdir/build-aux/cs_config_test.py" pythonpath_filter _omniidlmodule.so _omniidlmodule.so)
    fi
    if test "x$OMNIIDLLDLIBPATH" = "x"; then
      OMNIIDLLDLIBPATH=$(eval $SALOMEENVCMD ; python -B "$srcdir/build-aux/cs_config_test.py" ld_library_path_filter libpython*)
    fi
  else
    if test "x$OMNIIDLPYTHONPATH" = "x"; then
     OMNIIDLPYTHONPATH=$(python -B "$srcdir/build-aux/cs_config_test.py" pythonpath_filter _omniidlmodule.so _omniidlmodule.so)
    fi
    if test "x$OMNIIDLLDLIBPATH" = "x"; then
      OMNIIDLLDLIBPATH=$(python -B "$srcdir/build-aux/cs_config_test.py" ld_library_path_filter libpython*)
    fi
  fi

fi

AC_SUBST(OMNIIDLPYTHONPATH)
AC_SUBST(OMNIIDLLDLIBPATH)

CS_AC_TEST_SALOME_KERNEL
CS_AC_TEST_SALOME_GUI

AC_LANG_SAVE

omniORB_ok=no # in case following test is not called

AS_IF([test $cs_have_salome_kernel = yes -o $cs_have_salome_gui = yes],
      [CS_AC_TEST_OMNIORB])

AC_LANG_RESTORE

AM_CONDITIONAL(HAVE_SALOME,
               test $cs_have_salome_kernel = yes \
                 -a $cs_have_salome_gui = yes \
                 -a $omniORB_ok = yes)

])dnl

# CS_AC_TEST_SALOME_KERNEL
#-------------------------
# modifies or sets cs_have_salome_kernel, SALOME_KERNEL_CPPFLAGS, SALOME_KERNEL_LDFLAGS
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
               with_salome_kernel=no
             else
               with_salome_kernel=$KERNEL_ROOT_DIR
             fi])

if test "x$with_salome_kernel" != "xno" ; then

  cs_have_salome_kernel=yes

  if test x"$with_salome_kernel" != xyes -a x"$with_salome_kernel" != xcheck ; then
    SALOME_KERNEL="$with_salome_kernel"
  else
    SALOME_KERNEL="/usr"
  fi

  SALOME_KERNEL_CPPFLAGS="-I$SALOME_KERNEL/include/salome"
  SALOME_KERNEL_IDL="-I$SALOME_KERNEL/idl/salome"
  SALOME_KERNEL_LDFLAGS="-L$SALOME_KERNEL/lib/salome"

fi

AC_SUBST(SALOME_KERNEL)
AC_SUBST(SALOME_KERNEL_CPPFLAGS)
AC_SUBST(SALOME_KERNEL_IDL)
AC_SUBST(SALOME_KERNEL_LDFLAGS)

AM_CONDITIONAL(HAVE_SALOME_KERNEL, test x$cs_have_salome_kernel = xyes)

])dnl


# CS_AC_TEST_SALOME_GUI
#----------------------
# modifies or sets cs_have_salome_gui, SALOME_GUI_CPPFLAGS, SALOME_GUI_LDFLAGS
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
               with_salome_gui=no
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

else
  cs_have_salome_gui=no
fi

AC_SUBST(cs_have_salome_gui)
AC_SUBST(SALOME_GUI)
AC_SUBST(SALOME_GUI_CPPFLAGS)
AC_SUBST(SALOME_GUI_IDL)
AC_SUBST(SALOME_GUI_LDFLAGS)

AM_CONDITIONAL(HAVE_SALOME_GUI, test x$cs_have_salome_gui = xyes)

])dnl

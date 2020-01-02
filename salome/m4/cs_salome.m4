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

# SALOME_CFD_AC_SALOME_ENV
#-------------------------
# Macro determining the presence of a Salome environment

AC_DEFUN([SALOME_CFD_AC_SALOME_ENV], [

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

  # Check for Salome modules
  #-------------------------

  salome_env_modules=`${srcdir}/build-aux/list_salome_modules.py ${ROOT_SALOME}`
  AC_SUBST(salome_env_modules)

  if test "$xsalome_env_modules" != "x" ; then
     echo "Modules required by Salome: ${salome_env_modules}"
     echo
     echo "These modules should be loaded before running configure and make."
  fi

  # Recommended environment file for salome-platform.org installer builds
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 2 -name salome.sh | tail -1 2>/dev/null)
    if test "x$salome_env" != "x"; then
      SALOMEENVCMD=". $salome_env"
    fi
  fi

  # Environment for EDF or "universal binaries" type build for Salome
  # Note that ROOT_SALOME is required but not exported in all cases.
  if test "x$SALOMEENVCMD" = "x"; then
    salome_env=$(find $with_salome -maxdepth 1 -name salome_modules.sh 2>/dev/null)
    if test "x$salome_env" != "x"; then
      salome_pre=$(find $with_salome -maxdepth 1 -name salome_prerequisites.sh 2>/dev/null)
      if test "x$salome_pre" != "x"; then
        SALOMEENVCMD=". $salome_pre; export ROOT_SALOME=$with_salome; . $salome_env"
      fi
    fi
  fi

  unset salome_pre
  unset salome_env

  AC_ARG_VAR([SALOMEENVCMD], [SALOME environment setting commands])

fi

])dnl

# SALOME_CFD_AC_TEST_SALOME
#--------------------------
# Macro function grouping the different tests
# Requires SALOME_CFD_AC_TEST_SALOME_ENV be called first

AC_DEFUN([SALOME_CFD_AC_TEST_SALOME], [

AC_REQUIRE([SALOME_CFD_AC_SALOME_ENV])

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
     OMNIIDLPYTHONPATH=$(eval $SALOMEENVCMD ; python -B "$srcdir/config/cs_config_test.py" pythonpath_filter _omniidlmodule.so _omniidlmodule.so)
    fi
    if test "x$OMNIIDLLDLIBPATH" = "x"; then
      OMNIIDLLDLIBPATH=$(eval $SALOMEENVCMD ; python -B "$srcdir/config/cs_config_test.py" ld_library_path_filter libpython*)
    fi
  else
    if test "x$OMNIIDLPYTHONPATH" = "x"; then
     OMNIIDLPYTHONPATH=$(python -B "$srcdir/config/cs_config_test.py" pythonpath_filter _omniidlmodule.so _omniidlmodule.so)
    fi
    if test "x$OMNIIDLLDLIBPATH" = "x"; then
      OMNIIDLLDLIBPATH=$(python -B "$srcdir/config/cs_config_test.py" ld_library_path_filter libpython*)
    fi
  fi

fi

AC_SUBST(OMNIIDLPYTHONPATH)
AC_SUBST(OMNIIDLLDLIBPATH)

if test "x$OMNIIDL" = "x" ; then
  AC_MSG_WARN([omniIDL not found: SALOME KERNEL and GUI support cannot be built])
  with_salome_kernel=no
  with_salome_gui=no
fi

SALOME_CFD_AC_TEST_SALOME_KERNEL
SALOME_CFD_AC_TEST_SALOME_GUI

AC_LANG_SAVE

omniORB_ok=no # in case following test is not called

AS_IF([test $salome_cfd_have_kernel = yes -o $salome_cfd_have_gui = yes],
      [CS_AC_TEST_OMNIORB])

AC_LANG_RESTORE

AM_CONDITIONAL(HAVE_SALOME,
               test $salome_cfd_have_kernel = yes \
                 -a $salome_cfd_have_gui = yes \
                 -a $omniORB_ok = yes)

])dnl

# SALOME_CFD_AC_TEST_SALOME_KERNEL
#---------------------------------
# modifies or sets salome_cfd_have_kernel, SALOME_KERNEL_CPPFLAGS, SALOME_KERNEL_LDFLAGS
# depending on libraries found

AC_DEFUN([SALOME_CFD_AC_TEST_SALOME_KERNEL], [

salome_cfd_have_kernel=no

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

  salome_cfd_have_kernel=yes

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

AM_CONDITIONAL(HAVE_SALOME_KERNEL, test x$salome_cfd_have_kernel = xyes)

])dnl


# SALOME_CFD_AC_TEST_SALOME_GUI
#----------------------
# modifies or sets salome_cfd_have_gui, SALOME_GUI_CPPFLAGS, SALOME_GUI_LDFLAGS
# depending on libraries found

AC_DEFUN([SALOME_CFD_AC_TEST_SALOME_GUI], [

salome_cfd_have_gui=no

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

  salome_cfd_have_gui=yes

  if test x"$with_salome_gui" != xyes -a x"$with_salome_gui" != xcheck ; then
    SALOME_GUI="$with_salome_gui"
  else
    SALOME_GUI="/usr"
  fi

  SALOME_GUI_CPPFLAGS="-I$SALOME_GUI/include/salome"
  SALOME_GUI_IDL="-I$SALOME_GUI/idl/salome"
  SALOME_GUI_LDFLAGS="-L$SALOME_GUI/lib/salome"

else
  salome_cfd_have_gui=no
fi

AC_SUBST(salome_cfd_have_gui)
AC_SUBST(SALOME_GUI)
AC_SUBST(SALOME_GUI_CPPFLAGS)
AC_SUBST(SALOME_GUI_IDL)
AC_SUBST(SALOME_GUI_LDFLAGS)

AM_CONDITIONAL(HAVE_SALOME_GUI, test x$salome_cfd_have_gui = xyes)

])dnl

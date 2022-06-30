dnl Copyright (C) 2005-2022 EDF
dnl
dnl This file is part of the PLE software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl PLE source distribution.

# PLE_AC_TEST_PYTHON_DEV
# ----------------------
# Check for the presence of Python devel .h files

AC_DEFUN([PLE_AC_PYTHON_DEV],[
  AC_ARG_VAR([PYTHON_VERSION],[])
  AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
  if test -z "$PYTHON"; then
    AC_MSG_ERROR([User asked for python $PYTHON_VERSION, but it was not found.])
  fi

  # Check for distutils
  AC_MSG_CHECKING([if distutils is available])
  ac_distutils_av=`$PYTHON -c "import distutils" 2>&1`
  if test $? -eq 0; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR(["distutils" package is missing])
  fi

  # Find the python include path
  AC_MSG_CHECKING([for Python *.h path])
  # Get generic include path
  generic_python_inc=`$PYTHON -c "import distutils.sysconfig as dsyscfg; \
    print(dsyscfg.get_python_inc());"`
  # Get platform specific include path
  pf_python_inc=`$PYTHON -c "import distutils.sysconfig as dsyscfg; \
    print(dsyscfg.get_python_inc(plat_specific=1));"`

  # First check for generic path
  if test ! -f "${generic_python_inc}/Python.h"; then
    # If not there, then check for platform specific path
    if test ! -f "${pf_python_inc}/Python.h"; then
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([could not find Python.h file. Please check that python-devel is installed.])
    else
      AC_MSG_RESULT([yes])
    fi
  else
    AC_MSG_RESULT([yes])
  fi


  AC_MSG_CHECKING([for python cppflags update if needed])
  # Update the cppflags if needed.
  if test -z "$PYTHON_CPPFLAGS"; then
    if test -n "${generic_python_inc}"; then
      if test "${pf_python_inc}" != "${generic_python_inc}"; then
        python_inc_path="-I$generic_python_inc -I$pf_python_inc"
      else
        python_inc_path="-I$generic_python_inc"
      fi
    fi
    PYTHON_CPPFLAGS=$python_inc_path
  fi
  AC_MSG_RESULT([$PYTHON_CPPFLAGS])
  AC_SUBST([PYTHON_CPPFLAGS])

])


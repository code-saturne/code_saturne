dnl Copyright (C) 2005-2022 EDF
dnl
dnl This file is part of the PLE software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl PLE source distribution.

# PLE_AC_TEST_MPI4PY
#----------------------
# modifies or sets ple_have_mpi4py

AC_ARG_WITH(mpi4py,
            [AS_HELP_STRING([--with-mpi4py=PATH],
                            [specify directory for mpi4py])],
            [if test "x$withval" = "x"; then
               with_mpi4py=yes
             fi],
            [with_mpi4py=check])

if test "$with_mpi4py" != no; then
  if test -d "$with_mpi4py" ; then
    MPI4PY="$with_mpi4py"
    AC_SUBST(MPI4PY)
  else
    AC_MSG_FAILURE([directory specified by --with-mpi4py=$with_mpi4py does not exist!])
  fi
fi

AC_DEFUN([PLE_AC_TEST_MPI4PY], [
  if test -z $PYTHON;
  then
    if test -z "$2";
    then
      PYTHON="python3"
    else
      PYTHON="$2"
    fi
  fi
  PYTHON_NAME=`$PYTHON`
  AC_MSG_CHECKING("$PYTHON_NAME module: mpi4py")

  if test "$with_mpi4py" != no;
  then
    mpi4py_import_cmd="import sys; sys.path.insert(0, '$with_mpi4py'); import mpi4py"
  else
    mpi4py_import_cmd="import mpi4py"
  fi
  mpi4py_inc_cmd=""$mpi4py_import_cmd"; print(mpi4py.get_include())"

  $PYTHON -c "$mpi4py_import_cmd" 2>/dev/null

  if test $? -eq 0;
  then
    AC_MSG_RESULT(yes)
    eval AS_TR_CPP(HAVE_PYMOD_mpi4py)=yes
    mpi4py_includes="$($PYTHON -c "$mpi4py_inc_cmd")"
    MPI4PY_CFLAGS="-I${mpi4py_includes}"
    if test "$with_mpi4py" != no; then
      MPI4PY_PATH="${with_mpi4py}"
      AC_SUBST(MPI4PY_PATH)
    else
      MPI4PY_PATH="None"
    fi
  else
    AC_MSG_RESULT(no)
    eval AS_TR_CPP(HAVE_PYMOD_mpi4py)=no
    if test -n "$1"
    then
      AC_MSG_ERROR(failed to find required module mpi4py)
      exit 1
    fi
  fi
  AC_SUBST(MPI4PY_CFLAGS)
])

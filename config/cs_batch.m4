dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2011 EDF S.A., France
dnl
dnl   The Code_Saturne Kernel is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne Kernel is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne Preprocessor; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_TEST_BATCH
#-----------------
# Tries to determine if a batch system is available

AC_DEFUN([CS_AC_TEST_BATCH], [

cs_batch_template=no

AC_ARG_WITH(batch,
            [AS_HELP_STRING([--with-batch=TYPE],
                            [specify batch type or template path])],
            [with_batch=$withval],
            [with_batch=check])

# Attempt at auto-detection

if test "x$with_batch" = "xcheck" ; then

  AC_MSG_CHECKING([for batch system])

  # Check for available batch types

  cmd_prev=no

  for cmd in qsub ccc_msub bsub llsubmit sbatch; do

    which $cmd > /dev/null 2>&1
    if test $? = 0 ; then

      if test "$cmd_prev" != "no" ; then
        AC_MSG_ERROR([
At least 2 batch submission commands found ($cmd_prev and $cmd);
use --with-batch=<value> to specify either batch type (PBS, CCC, LSF,
LOADL, SGE, or SLURM), absolute path to batch template file, or "no")])
      else
        cmd_prev=$cmd
      fi
    fi

  done

  case $cmd_prev in
    qsub)
      # qsub may be either PBS (PBS-PRO, TORQUE) or Sun Grid Engine
      # but Sun Grid Engine should have qmon command.
       which qmon > /dev/null 2>&1
       if test $? = 0 ; then
         cs_batch_template=SGE
       else
         cs_batch_template=PBS
       fi
       ;;

    ccc_msub) cs_batch_template=CCC ;; 
    bsub) cs_batch_template=LSF ;;
    llsubmit) cs_batch_template=LOADL ;;
    sbatch) cs_batch_template=SLURM ;;
    *) cs_batch_template=no ;;
  esac

  AC_MSG_RESULT($cs_batch_template)

elif test "x$with_batch" != "xno" ; then
  cs_batch_template=$with_batch
fi

AC_SUBST(cs_batch_template)

])dnl


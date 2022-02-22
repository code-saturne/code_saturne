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

# CS_AC_TEST_FC_MOD([use_threads])
#------------------
# checks how the Fortran compiler handles modules

AC_DEFUN([CS_AC_TEST_FC_MOD], [

AC_LANG_PUSH(Fortran)

# Create temporary directory

i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ;
do
  i=`expr $i + 1`
done
mkdir tmpdir_$i

# Compile module in temporary directory and check extension
# (some compilers put module filename in uppercase letters,
# so also check this)

cd tmpdir_$i
AC_COMPILE_IFELSE([
      module conftest_module
      contains
      subroutine conftest_routine
      write(*,'(a)') 'conftest'
      end subroutine conftest_routine
      end module conftest_module
  ],
                  [cs_fc_modext=`ls | sed -n 's,conftest_module\.,,p'`
                   if test x$cs_fc_modext = x ; then
                     cs_fc_modext=`ls | sed -n 's,CONFTEST_MODULE\.,,p'`
                     if test x$cs_fc_modext = x ; then
                       cs_fc_modext=""
                     fi
                   fi
                  ],
                  [cs_fc_modext=""])

# Go up one level and search for module using probable flags

cd ..

for cs_fc_flag in "-I" "-I " "-J" "-M" "-p"; do
  if test "x$cs_fc_modflag" = "x" ; then
    save_FCFLAGS="$FCFLAGS"
    FCFLAGS="$save_FCFLAGS ${cs_fc_flag}tmpdir_$i"
    AC_COMPILE_IFELSE([
      program conftest_program
      use conftest_module
      call conftest_routine
      end program conftest_program
      ],
                      [cs_fc_modflag="$cs_fc_flag"],
                      [])
    FCFLAGS="$save_FCFLAGS"
  fi
done

# Now remove temporary directory and finish

rm -fr tmpdir_$i
AC_LANG_POP(Fortran)

FCMODEXT=$cs_fc_modext
FCMODINCLUDE=$cs_fc_modflag

AC_SUBST(FCMODEXT)
AC_SUBST(FCMODINCLUDE)

])dnl


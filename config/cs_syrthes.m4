dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2009 EDF S.A., France
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

# CS_AC_TEST_SYRTHES
#-------------------
# Check for SYRTHES version ;
# defines SYRTHES_CPPFLAGS, SYRTHES_LDFLAGS, and SYRTHES_LIBS locally
# (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CS_AC_TEST_SYRTHES], [

AC_ARG_WITH(syrthes, [AS_HELP_STRING([--with-syrthes=PATH], [specify prefix directory for SYRTHES])])

if test -f $with_syrthes/bin/syrthes.profile ; then
  have_syrthes=yes
  AC_MSG_NOTICE([sourcing $with_syrthes/bin/syrthes.profile])
  . "$with_syrthes/bin/syrthes.profile"
else
  have_syrthes=no
  AC_MSG_WARN([cannot find syrthes.profile])
fi

if test "x$have_syrthes" = "xyes"; then

# Get SYRTHES compilers

outfile=makefile-tmp
cp $with_syrthes/bin/Makefile $outfile

cat >> $outfile <<\_______EOF

version:
	@echo $(VERSION)
cc:
	@echo $(CC)
fc:
	@echo $(FC)
cflags:
	@echo $(CFLAGS)
fcflags:
	@echo $(FCFLAGS)
_______EOF

SYRTHES_NOM_ARCH=${NOM_ARCH}
SYRTHES_VERSION=`make -f $outfile version`
SYRTHES_CC=`make -f $outfile cc`
SYRTHES_FC=`make -f $outfile fc`
SYRTHES_CFLAGS=`make -f $outfile cflags`
SYRTHES_FCFLAGS=`make -f $outfile fcflags`

rm -f $outfile

# Get mandatory Fortran libs for linking stage
# One assumes that SYRTHES is writtent in Fortran 77
F77=$SYRTHES_FC
AC_F77_LIBRARY_LDFLAGS
SYRTHES_FCLIBS=$FLIBS

SYRTHES_CPPFLAGS="-I$with_syrthes/include"
SYRTHES_LDFLAGS="-L$with_syrthes/lib/${SYRTHES_NOM_ARCH}"
SYRTHES_LIBS="-lsatsyrthes${SYRTHES_VERSION}_${SYRTHES_NOM_ARCH} -lsyrthes${SYRTHES_VERSION}_${SYRTHES_NOM_ARCH}"

if test "x$cs_have_mpi" = "xyes"; then
   FVM_COUPL_LDFLAGS="$FVM_LDFLAGS"
   FVM_COUPL_LIBS="-lfvm_coupl"
else
   FVM_COUPL_LDFLAGS=
   FVM_COUPL_LIBS=
fi

SYRTHES_LDFLAGS="$FVM_COUPL_LDFLAGS $BFT_LDFLAGS $MPI_LDFLAGS $SYRTHES_LDFLAGS"
SYRTHES_LIBS="$FVM_COUPL_LIBS $BFT_LIBS $MPI_LIBS $SYRTHES_LIBS"

AC_SUBST(SYRTHES_CC)
AC_SUBST(SYRTHES_FC)
AC_SUBST(SYRTHES_CFLAGS)
AC_SUBST(SYRTHES_FCFLAGS)
AC_SUBST(SYRTHES_FCLIBS)
AC_SUBST(SYRTHES_CPPFLAGS)
AC_SUBST(SYRTHES_LDFLAGS)
AC_SUBST(SYRTHES_LIBS)

fi

AM_CONDITIONAL(HAVE_SYRTHES, [test "x$have_syrthes" = "xyes"])

])dnl

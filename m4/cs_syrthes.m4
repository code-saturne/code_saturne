dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2011 EDF S.A.
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

# CS_AC_TEST_SYRTHES
#-------------------
# Check for SYRTHES version ;
# defines SYRTHES_CPPFLAGS, SYRTHES_LDFLAGS, and SYRTHES_LIBS locally
# (i.e. as simple variables, not AC_SUBST)

AC_DEFUN([CS_AC_TEST_SYRTHES], [

AC_ARG_WITH(syrthes, [AS_HELP_STRING([--with-syrthes=PATH], [specify prefix directory for SYRTHES])])

# First try with syrthes.profile, second try for Debian-like packaging
# of SYRTHES where there is no syrthes.profile

AC_MSG_CHECKING([for SYRTHES support])
if test -f $with_syrthes/bin/syrthes.profile ; then
  cs_have_syrthes=yes
  syrthes_prefix=$with_syrthes
  AC_MSG_RESULT([sourcing $syrthes_prefix/bin/syrthes.profile])
  . "$syrthes_prefix/bin/syrthes.profile"
elif test -f $with_syrthes/bin/syrthes_create_case ; then
  cs_have_syrthes=yes
  syrthes_prefix=$with_syrthes
  AC_MSG_RESULT([found $syrthes_prefix/bin/syrthes_create_case])
else
  cs_have_syrthes=no
  AC_MSG_WARN([cannot find syrthes.profile])
fi

if test "x$cs_have_syrthes" = "xyes"; then

# Get SYRTHES compilers

outfile=makefile-tmp
cp $with_syrthes/bin/Makefile $outfile

cat >> $outfile <<\_______EOF

syr_info:
	@echo $(NOM_ARCH) > syr-nomarch-tmp
	@echo $(VERSION) > syr-version-tmp
	@echo $(CC) > syr-cc-tmp
	@echo $(FC) > syr-fc-tmp
	@echo $(CFLAGS) > syr-cflags-tmp
	@echo $(FCFLAGS) > syr-fcflags-tmp
_______EOF

make -f $outfile syr_info > /dev/null

SYRTHES_NOM_ARCH=`cat syr-nomarch-tmp`
SYRTHES_VERSION=`cat syr-version-tmp`
SYRTHES_CC=`cat syr-cc-tmp`
SYRTHES_FC=`cat syr-fc-tmp`
SYRTHES_CFLAGS=`cat syr-cflags-tmp`
SYRTHES_FCFLAGS=`cat syr-fcflags-tmp`

rm -f $outfile syr-nomarch-tmp syr-version-tmp syr-cc-tmp syr-fc-tmp syr-cflags-tmp syr-fcflags-tmp

# Get mandatory Fortran libs for linking stage
# We assume that only SYRTHES user subroutines written in Fortran 77
# will need to be compiled (as versions older than 3.4 are not
# supported, and version 4 will be coupled directly through PLE,
# not through the syrcs wrapper).

F77=$SYRTHES_FC
AC_F77_LIBRARY_LDFLAGS
SYRTHES_FCLIBS=$FLIBS

SYRTHES_CPPFLAGS="-I$with_syrthes/include"
SYRTHES_LDFLAGS="-L$with_syrthes/lib/${SYRTHES_NOM_ARCH}"
SYRTHES_LIBS="-lsatsyrthes${SYRTHES_VERSION}_${SYRTHES_NOM_ARCH} -lsyrthes${SYRTHES_VERSION}_${SYRTHES_NOM_ARCH}"

SYRTHES_LDFLAGS="$PLE_LDFLAGS $MPI_LDFLAGS $SYRTHES_LDFLAGS"
SYRTHES_LIBS="$PLE_LIBS $MPI_LIBS $SYRTHES_LIBS"

AC_SUBST(syrthes_prefix)

AC_SUBST(SYRTHES_CC)
AC_SUBST(SYRTHES_FC)
AC_SUBST(SYRTHES_CFLAGS)
AC_SUBST(SYRTHES_FCFLAGS)
AC_SUBST(SYRTHES_FCLIBS)
AC_SUBST(SYRTHES_CPPFLAGS)
AC_SUBST(SYRTHES_LDFLAGS)
AC_SUBST(SYRTHES_LIBS)

fi

])dnl

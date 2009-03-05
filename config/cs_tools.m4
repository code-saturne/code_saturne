dnl Copyright (C) 2009 EDF S.A., France
dnl
dnl This file is part of the Code_Saturne Preprocessor.  For license
dnl information, see the COPYING file in the top level directory of the
dnl Code_Saturne Preprocessor source distribution.

# CS_AC_TEST_PREPRO(Minimal Release string, [Maximal Release string])
#--------------------------------------------------------------------
# Check for the Preprocessor executable ; defines CSPP_HOME

AC_DEFUN([CS_AC_TEST_PREPRO], [

AC_ARG_WITH(prepro, [AS_HELP_STRING([--with-prepro=PATH], [specify prefix directory for the Preprocessor])])

AC_CHECK_PROG(CSPP_HOME, ecs, $with_prepro, , $with_prepro/bin, )
AC_SUBST(CSPP_HOME)

])dnl


# CS_AC_TEST_GUI(Minimal Release string, [Maximal Release string])
#-----------------------------------------------------------------
# Check for the GUI executable ; defines CSGUIHOME

AC_DEFUN([CS_AC_TEST_GUI], [

AC_ARG_WITH(gui, [AS_HELP_STRING([--with-gui=PATH], [specify prefix directory for the GUI])])

AC_CHECK_PROG(CSGUI_HOME, ics, $with_gui, , $with_gui, )
AC_SUBST(CSGUI_HOME)

])dnl


# CS_AC_TEST_SYRCS(Minimal Release string, [Maximal Release string])
#-------------------------------------------------------------------
# Check for the SYR_CS library ; defines SYRCS_HOME

AC_DEFUN([CS_AC_TEST_SYRCS], [

AC_ARG_WITH(syrcs, [AS_HELP_STRING([--with-syrcs=PATH], [specify prefix directory for Syr_CS library])])

if test "x$with_syrcs" != "x" ; then
  SYRCS_LDFLAGS="-L$with_syrcs/lib/Linux"
fi

saved_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $SYRCS_LDFLAGS"

AC_CHECK_LIB(syr_cs, syr_coupling_initialize, SYRCS_HOME=$with_syrcs, , )

LDFLAGS="$saved_LDFLAGS"


AC_SUBST(SYRCS_HOME)



])dnl

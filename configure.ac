# Process this file with autoconf to produce a configure script

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([cs_licence_c_comment],
[/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/])

m4_define([cs_major_version], [2])
m4_define([cs_minor_version], [2])
m4_define([cs_release_version], [2])
m4_define([cs_version_extra], [])
m4_define([cs_version_string],
          [cs_major_version.cs_minor_version.cs_release_version@&t@cs_version_extra])

AC_INIT([code_saturne],[cs_version_string],[saturne-support@edf.fr],[],[http://www.code-saturne.org])
AC_CONFIG_SRCDIR([src/apps/cs_solver.c])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([cs_config.h])

AM_GNU_GETTEXT_VERSION([0.18.1])

# Initialize automake with the following options:
# - foreign, so as to handle the absence of ChangeLog (automatically generated)
# - tar-pax, so as to handle long lines (> 99 characters) in tar archives
# Warnings can be activated at bootstrap with 'autoreconf -vi --warnings=all'
AM_INIT_AUTOMAKE([foreign tar-pax])

# Enable maintainer mode by default for a developper checkout
AS_IF([test -d ${srcdir}/.svn -o -d ${srcdir}/.git],
      [AM_MAINTAINER_MODE([enable])],
      [AM_MAINTAINER_MODE([disable])])

# Enable silent rules with "./configure --enable-silent-rules" or "make V=0"
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

# Check for environment modules
#------------------------------

CS_AC_TEST_ENV_MODULES

# Check for C and Fortran compilers
#----------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_CXXFLAGS=$CXXFLAGS
user_FCFLAGS=$FCFLAGS
user_LDFLAGS=$LDFLAGS
user_LIBS=$LIBS

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC

# Search for Fortran 77 compilers just in case we have SYRTHES coupling
# To be removed as soon SYRTHES 4 becomes the official version
AC_PROG_F77

AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_CC_C_O
AC_PROG_SED
AM_PROG_LEX
AC_PROG_YACC

AC_FC_LIBRARY_LDFLAGS

# Debug or production compilation mode (debug by default) ?

AC_ARG_ENABLE(debug,
  [AS_HELP_STRING([--enable-debug], [enable debugging (reduces optimization)])],
  [
    case "${enableval}" in
      yes) debug=yes ;;
      no)  debug=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac
  ],
  [ debug=no ]
)
AC_SUBST(debug)

# Optionally deactivate automatic determination of flags on known systems

AC_ARG_ENABLE(auto-flags,
  [AS_HELP_STRING([--disable-auto-flags], [do not define *FLAGS on known systems])],
  [
    case "${enableval}" in
      yes) auto_flags=yes ;;
      no)  auto_flags=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-auto-flags]) ;;
    esac
  ],
  [ auto_flags=yes ]
)

# Optionally activate French translation of Fortran logs (compile-time)

AC_ARG_ENABLE(french,
  [AS_HELP_STRING([--enable-french], [enable French translation of Fortran logs])],
  [
    case "${enableval}" in
      yes) french=yes ;;
      no)  french=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-french]) ;;
    esac
  ],
  [ french=no ]
)

# Optionnaly install the code with relocatable features

AC_ARG_ENABLE(relocatable,
  [AS_HELP_STRING([--enable-relocatable], [enable relocatable installation])],
  [
    case "${enableval}" in
      yes) relocatable=yes ;;
      no)  relocatable=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-relocatable]) ;;
    esac
  ],
  [ relocatable=no ]
)
AC_SUBST(relocatable)
AM_CONDITIONAL(HAVE_RELOCATABLE, [test "${relocatable}" = yes])

# Default compiler options (may be modified
# by defining CFLAGS in the environment)
#------------------------------------------

if test "x$auto_flags" = "xyes" ; then

  # Source associated recommended compiler options

  if test -f "$srcdir/config/cs_auto_flags.sh" ; then
    AC_MSG_NOTICE([sourcing config/cs_auto_flags.sh (test for known compilers)])
    . "$srcdir/config/cs_auto_flags.sh"
  else
    AC_MSG_WARN([config/cs_auto_flags.sh default configuration file not found])
  fi

  # Default flags
  CPPFLAGS="$cppflags_default $user_CPPFLAGS"
  CFLAGS="$cflags_default $user_CFLAGS"
  CXXFLAGS="$cxxflags_default $user_CXXFLAGS"
  FCFLAGS="$fcflags_default $user_FCFLAGS"
  LDFLAGS="$ldflags_default $user_LDFLAGS"
  LIBS="$libs_default $user_LIBS"
  LDRPATH="$ldflags_rpath"

  if test "x$debug" = xyes; then
    # Debug flags
    CFLAGS_DBG="$cflags_default_dbg"
    CXXFLAGS_DBG="$cxxflags_default_dbg"
    FCFLAGS_DBG="$fcflags_default_dbg"
    # Add debug flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_dbg"
  else
    # Normal optimization flags
    CFLAGS_OPT="$cflags_default_opt"
    CXXFLAGS_OPT="$cxxflags_default_opt"
    FCFLAGS_OPT="$fcflags_default_opt"
    # Hot optimization flags
    CFLAGS_HOT="$cflags_default_hot"
    CXXFLAGS_HOT="$cxxflags_default_hot"
    FCFLAGS_HOT="$fcflags_default_hot"
    # Add optimization flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_opt"
  fi
  CFLAGS_EXT="$cflags_default_ext"

  AC_SUBST(CFLAGS_DBG)
  AC_SUBST(CFLAGS_OPT)
  AC_SUBST(CFLAGS_HOT)
  AC_SUBST(CFLAGS_EXT)

  AC_SUBST(CXXFLAGS_DBG)
  AC_SUBST(CXXFLAGS_OPT)
  AC_SUBST(CXXFLAGS_HOT)

  AC_SUBST(FCFLAGS_DBG)
  AC_SUBST(FCFLAGS_OPT)
  AC_SUBST(FCFLAGS_HOT)

  AC_SUBST(LDRPATH)
fi

# Add some preprocessor flags to include additional system features
#------------------------------------------------------------------

case "$host_os" in
  linux*)
    CPPFLAGS="${CPPFLAGS} -D_POSIX_SOURCE"
    ;;
  *bsd* | solaris2.* | aix*)
    ;;
  darwin*)
    CPPFLAGS="${CPPFLAGS} -D_DARWIN_C_SOURCE"
    ;;
  osf*)
    CPPFLAGS="${CPPFLAGS} -D_OSF_SOURCE"
    ;;
  hpux*)
    CPPFLAGS="${CPPFLAGS} -D_HPUX_SOURCE"
    ;;
  *)
    ;;
esac

# Preprocessor flags for debugging purposes

if test "x$debug" = "xyes"; then
  CPPFLAGS="${CPPFLAGS} -DDEBUG"
else
  CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi
  
# Fortran flags for French translation

if test "x$french" = "xyes"; then
  if test "x$cs_ibm_bg_type" != "x" ; then
    FCFLAGS="${FCFLAGS} -WF,-D_CS_LANG_FR"
  else
    FCFLAGS="${FCFLAGS} -D_CS_LANG_FR"
  fi
fi

# We may only turn on processing for libtool now that the basic compiler
# and linker flags are set (to avoid issues with linkers with different
# modes such as 32 and 64 bit which may be modified by compiler or
# linker flags).

# AC_DISABLE_SHARED may not be used inside a test, as it seems to
# take effect whether the result is true or not (maybe due to a bug),
# so we directly use autoconf's enable_shared variable.

if test "x$cs_ibm_bg_type" = "xL" -o "x$cs_ibm_bg_type" = "xP" \
     -o "$host_os" = cygwin -o "$host_os" = mingw32 ; then
   enable_shared=no
fi
AC_SUBST(enable_shared)

# Initialize libtool.

LT_INIT

# On some systems, such as Mac OS X in static mode, limitations
# of the linker require special handling for user user subroutines
# (due to not handling multiple definitions in the case of Mac OS X).
# Detect this here, so as to transfer the name of a specific solution
# (if necessary) to the Python package.

cs_special_user_link=
case "$host_os" in
  darwin*)
    if test "x$enable_shared" = "xno" ; then
      cs_special_user_link='ar_x'
    fi
    ;;
esac

AC_SUBST(cs_special_user_link)

#------------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------------

AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_PID_T

AC_HEADER_TIME
AC_HEADER_STDBOOL

AC_CHECK_TYPES([long long, unsigned long long])
AC_CHECK_TYPES([int32_t])
AC_CHECK_TYPES([int64_t])
AC_CHECK_TYPES([uint32_t])
AC_CHECK_TYPES([uint64_t])
AC_CHECK_TYPES([ptrdiff_t])

AC_CHECK_SIZEOF([void *])
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])

# Use long global numbers ?
AC_ARG_ENABLE(long-gnum,
  [AS_HELP_STRING([--enable-long-gnum],[use long global numbers])],
  [
    case "${enableval}" in
      yes) cs_have_long_gnum=yes ;;
      no)  cs_have_long_gnum=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-long-gnum]) ;;
    esac
  ],
  [ cs_have_long_gnum=no ]
)

if test "x$cs_have_long_gnum" = "xyes"; then
  AC_DEFINE([HAVE_LONG_GNUM],1,[Use 64-bit type if available for fvm_gnum_t.])
fi
AC_SUBST(cs_have_long_gnum)

# Check for Fortran 2003 flush feature
CS_AC_TEST_FC_FLUSH

# Check for Fortran module generation
CS_AC_TEST_FC_MOD

#------------------------------------------------------------------------------
# Determine OpenMP support
#------------------------------------------------------------------------------

# It seems that we may not use AC_OPENMP both for C and Fortran, so
# we use our own method here, based on flags already set by default
# or by the user.

# From this point, no runtime tests should be run, so no issues should
# arise due to missing RPATH or LD_LIBRARY_PATH values in case of non-standard
# compiler install paths, such as may happen with gcc's libgomp.

cs_have_openmp=no

AC_ARG_ENABLE(openmp,
  [AS_HELP_STRING([--disable-openmp], [disable OpenMP support])],
  [
    case "${enableval}" in
      yes) cs_have_openmp=yes ;;
      no)  cs_have_openmp=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
    esac
  ],
  [ cs_have_openmp=no ]
)

if test "x$cs_have_openmp" = "xyes" ; then

  saved_CFLAGS="$CFLAGS"
  saved_FCFLAGS="$FCFLAGS"
  saved_LDFLAGS="$LDFLAGS"

  CFLAGS="${CFLAGS} ${cflags_default_omp}"
  FCFLAGS="${FCFLAGS} ${fcflags_default_omp}"
  LDFLAGS="${LDFLAGS} ${cflags_default_omp}"

  AC_MSG_CHECKING([for OpenMP (C)])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <omp.h>]],
                 [[ omp_get_num_threads();]])],
                 [cs_have_openmp=yes],
                 [cs_have_openmp=no])
  AC_MSG_RESULT($cs_have_openmp)
  if test "x$cs_have_openmp" = "xyes" ; then
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([for OpenMP (Fortran)])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                   [[
      integer n
      n = omp_get_num_threads()
]])],
                   [cs_have_openmp=yes],
                   [cs_have_openmp=no])
    AC_MSG_RESULT($cs_have_openmp)
    AC_LANG_POP([Fortran])
  fi
  if test "x$cs_have_openmp" = "xyes" ; then
    AC_DEFINE([HAVE_OPENMP], 1, [openmp support])
  else
    cs_have_openmp=no
    CFLAGS="$saved_CFLAGS"
    FCFLAGS="$saved_FCFLAGS"
    LDFLAGS="$saved_LDFLAGS"
  fi
fi
AC_SUBST(cs_have_openmp)

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

CS_AC_TEST_PLE(["1.0.0"])

AM_CONDITIONAL(HAVE_INTERNAL_PLE, test x$cs_have_internal_ple = xyes)
if test "x$cs_have_internal_ple" = xyes; then
  AC_CONFIG_SUBDIRS([libple])
fi

CS_AC_TEST_LIBXML2
CS_AC_TEST_BLAS([$cs_have_openmp])
CS_AC_TEST_MPI
CS_AC_TEST_METIS
CS_AC_TEST_SCOTCH
CS_AC_TEST_HDF5
CS_AC_TEST_CGNS
CS_AC_TEST_MED
CS_AC_TEST_ADF
CS_AC_TEST_CCM
CS_AC_TEST_ZLIB
CS_AC_TEST_SYRTHES
CS_AC_TEST_SALOME
CS_AC_TEST_MEDCOUPLING

CS_AC_TEST_DOCS

cs_have_partitioner=no
if test x$cs_have_parmetis = xyes -o x$cs_have_metis = xyes ; then
  cs_have_partitioner=yes
fi
if test x$cs_have_ptscotch = xyes -o x$cs_have_scotch = xyes ; then
  cs_have_partitioner=yes
fi

AM_CONDITIONAL(HAVE_PARTITIONER, test x$cs_have_partitioner = xyes)

unset cs_have_partitioner

cs_have_libsyrcs=no
if test x$cs_have_syrthes = xyes -a x$cs_have_mpi = xyes ; then
  cs_have_libsyrcs=yes
fi

AM_CONDITIONAL(HAVE_LIBSYRCS, test x$cs_have_libsyrcs = xyes)

unset cs_have_libsyrcs

# Check if we need to force link with C++
cs_have_link_cxx=no
if test x$cs_have_med_link_cxx = xyes ; then
  cs_have_link_cxx=yes
fi
AM_CONDITIONAL(HAVE_LINK_CXX, test x$cs_have_link_cxx = xyes)

#------------------------------------------------------------------------------
# Enables Code_Saturne front-end 
#------------------------------------------------------------------------------

AC_ARG_ENABLE(frontend, [AS_HELP_STRING([--disable-frontend],
                                        [disable front-end elements])],
  [
    case "${enableval}" in
      yes) cs_have_frontend=yes ;;
      no)  cs_have_frontend=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-frontend]) ;;
    esac
  ],
  [ cs_have_frontend=yes ]
)
AM_CONDITIONAL(HAVE_FRONTEND, test x$cs_have_frontend = xyes)
AC_SUBST(cs_have_frontend)

#------------------------------------------------------------------------------
# Checks for Python support.
#------------------------------------------------------------------------------

AC_ARG_VAR([PYTHON], [the Python interpreter])
AM_PATH_PYTHON(["2.4"])

#------------------------------------------------------------------------------
# Enables Code_Saturne graphical user interface
#------------------------------------------------------------------------------

AC_ARG_ENABLE(gui, [AS_HELP_STRING([--disable-gui],
                                   [disable the Graphical User Interface])],
  [
    case "${enableval}" in
      yes) cs_have_gui=yes ;;
      no)  cs_have_gui=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-gui]) ;;
    esac
  ],
  [ cs_have_gui=yes ]
)

if test "x$cs_have_frontend" = "xno" ; then
  cs_have_gui=no
fi

#------------------------------------------------------------------------------
# Checks for PyQt support.
#------------------------------------------------------------------------------

# Test for PyQt4
if test "x$cs_have_gui" = "xyes" ; then

  AC_MSG_CHECKING([for PyQt4 version >= 4.3])

  prog="[import sys, string
try: import PyQt4
except ImportError: sys.exit(1)
from PyQt4.QtCore import *
if map(int, string.split(QT_VERSION_STR, '.')) < [4,3,0]: sys.exit(1)
if map(int, string.split(PYQT_VERSION_STR, '.')) < [4,3,0]: sys.exit(1)
sys.exit(0)]"

  ${PYTHON} -c "${prog}"
  retval=$?

  if test $retval -ne 0 ; then
    AC_MSG_RESULT([no])
    AC_MSG_FAILURE([cannot find PyQt4 support (>= 4.3), Graphical User Interface cannot be installed])
    cs_have_gui=no
  else
    AC_MSG_RESULT([yes])
  fi

fi

# Test for PyQt4 developper tools
if test "x$cs_have_gui" = "xyes" ; then

  # Try to find PyQt4 tools through Python interpreter
  prog="[import os, sys
from PyQt4.QtCore import QCoreApplication
app = QCoreApplication([])
path = app.applicationDirPath()
sys.stdout.write(path)]"

  PYPATH=`${PYTHON} -c "${prog}"`

  AC_ARG_VAR([PYUIC4], [PyQt4 user interfaces compiler])
  AC_ARG_VAR([PYRCC4], [PyQt4 resources compiler])

  AC_PATH_PROG([PYUIC4], [pyuic4], [], [${PYPATH}:$PATH])
  AC_PATH_PROG([PYRCC4], [pyrcc4], [], [${PYPATH}:$PATH])

  if test -z "$ac_cv_path_PYUIC4" -o -z "$ac_cv_path_PYRCC4"; then
    AC_MSG_FAILURE([cannot find PyQt4 dev tools, Graphical User Interface cannot be installed])
    cs_have_gui=no
  fi

fi

#------------------------------------------------------------------------------
# Checks for Qt tools.
#------------------------------------------------------------------------------

# Test for Qt4 developper tools
if test "x$cs_have_gui" = "xyes" ; then

  AC_ARG_VAR([LRELEASE], [Qt message files generator])
  AC_ARG_VAR([MOC], [Qt meta objects compiler])

  AC_PATH_PROGS([LRELEASE], [lrelease lrelease-qt4])
  AC_PATH_PROGS([MOC], [moc moc-qt4])

  if test "$LRELEASE" = : -o "$MOC" = : ; then
    AC_MSG_FAILURE([cannot find Qt4 dev tools, Graphical User Interface cannot be installed])
    cs_have_gui=no
  fi

fi

AM_CONDITIONAL(HAVE_GUI, [test "x$cs_have_gui" = "xyes"])
AC_SUBST(cs_have_gui)

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([sys/types.h sys/utsname.h sys/stat.h dirent.h stddef.h])
AC_CHECK_HEADERS([malloc.h unistd.h fcntl.h sys/types.h sys/signal.h])
AC_CHECK_HEADERS([sys/procfs.h sys/sysinfo.h sys/resource.h])
AC_CHECK_HEADERS([float.h string.h sys/time.h])

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([mkdir chdir stat dup2])
AC_CHECK_FUNCS([snprintf])
AC_CHECK_FUNCS([getcwd sleep])
AC_CHECK_FUNCS([getpwuid geteuid])
AC_CHECK_FUNCS([uname])
AC_CHECK_FUNCS([clock_gettime getrusage gettimeofday sbrk sysinfo])
AC_CHECK_FUNCS([posix_memalign])
AC_CHECK_FUNCS([memset pow modf])

AC_FUNC_STRTOD

#------------------------------------------------------------------------------
# Checks for system services.
#------------------------------------------------------------------------------

# Largefile support (may be disabled in case of conflicts)
#------------------

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Test for size of off_t once largefile support is activated

AC_CHECK_SIZEOF([off_t])

# Optional malloc hooks support (for glibc)
#------------------------------

cs_have_malloc_hooks=no

AC_ARG_ENABLE(malloc_hooks,
  [AS_HELP_STRING([--enable-malloc-hooks], [use malloc hooks when available])],
  [
    case "${enableval}" in
      yes) malloc_hooks=yes ;;
      no)  malloc_hooks=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-malloc-hooks]) ;;
    esac
  ],
  [ malloc_hooks=no ]
)

# If we use malloc hooks, we may need to add a -lmalloc link flag.

if test "x$malloc_hooks" = "xyes" ; then

  # Basic test
  AC_MSG_CHECKING([for malloc hooks])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],
                 [[ int __malloc_hook(); __malloc_hook(); ]])],
                 [cs_have_malloc_hooks=yes],
                 [cs_have_malloc_hooks=no])
  AC_MSG_RESULT($cs_have_malloc_hooks)

fi

if test "x$cs_have_malloc_hooks" = "xyes"; then
  AC_DEFINE(HAVE_MALLOC_HOOKS, 1, HAVE_MALLOC_HOOKS)
fi

# Backtrace support (for glibc)
#------------------

# Basic backtrace functionnality

cs_have_glibc_backtrace=no

AC_MSG_CHECKING([for glibc backtrace])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#if defined(__GNUC__)
#define _GNU_SOURCE
#include <memory.h>
#include <execinfo.h>
#endif]],
               [[ void *tr_array[10]; backtrace(tr_array, 10); ]])],
               [cs_have_glibc_backtrace=yes],
               [cs_have_glibc_backtrace=no])
AC_MSG_RESULT($cs_have_glibc_backtrace)

if test "x$cs_have_glibc_backtrace" = "xyes"; then
  AC_DEFINE(HAVE_GLIBC_BACKTRACE, 1, HAVE_GLIBC_BACKTRACE)
fi

# C++ demangling

cs_have_cplus_demangle=no

AC_MSG_CHECKING([for glibc cplus_demangle])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <demangle.h>]],
               [[ const char *s; cplus_demangle(s, auto_demangling); ]])],
               [cs_have_cplus_demangle=yes],
               [cs_have_cplus_demangle=no])
AC_MSG_RESULT($cs_have_cplus_demangle)

if test "x$cs_have_cplus_demangle" = "xyes"; then
  AC_DEFINE(HAVE_CPLUS_DEMANGLE, 1, HAVE_CPLUS_DEMANGLE)
fi

# Optional shared library support
#--------------------------------

cs_have_dlloader=no

AC_ARG_ENABLE(dlloader,
  [AS_HELP_STRING([--disable-dlloader], [disable dynamic shared library loading])],
  [
    case "${enableval}" in
      yes) cs_have_dlloader=yes ;;
      no)  cs_have_dlloader=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-dlloader]) ;;
    esac
  ],
  [ cs_have_dlloader=yes ]
)

if test "x$cs_have_dlloader" = "xyes" ; then

  saved_LIBS="$LIBS"
  LIBS="$LIBS -ldl"

  AC_MSG_CHECKING([for dlopen])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <dlfcn.h>]],
                 [[ dlopen("none.so", RTLD_LAZY);]])],
                 [cs_have_dlloader=dlopen],
                 [cs_have_dlloader=no])
  AC_MSG_RESULT($cs_have_dlloader)
  if test "x$cs_have_dlloader" = "xdlopen" ; then
    AC_DEFINE([HAVE_DLOPEN], 1, [dlopen support])
  else
    LIBS="$saved_LIBS"
    cs_have_dlloader=no
  fi

fi

# Plugin modules support
#-----------------------

cs_have_plugins=no
if test "x$cs_have_dlloader" = "xdlopen" ; then
  # Test for different plugins here
  if test x$cs_have_medcoupling = xyes ; then
    cs_have_plugins=yes
  fi
fi

AM_CONDITIONAL(HAVE_PLUGINS, test x$cs_have_plugins = xyes)

#------------------------------------------------------------------------------
# Internationalization
#------------------------------------------------------------------------------

AM_GNU_GETTEXT([external])

#------------------------------------------------------------------------------
# Pass compiler options to automake files and cs_config.py
#------------------------------------------------------------------------------

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(FCFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------

# Main configuration

AC_CONFIG_FILES([Makefile po/Makefile.in
                 preprocessor/Makefile
                 preprocessor/appli/Makefile preprocessor/base/Makefile
                 preprocessor/pre-post/Makefile preprocessor/util/Makefile
                 preprocessor/tools/Makefile
                 salome/Makefile
                 src/Makefile      src/apps/Makefile src/gui/Makefile
                 src/bft/Makefile  src/fvm/Makefile  src/mei/Makefile
                 src/base/Makefile src/atmo/Makefile src/cfbl/Makefile
                 src/cogz/Makefile src/cplv/Makefile src/ctwr/Makefile
                 src/elec/Makefile src/comb/Makefile src/pprt/Makefile
                 src/lagr/Makefile src/rayt/Makefile src/turb/Makefile
                 src/alge/Makefile src/mesh/Makefile src/user/Makefile
                 src/user_examples/Makefile libsyrcs/Makefile
                 gui/Makefile gui/Base/Makefile gui/Pages/Makefile
                 salome/fsi_coupling/Makefile
                 salome/fsi_coupling/idl/Makefile
                 salome/fsi_coupling/libmilieu/Makefile
                 salome/fsi_coupling/resources/Makefile
                 salome/fsi_coupling/src/Makefile
                 salome/fsi_coupling/src/FSI_ASTER/Makefile
                 salome/fsi_coupling/src/FSI_MILIEU/Makefile
                 salome/fsi_coupling/src/FSI_SATURNE/Makefile
                 salome/cfd_study/Makefile
                 salome/cfd_study/doc/Makefile
                 salome/cfd_study/idl/Makefile
                 salome/cfd_study/resources/Makefile
                 salome/cfd_study/src/Makefile
                 salome/cfd_study/src/CFDSTUDY/Makefile
                 salome/cfd_study/src/CFDSTUDYGUI/Makefile
                 docs/Makefile docs/doxygen/Makefile docs/style/csvers.tex
                 docs/autovnv/Makefile docs/autovnv/graphics/Makefile
                 docs/developper/Makefile docs/developper/graphics/Makefile
                 docs/install/Makefile
                 docs/refcard/Makefile docs/refcard/graphics/Makefile
                 docs/theory/Makefile docs/theory/graphics/Makefile
                 docs/tutorial/Makefile docs/tutorial/graphics/Makefile
                 docs/user/Makefile docs/user/graphics/Makefile
                 tests/Makefile
                 extras/salome/fsi_appli_config.xml
                 bin/cs_config.py bin/SaturneGUI])

#------------------------------------------------------------------------------

AC_OUTPUT

#------------------------------------------------------------------------------

# Libtool configuration bug workaround (problem encountered with compilers
# which do not recognize "-soname", so that "-Wl,-soname" is needed to pass
# flags to the linker; Manually setting $wl to "-Wl" in libtool corrects
# the problem, and so does this workaround borrowed from
# the HDF5 configure.in, which does it automatically).

if test "x$GCC" != "xyes"; then
  echo 'fixing $wl in libtool'
ed - libtool <<EOF 2> /dev/null
/^wl=""/s//wl="-Wl,"/
w
q
EOF
fi

#------------------------------------------------------------------------------

echo
echo "Configuration options:"
echo " use debugging code: "$debug""
echo " use malloc hooks: "$cs_have_malloc_hooks""
if test x$cs_have_frontend = xno ; then
  echo " use frontend: no"
else
  echo " use graphical user interface: "$cs_have_gui""
fi
echo " use long integers: "$cs_have_long_gnum""
echo " Zlib (gzipped file) support: "$cs_have_zlib""
echo " MPI (Message Passing Interface) support: "$cs_have_mpi""
if (test x$cs_have_mpi = xyes) ; then
  echo "   MPI I/O support: "$cs_have_mpi_io""
  echo "   MPI2 one-sided communication support: "$cs_have_mpi_one_sided""
fi
echo " OpenMP support: "$cs_have_openmp""
echo " BLAS (Basic Linear Algebra Subprograms) support: "$cs_have_blas""
echo " Libxml2 (XML Reader) support: "$cs_have_libxml2""
echo " ParMETIS (Parallel Graph Partitioning) support: "$cs_have_parmetis""
if test x$cs_have_parmetis = xno ; then
  echo " METIS (Graph Partitioning) support: "$cs_have_metis""
fi
echo " PT-SCOTCH (Parallel Graph Partitioning) support: "$cs_have_ptscotch""
if test x$cs_have_ptscotch = xno ; then
  echo " SCOTCH (Graph Partitioning) support: "$cs_have_scotch""
fi
echo " ADF support: "$cs_have_adf""
echo " CCM support: "$cs_have_ccm""
echo " HDF (Hierarchical Data Format) support: "$cs_have_hdf5""
echo " CGNS (CFD General Notation System) support: "$cs_have_cgns""
echo " MED (Model for Exchange of Data) support: "$cs_have_med""
if test x$cs_have_med = xyes ; then
  echo "   MED MPI I/O support: "$cs_have_med_mpi""
fi
echo " MEDCoupling support: "$cs_have_medcoupling""
if test x$cs_have_medcoupling = xyes ; then
  echo " ParaMEDMEM support: "$cs_have_paramedmem""
fi
echo " SALOME GUI support: "$cs_have_salome_gui""
echo " SALOME Kernel support: "$cs_have_salome_kernel""
echo " Dynamic loader support (for YACS): "$cs_have_dlloader""
echo " SYRTHES 3 coupling support: "$cs_have_syrthes""
echo
if test x$cs_have_med2_headers = xyes ; then
  echo " Warning: MED version used is 2.3, which is deprecated."
  echo "          in future versions, MED >= 3.0 will be required."
  echo
fi
if test x$cs_have_cgns2_headers = xyes ; then
  echo " Warning: CGNS version used is 2.5, which is deprecated."
  echo "          in future versions, CGNS >= 3.1 will be required."
  echo
fi
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the Code_Saturne kernel"
echo


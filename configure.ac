# Process this file with autoconf to produce a configure script.

#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009-2010 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([cs_licence_c_comment],
[/*
   This file is part of the Code_Saturne Kernel, element of the
   Code_Saturne CFD tool.

   Copyright (C) 2009-2010 EDF S.A., France

   The Code_Saturne Kernel is free software; you can redistribute it
   and/or modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2 of
   the License, or (at your option) any later version.

   The Code_Saturne Kernel is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public Licence
   along with the Code_Saturne Preprocessor; if not, write to the
   Free Software Foundation, Inc.,
   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/])

m4_define([cs_major_version], [2])
m4_define([cs_minor_version], [0])
m4_define([cs_release_version], [0])
m4_define([cs_version_extra], [-rc1])
m4_define([cs_version_string],
          [cs_major_version.cs_minor_version.cs_release_version@&t@cs_version_extra])

CS_GET_VERSION

AC_INIT([code_saturne],[cs_version_string],[saturne-support@edf.fr])
AC_CONFIG_SRCDIR([src/apps/cs_solver.c])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_SUBST(ac_aux_dir)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([cs_config.h])

AM_GNU_GETTEXT_VERSION([0.17])

AM_INIT_AUTOMAKE([-Wall])

AM_MAINTAINER_MODE

# Enable silent rules with "./configure --enable-silent-rules" or "make V=0"
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

# Remember the configure command line
CONFIG_CMD=`echo "$0 $ac_configure_args"`

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

# Check for C and Fortran compilers
#----------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_FCFLAGS=$FCFLAGS
user_LDFLAGS=$LDFLAGS
user_LIBS=$LIBS

AC_PROG_CC
#if test "x$ac_cv_prog_ac_ct_CC" = x; then
#  AC_MSG_ERROR([no C compiler found, please install one])
#fi

AC_PROG_FC
#if test "x$ac_cv_prog_ac_ct_FC" = x; then
#  AC_MSG_ERROR([no Fortran 90 compiler found, please install one])
#fi

AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_CC_C_O
AC_PROG_SED

AC_FC_LIBRARY_LDFLAGS

# Debug or production compilation mode (debug by default) ?

AC_ARG_ENABLE(debug,
  [AS_HELP_STRING([--enable-debug], [enable debugging (reduces optimization)])],
  [
    case "${enableval}" in
      yes) debug=true ;;
      no)  debug=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac
  ],
  [ debug=false ]
)

# Optionally deactivate automatic determination of flags on known systems

AC_ARG_ENABLE(auto-flags,
  [AS_HELP_STRING([--disable-auto-flags], [do not define *FLAGS on known systems])],
  [
    case "${enableval}" in
      yes) auto_flags=true ;;
      no)  auto_flags=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-auto-flags]) ;;
    esac
  ],
  [ auto_flags=true ]
)

# Optionally activate French translation of Fortran logs (compile-time)

AC_ARG_ENABLE(french,
  [AS_HELP_STRING([--enable-french], [enable French translation of Fortran logs])],
  [
    case "${enableval}" in
      yes) french=true ;;
      no)  french=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-french]) ;;
    esac
  ],
  [ french=false ]
)

# Default compiler options (may be modified
# by defining CFLAGS in the environment)
#------------------------------------------

if test "x$auto_flags" = "xtrue" ; then

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
  FCFLAGS="$fcflags_default $user_FCFLAGS"
  LDFLAGS="$ldflags_default $user_LDFLAGS"
  LIBS="$libs_default $user_LIBS"
  LDRPATH="$ldflags_rpath"

  if test "x$debug" = xtrue; then
    # Debug flags
    CFLAGS_DBG="$cflags_default_dbg"
    FCFLAGS_DBG="$fcflags_default_dbg"
    # Add debug flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_dbg"
  else
    # Normal optimization flags
    CFLAGS_OPT="$cflags_default_opt"
    FCFLAGS_OPT="$fcflags_default_opt"
    # Hot optimization flags
    CFLAGS_HOT="$cflags_default_hot"
    FCFLAGS_HOT="$fcflags_default_hot"
    # Add optimization flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_opt"
  fi
  CFLAGS_EXT="$cflags_default_ext"

  AC_SUBST(CFLAGS_DBG)
  AC_SUBST(CFLAGS_OPT)
  AC_SUBST(CFLAGS_HOT)
  AC_SUBST(CFLAGS_EXT)

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

if test "x$debug" = "xtrue"; then
  CPPFLAGS="${CPPFLAGS} -DDEBUG"
else
  CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi
  
# Fortran flags for French translation

if test "x$french" = "xtrue"; then
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

if test "x$cs_ibm_bg_type" != "x" ; then
   enable_shared=no
fi

LT_INIT

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
      yes) have_long_gnum=yes ;;
      no)  have_long_gnum=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-long-gnum]) ;;
    esac
  ],
  [ have_long_gnum=no ]
)

if test "x$have_long_gnum" = "xyes"; then
  AC_DEFINE([HAVE_LONG_GNUM],1,[Use 64-bit type if available for fvm_gnum_t.])
fi

# Use long int for mesh reading ?
AC_ARG_ENABLE(long-int,
  [AS_HELP_STRING([--enable-long-int], [use long integers for connectivity])],
  [
    case "${enableval}" in
      yes) long_int=yes ;;
      no)  long_int=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-long-int]) ;;
    esac
  ],
  [ long_int=no ]
)

if test "x$long_int" = "xyes"; then
  AC_DEFINE([USE_LONG_INT], [1], [Use long int for connectivity])
fi

# Check for Fortran 2003 flush feature
CS_AC_TEST_FC_FLUSH

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

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

CS_AC_TEST_PLE(["1.0.0"])

AM_CONDITIONAL(HAVE_INTERNAL_PLE, test x$have_internal_ple = xyes)
if test "x$have_internal_ple" = xyes; then
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

have_partitioner=no
if test x$have_parmetis = xyes -o x$have_metis = xyes ; then
  have_partitioner=yes
fi
if test x$have_ptscotch = xyes -o x$have_scotch = xyes ; then
  have_partitioner=yes
fi

AM_CONDITIONAL(HAVE_PARTITIONER, test x$have_partitioner = xyes)

unset have_partitioner

#------------------------------------------------------------------------------
# Enables Code_Saturne preprocessor
#------------------------------------------------------------------------------

AC_ARG_ENABLE(prepro, [AS_HELP_STRING([--disable-prepro],
                                      [disable the Preprocessor])],
  [
    case "${enableval}" in
      yes) have_prepro=yes ;;
      no)  have_prepro=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-prepro]) ;;
    esac
  ],
  [ have_prepro=yes ]
)
AM_CONDITIONAL(HAVE_PREPROCESSOR, test x$have_prepro = xyes)

#------------------------------------------------------------------------------
# Enables Code_Saturne graphical user interface
#------------------------------------------------------------------------------

AC_ARG_ENABLE(gui, [AS_HELP_STRING([--disable-gui],
                                   [disable the Graphical User Interface])],
  [
    case "${enableval}" in
      yes) have_gui=yes ;;
      no)  have_gui=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-gui]) ;;
    esac
  ],
  [ have_gui=yes ]
)

#------------------------------------------------------------------------------
# Checks for Python support.
#------------------------------------------------------------------------------

AC_ARG_WITH(python, [AS_HELP_STRING([--with-python=PYTHON],
                                    [specify a Python interpreter])],
                         [cs_python=${withval}],
                         [cs_python=python])dnl

cs_python_basename=`basename ${cs_python}`
if test "x${cs_python_basename}" != "x${cs_python}" ; then
  cs_python_dirname=`dirname ${cs_python}`
  AC_PATH_PROG([PYTHON], [${cs_python_basename}], ["no"], [${cs_python_dirname}])
  unset cs_python_dirname
else
  AC_PATH_PROG([PYTHON], [${cs_python}], ["no"])
fi
unset cs_python_basename

cs_python=${PYTHON}
AC_SUBST(cs_python)

if test "x$have_gui" = xyes; then
  cs_python_version="2.4"
else
  cs_python_version="2.3"
fi

AM_PATH_PYTHON([$cs_python_version])

#------------------------------------------------------------------------------
# Checks for PyQt support.
#------------------------------------------------------------------------------

AC_ARG_WITH(pyqt4-exec, [AS_HELP_STRING([--with-pyqt4-exec=PATH],
                                        [specify prefix directory for PyQt4 developper tools (pyuic4 and pyrcc4)])],
                        [PYQT4_EXE=${withval}],
                        [PYQT4_EXE=none])dnl

# Test for PyQt4
if test "x$have_gui" = "xyes" ; then

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
    have_gui=no
  else
    AC_MSG_RESULT([yes])
  fi

fi

# Test for PyQt4 developper tools
if test "x$have_gui" = "xyes" ; then

  if test "x${PYQT4_EXE}" = "xnone" ; then

    # Try to find PyQt4 tools through Python interpreter
    prog="[import os, sys
from PyQt4.QtCore import QCoreApplication
app = QCoreApplication([])
path = app.applicationDirPath()
sys.stdout.write(path)]"

    PYPATH=`${PYTHON} -c "${prog}"`

  else
    PYPATH=${PYQT4_EXE}
  fi

  AC_PATH_PROG([PYUIC4], [pyuic4], ["no"], [${PYPATH}:$PATH])
  AC_PATH_PROG([PYRCC4], [pyrcc4], ["no"], [${PYPATH}:$PATH])

  if test "x$PYUIC4" = "xno" -o "x$PYRCC4" = "xno" ; then
    AC_MSG_FAILURE([cannot find PyQt4 dev tools, Graphical User Interface cannot be installed])
    have_gui=no
  fi

fi

AM_CONDITIONAL(HAVE_GUI, [test "x$have_gui" = "xyes"])

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([sys/types.h sys/utsname.h sys/stat.h stddef.h])
AC_CHECK_HEADERS([malloc.h unistd.h fcntl.h sys/types.h sys/signal.h])
AC_CHECK_HEADERS([sys/procfs.h sys/sysinfo.h sys/resource.h])

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([mkdir chdir stat dup2])
AC_CHECK_FUNCS([snprintf])
AC_CHECK_FUNCS([getcwd sleep])
AC_CHECK_FUNCS([getpwuid geteuid])
AC_CHECK_FUNCS([uname])
AC_CHECK_FUNCS([getrusage gettimeofday sbrk sysinfo])
AC_CHECK_FUNCS([posix_memalign])
AC_CHECK_FUNCS([memset pow])

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

# Optional socket support
#------------------------

have_socket=no

AC_ARG_ENABLE(sockets,
  [AS_HELP_STRING([--disable-sockets], [disable communications through IP sockets])],
  [
    case "${enableval}" in
      yes) have_socket=yes ;;
      no)  have_socket=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-sockets]) ;;
    esac
  ],
  [ have_socket=yes ]
)
if test "x$have_socket" = "xyes" ; then
  AC_DEFINE([HAVE_SOCKET], 1, [sockets support])
fi

# Optional shared library support
#--------------------------------

have_dlloader=no

AC_ARG_ENABLE(dlloader,
  [AS_HELP_STRING([--disable-dlloader], [disable dynamic shared library loading])],
  [
    case "${enableval}" in
      yes) have_dlloader=yes ;;
      no)  have_dlloader=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-dlloader]) ;;
    esac
  ],
  [ have_dlloader=yes ]
)

if test "x$have_dlloader" = "xyes" ; then

  saved_LIBS="$LIBS"
  LIBS="$LIBS -ldl"

  AC_MSG_CHECKING([for dlopen])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <dlfcn.h>]],
                 [[ dlopen("none.so", RTLD_LAZY);]])],
                 [have_dlloader=yes],
                 [have_dlloader=no])
  AC_MSG_RESULT($have_dlloader)
  if test "x$have_dlloader" = "xyes" ; then
    AC_DEFINE([HAVE_DLOPEN], 1, [dlopen support])
  else
    LIBS="$saved_LIBS"
    have_dlloader=no
  fi

fi

# Optional malloc hooks support (for glibc)
#------------------------------

have_malloc_hooks=no

AC_ARG_ENABLE(malloc_hooks,
  [AS_HELP_STRING([--disable-malloc-hooks], [Do not use malloc hooks when available])],
  [
    case "${enableval}" in
      yes) malloc_hooks=true ;;
      no)  malloc_hooks=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-malloc-hooks]) ;;
    esac
  ],
  [ malloc_hooks=true ]
)

# If we use malloc hooks, we may need to add a -lmalloc link flag.

if test "x$malloc_hooks" = "xtrue" ; then

  # Basic test
  AC_MSG_CHECKING([for malloc hooks])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[]],
                 [[ int __malloc_hook(); __malloc_hook(); ]])],
                 [have_malloc_hooks=yes],
                 [have_malloc_hooks=no])
  AC_MSG_RESULT($have_malloc_hooks)

fi

if test "x$have_malloc_hooks" = "xyes"; then
  AC_DEFINE(HAVE_MALLOC_HOOKS, 1, HAVE_MALLOC_HOOKS)
fi

# Backtrace support (for glibc)
#------------------

# Basic backtrace functionnality

have_glibc_backtrace=no

AC_MSG_CHECKING([for glibc backtrace])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#if defined(__GNUC__)
#define _GNU_SOURCE
#include <memory.h>
#include <execinfo.h>
#endif]],
               [[ void *tr_array[10]; backtrace(tr_array, 10); ]])],
               [have_glibc_backtrace=yes],
               [have_glibc_backtrace=no])
AC_MSG_RESULT($have_glibc_backtrace)

if test "x$have_glibc_backtrace" = "xyes"; then
  AC_DEFINE(HAVE_GLIBC_BACKTRACE, 1, HAVE_GLIBC_BACKTRACE)
fi

# C++ demangling

have_cplus_demangle=no

AC_MSG_CHECKING([for glibc cplus_demangle])
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <demangle.h>]],
               [[ const char *s; cplus_demangle(s, auto_demangling); ]])],
               [have_cplus_demangle=yes],
               [have_cplus_demangle=no])
AC_MSG_RESULT($have_cplus_demangle)

if test "x$have_cplus_demangle" = "xyes"; then
  AC_DEFINE(HAVE_CPLUS_DEMANGLE, 1, HAVE_CPLUS_DEMANGLE)
fi

#------------------------------------------------------------------------------
# Internationalization
#------------------------------------------------------------------------------

AC_DEFINE_UNQUOTED(LOCALEDIR, "${localedir}", LOCALEDIR)
AM_GNU_GETTEXT([external])

#------------------------------------------------------------------------------
# Pass compiler options to automake files and cs_config.py
#------------------------------------------------------------------------------

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(FCFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------

# Output configuration info

CS_AC_CONFIG_INFO_INIT([cs_config_comments_include])

CS_AC_CONFIG_INFO_EXTRA([
# use debugging code:                              "$debug"
# use malloc hooks:                                "$have_malloc_hooks"
# use graphical user interface:                    "$have_gui"
# MPI (Message Passing Interface) support:         "$cs_have_mpi"
# OpenMP support:                                  "$cs_have_openmp"
# BLAS (Basic Linear Algebra Subprograms) support: "$cs_have_blas"
# Libxml2 (XML Reader) support:                    "$have_libxml2"
# Zlib (gzipped file) support:                     "$have_zlib"
# ADF support:                                     "$have_adf"
# CCM support:                                     "$have_ccm"
# HDF5 (Hierarchical Data Format) support:         "$have_hdf5"
# CGNS (CFD General Notation System) support:      "$have_cgns"
# MED (Model for Exchange of Data) support:        "$have_med"
# IP socket support (for SYRTHES 3 or CFD_Proxy):  "$have_socket"
# Dynamic loader support (for YACS):               "$have_dlloader"
])

CS_AC_CONFIG_INFO_CC([${cs_ac_cc_version}], [${cs_ac_cc_version_full}])
CS_AC_CONFIG_INFO_FC([${cs_ac_fc_version}], [${cs_ac_fc_version_full}])

CS_AC_CONFIG_INFO_FINALIZE

AC_SUBST_FILE([cs_config_comments])
cs_config_comments=cs_config_comments_include

# Main configuration

AC_CONFIG_FILES([Makefile po/Makefile.in
                 preprocessor/Makefile
                 preprocessor/appli/Makefile preprocessor/base/Makefile
                 preprocessor/pre-post/Makefile preprocessor/util/Makefile
                 preprocessor/tools/Makefile
                 src/Makefile      src/apps/Makefile
                 src/bft/Makefile  src/fvm/Makefile  src/mei/Makefile
                 src/base/Makefile src/atmo/Makefile src/cfbl/Makefile
                 src/cogz/Makefile src/cplv/Makefile src/ctwr/Makefile
                 src/elec/Makefile src/fuel/Makefile src/mati/Makefile
                 src/pprt/Makefile src/lagr/Makefile src/rayt/Makefile
                 libsyrcs/Makefile
                 gui/Makefile gui/Base/Makefile gui/Pages/Makefile
                 doc/Makefile doc/style/csvers.tex
                 doc/howto/Makefile doc/howto/graphics/Makefile
                 doc/developper/Makefile doc/developper/graphics/Makefile
                 doc/install/Makefile
                 doc/refcard/Makefile doc/refcard/graphics/Makefile
                 doc/theory/Makefile doc/theory/graphics/Makefile
                 doc/tutorial/Makefile doc/tutorial/graphics/Makefile
                 doc/user/Makefile doc/user/graphics/Makefile
                 tests/Makefile
                 bin/cs_config_build.py bin/SaturneGUI])

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

rm -f cs_config_comments_include

echo
echo "Configuration options:"
echo " use debugging code: "$debug""
echo " use malloc hooks: "$have_malloc_hooks""
echo " use graphical user interface: "$have_gui""
echo " use preprocessor: "$have_prepro""
echo " use long integers: "$long_int""
echo " Zlib (gzipped file) support: "$have_zlib""
echo " MPI (Message Passing Interface) support: "$cs_have_mpi""
if (test x$cs_have_mpi = xyes) ; then
  echo "   MPI I/O support: "$cs_have_mpi_io""
  echo "   MPI2 one-sided communication support: "$cs_have_mpi_one_sided""
fi
echo " OpenMP support: "$cs_have_openmp""
echo " BLAS (Basic Linear Algebra Subprograms) support: "$cs_have_blas""
echo " Libxml2 (XML Reader) support: "$have_libxml2""
echo " ParMETIS (Parallel Graph Partitioning) support: "$have_parmetis""
if test x$have_parmetis = xno ; then
  echo " METIS (Graph Partitioning) support: "$have_metis""
fi
echo " PT-SCOTCH (Parallel Graph Partitioning) support: "$have_ptscotch""
if test x$have_ptscotch = xno ; then
  echo " SCOTCH (Graph Partitioning) support: "$have_scotch""
fi
echo " ADF support: "$have_adf""
echo " CCM support: "$have_ccm""
echo " HDF (Hierarchical Data Format) support: "$have_hdf5""
echo " CGNS (CFD General Notation System) support: "$have_cgns""
echo " MED (Model for Exchange of Data) support: "$have_med""
echo " SYRTHES 3 coupling support: "$have_syrthes""
echo " IP socket support (for SYRTHES 3 or CFD_Proxy): "$have_socket""
echo " Dynamic loader support (for YACS): "$have_dlloader""
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the Code_Saturne kernel"
echo


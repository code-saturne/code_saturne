# Process this file with autoconf to produce a configure script

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

  Copyright (C) 1998-2018 EDF S.A.

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

# Note: proper determination of the version string requires the NEWS file
#       to be maintained properly (as it should be anyways).

m4_define([cs_version_string],
          [m4_normalize(esyscmd([build-aux/cs_version.py]))])

AC_INIT([code_saturne],[cs_version_string],[saturne-support@edf.fr],[],[http://code-saturne.org])
AC_CONFIG_SRCDIR([src/apps/cs_solver.c])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([cs_config.h])

AM_GNU_GETTEXT_VERSION([0.19])

# Initialize automake with the following options:
# - foreign, so as to handle the absence of ChangeLog (automatically generated)
# - tar-pax, so as to handle long lines (> 99 characters) in tar archives
# Warnings can be activated at bootstrap with 'autoreconf -vi --warnings=all'
AM_INIT_AUTOMAKE([foreign tar-pax subdir-objects])

# Enable maintainer mode by default for a developer checkout
AS_IF([test -d ${srcdir}/.git],
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

CS_AC_SALOME_ENV

# Check for C and Fortran compilers
#----------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_CFLAGS_HOT=$CFLAGS_HOT
user_CXXFLAGS=$CXXFLAGS
user_FCFLAGS=$FCFLAGS
user_FCFLAGS_HOT=$FCFLAGS_HOT
user_LDFLAGS=$LDFLAGS
user_LIBS=$LIBS
user_CS_LD=$CS_LD

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC

AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_CC_C_O
AC_PROG_SED

if test "x$CC" = "x" ; then
  AC_MSG_FAILURE([cannot find C compiler])
fi
if test "x$FC" = "x" ; then
  AC_MSG_FAILURE([cannot find Fortran compiler])
fi

# Needed to use Microsoft archiver lib.exe
# It copies the ar-lib script, similar to the compile script
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

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

AC_ARG_ENABLE(profile,
  [AS_HELP_STRING([--enable-profile], [enable profiling])],
  [
    case "${enableval}" in
      yes) profile=yes ;;
      no)  profile=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-profile]) ;;
    esac
  ],
  [ profile=no ]
)
AC_SUBST(profile)

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
if test "${relocatable}" = yes ; then
  AC_DEFINE([HAVE_RELOCATABLE], 1, [Relocatable installation])
fi

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
  CXXFLAGS_STD="$cxxflags_default_std"
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
    CFLAGS_HOT="$cflags_default_hot $user_CFLAGS_HOT"
    CXXFLAGS_HOT="$cxxflags_default_hot"
    FCFLAGS_HOT="$fcflags_default_hot $user_FCFLAGS_HOT"
    # Add optimization flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_opt"
    if test "x$profile" = xyes; then
      # Profiling flags
      CFLAGS_DBG="$cflags_default_prf"
      CXXFLAGS_DBG="$cxxflags_default_prf"
      FCFLAGS_DBG="$fcflags_default_prf"
      # Add profiling flags for linker
      LDFLAGS="$LDFLAGS $ldflags_default_prf"
    fi
  fi

  AC_SUBST(CFLAGS_DBG)
  AC_SUBST(CFLAGS_OPT)
  AC_SUBST(CFLAGS_HOT)

  AC_SUBST(CXXFLAGS_DBG)
  AC_SUBST(CXXFLAGS_OPT)
  AC_SUBST(CXXFLAGS_HOT)
  AC_SUBST(CXXFLAGS_STD)

  AC_SUBST(FCFLAGS_DBG)
  AC_SUBST(FCFLAGS_OPT)
  AC_SUBST(FCFLAGS_HOT)

  AC_SUBST(LDRPATH)
  AC_SUBST(ccversion, [${cs_cc_version}])

fi

# Add some preprocessor flags to include additional system features
#------------------------------------------------------------------

case "$host_os" in
  linux*)
    CPPFLAGS="${CPPFLAGS} -D_POSIX_SOURCE"
    ;;
  darwin*)
    CPPFLAGS="${CPPFLAGS} -D_DARWIN_C_SOURCE"
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

if test "$host_os" = mingw32 ; then
   enable_shared=no
fi
AC_SUBST(enable_shared)

# Initialize libtool.

LT_INIT([disable-static, win32-dll])

# In static mode, multiple definitions must be allowed for the linker.
# On some systems, such as Mac OS X in static mode, limitations
# of the linker require special handling for user user subroutines
# (due to not handling multiple definitions in the case of Mac OS X).
# Detect this here, so as to transfer the name of a specific solution
# (if necessary) to the Python package.

cs_special_user_link=
if test "x$enable_shared" = "xno" ; then
  case "$host_os" in
    linux* | *bsd* | cnk* | mingw32 | cygwin* )
      if test "x$auto_flags" = "xyes" ; then
        LDFLAGS="${LDFLAGS} -Wl,--allow-multiple-definition"
      fi
      ;;
    darwin*)
      cs_special_user_link='ar_x'
      ;;
  esac
else
  case "$host_os" in
    linux* )
      if test "x$auto_flags" = "xyes" ; then
        case "$cs_ac_cc_version" in
          Cray* )
            CFLAGS="${CFLAGS} -fPIC"
            ;;
          *)
            LDFLAGS="${LDFLAGS} -Wl,-export-dynamic"
            ;;
        esac
        case "$cs_ac_cxx_version" in
          Cray* )
            CXXFLAGS="${CXXFLAGS} -fPIC"
            ;;
          *)
            ;;
        esac
        case "$cs_ac_fc_version" in
          Cray* )
            FCFLAGS="${FCFLAGS} -fPIC"
            ;;
          *)
            ;;
        esac
      fi
      ;;
  esac
fi
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
  [ cs_have_long_gnum=yes ]
)

if test "x$cs_have_long_gnum" = "xyes"; then
  AC_DEFINE([HAVE_LONG_GNUM],1,[Use 64-bit type if available for cs_gnum_t.])
fi
AC_SUBST(cs_have_long_gnum)

# Check for Fortran module generation
CS_AC_TEST_FC_MOD

# First determination of Fortran libraries (before OpenMP test)
AC_FC_LIBRARY_LDFLAGS

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
cs_have_openmp_f=no

AC_ARG_ENABLE(openmp,
  [AS_HELP_STRING([--disable-openmp], [disable OpenMP support])],
  [
    case "${enableval}" in
      yes) cs_have_openmp=yes ;;
      no)  cs_have_openmp=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
    esac
  ],
  [ cs_have_openmp=yes ]
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
    AC_MSG_CHECKING([for OpenMP simd])
    # Note that as the m4 language does not allow escaping [ and ],
    # we use quadrigraphs @<:@ and @>:@ instead here.
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <omp.h>]],
                     [[double a@<:@128@:>@;
                      #pragma omp for simd safelen(4)
                      for (int i = 0; i < 128; i++) a@<:@i@:>@ = i;]])],
                     [cs_have_openmp_simd=yes],
                     [cs_have_openmp_simd=no])
    AC_MSG_RESULT($cs_have_openmp_simd)
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([for OpenMP (Fortran)])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
                   [[
      integer n
      n = omp_get_num_threads()
]])],
                   [cs_have_openmp_f=yes],
                   [cs_have_openmp_f=no])
    AC_MSG_RESULT($cs_have_openmp_f)
    AC_LANG_POP([Fortran])
  fi
  if test "x$cs_have_openmp" = "xyes" ; then
    AC_DEFINE([HAVE_OPENMP], 1, [openmp support])
    if test "x$cs_have_openmp_simd" = "xyes" ; then
      AC_DEFINE([HAVE_OPENMP_SIMD], 1, [openmp simd support])
    fi
    if test "x$cs_have_openmp_f" = "xno" ; then
      FCFLAGS="$saved_FCFLAGS"
    fi
  else
    cs_have_openmp=no
    CFLAGS="$saved_CFLAGS"
    FCFLAGS="$saved_FCFLAGS"
    LDFLAGS="$saved_LDFLAGS"
  fi
fi
AC_SUBST(cs_have_openmp)

# Now that compiler options are determined (relative notably to OpenMP),
# update Fortran libraries necessary to link.

if test "x$cs_have_openmp_f" = "xyes" ; then
  AC_FC_LIBRARY_LDFLAGS
fi

#------------------------------------------------------------------------------
# Checks for Python support.
#------------------------------------------------------------------------------

AC_ARG_VAR([PYTHON], [the Python interpreter])
AM_PATH_PYTHON(["2.6"])

# On MinGW hosts, Automake/Python mixes Windows-style and GNU-style paths
# It seems better to use the standard way of installing Python modules,
# as it is on Unix systems (whichever is chosen, Python modules will be
# freezed by cx_freeze for Windows packaging).
# So, we choose to override some Python paths.
if test "$host_os" = mingw32 ; then

  AC_MSG_NOTICE([Overriding Python paths on MinGW hosts])

  am_cv_python_pythondir=\${prefix}/lib/python$PYTHON_VERSION/site-packages
  AC_MSG_CHECKING([for $PYTHON script directory])
  pythondir=$am_cv_python_pythondir
  pkgpythondir=\${pythondir}/$PACKAGE
  AC_MSG_RESULT([$pythondir])

  am_cv_python_pyexecdir=\${exec_prefix}/lib/python$PYTHON_VERSION/site-packages
  AC_MSG_CHECKING([for $PYTHON extension module directory])
  pyexecdir=$am_cv_python_pyexecdir
  pkgpyexecdir=\${pyexecdir}/$PACKAGE
  AC_MSG_RESULT([$pyexecdir])
  AC_SUBST(pythonversion, [${PYTHON_VERSION}])
fi

win_prefix=""
if test "$host_os" = mingw32 ; then
  win_prefix=`cygpath --path --windows "${prefix}"`
  AC_SUBST(win_prefix, [${win_prefix}])
fi

#------------------------------------------------------------------------------
# Checks for additional version information.
#------------------------------------------------------------------------------

AC_ARG_ENABLE(cs_use_revision, [AS_HELP_STRING([--disable-revision],
                                               [disable Subversion or Git revision info])],
  [
    case "${enableval}" in
      yes) cs_use_revision=yes ;;
      no)  cs_use_revision=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-revision]) ;;
    esac
  ],
  [ cs_use_revision=yes ]
)

cs_src_revision=""
if test "x$cs_use_revision" = "xyes" ; then
  if test -d ${srcdir}/.svn -o -d ${srcdir}/.git ; then
    cs_src_revision=`$PYTHON ${srcdir}/build-aux/cs_version.py --revision-only`
  elif test -f ${srcdir}/build-aux/version ; then
    cs_src_revision=`cat ${srcdir}/build-aux/version`
  fi
fi

AC_DEFINE_UNQUOTED([CS_REVISION], ["$cs_src_revision"], [Subversion / Git revision string])
AC_SUBST(cs_src_revision)

# Define the different version numbers for later use

cs_version=`$PYTHON ${srcdir}/build-aux/cs_version.py --short`
cs_version_full=`$PYTHON ${srcdir}/build-aux/cs_version.py --full`
cs_version_short=`$PYTHON ${srcdir}/build-aux/cs_version.py --short`

AC_SUBST(cs_version)
AC_SUBST(cs_version_full)
AC_SUBST(cs_version_short)

#------------------------------------------------------------------------------
# Optional shared library and plugin support
#------------------------------------------------------------------------------

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

# Use long global numbers ?
AC_ARG_ENABLE(dlopen-rtld-global,
  [AS_HELP_STRING([--enable-dlopen-rtld-global], [add RTLD_GLOBAL to dlopen flags])],
  [
    case "${enableval}" in
      yes)  AC_DEFINE([CS_DLOPEN_USE_RTLD_GLOBAL], 1,
                       [use RTLD_GLOBAL in dlopen flags.]) ;;
      no) ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-dlopen-rtld-global]) ;;
    esac
  ],
  [ ]
)

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

CS_AC_TEST_PLE(["2.0.0"])

AM_CONDITIONAL(HAVE_INTERNAL_PLE, test x$cs_have_internal_ple = xyes)
if test "x$cs_have_internal_ple" = xyes; then
  AC_CONFIG_SUBDIRS([libple])
fi

CS_AC_TEST_SALOME

CS_AC_TEST_LIBXML2
CS_AC_TEST_BLAS([$cs_have_openmp])
CS_AC_TEST_MPI
CS_AC_TEST_METIS
CS_AC_TEST_SCOTCH
CS_AC_TEST_HDF5
CS_AC_TEST_CGNS
CS_AC_TEST_MED
CS_AC_TEST_CCM
CS_AC_TEST_EOS
CS_AC_TEST_FREESTEAM
CS_AC_TEST_COOLPROP
CS_AC_TEST_ZLIB
CS_AC_TEST_MEDCOUPLING
CS_AC_TEST_CATALYST
CS_AC_TEST_MELISSA
CS_AC_TEST_PETSC
CS_AC_TEST_VOFI

CS_AC_TEST_DOCS

cs_have_neptune=cfd=no
if test -f ${srcdir}/neptune_cfd/configure; then
  cs_have_neptune_cfd=yes
fi

# Check if we need to force link with C++
cs_have_link_cxx=no
if test x$cs_have_neptune_cfd = xyes ; then
  cs_have_link_cxx=yes
elif test x$cs_have_med_link_cxx = xyes ; then
  cs_have_link_cxx=yes
elif test x$cs_have_medcoupling = xyes -o x$cs_have_paramedmem = xyes ; then
  if test x$enable_shared = xno -o x$cs_have_plugin_medcoupling = xno ; then
    cs_have_link_cxx=yes
  fi
elif test x$cs_have_eos = xyes ; then
  cs_have_link_cxx=yes
elif test x$cs_have_coolprop = xyes ; then
  cs_have_link_cxx=yes
elif test x$cs_have_catalyst = xyes ; then
  if test x$enable_shared = xno -o x$cs_have_plugin_catalyst = xno ; then
    cs_have_link_cxx=yes
  fi
fi
if test x$user_CS_LD != "x" ; then
  CS_LD=$user_CS_LD
elif test $cs_have_link_cxx = yes ; then
  CS_LD=$CXX
else
  CS_LD=$CC
fi
AC_SUBST(CS_LD)
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
# Enables Code_Saturne back-end
#------------------------------------------------------------------------------

AC_ARG_ENABLE(backend, [AS_HELP_STRING([--disable-backend],
                                        [disable back-end elements])],
  [
    case "${enableval}" in
      yes) cs_have_backend=yes ;;
      no)  cs_have_backend=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-backend]) ;;
    esac
  ],
  [ cs_have_backend=yes ]
)
AM_CONDITIONAL(HAVE_BACKEND, test x$cs_have_backend = xyes)
AC_SUBST(cs_have_backend)

#------------------------------------------------------------------------------
# Enables Code_Saturne graphical user interface
#------------------------------------------------------------------------------

AC_ARG_ENABLE(gui, [AS_HELP_STRING([--disable-gui],
                                   [disable the Graphical User Interface])],
  [
    case "${enableval}" in
      yes) cs_have_gui=check ;;
      no)  cs_have_gui=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-gui]) ;;
    esac
  ],
  [ cs_have_gui=check ]
)

if test "x$cs_have_frontend" = "xno" ; then
  cs_have_gui=no
fi

#------------------------------------------------------------------------------
# Checks for PyQt support.
#------------------------------------------------------------------------------

# Test for PyQt4 or PyQt5
cs_have_gui_qt4=no
cs_have_gui_qt5=no

if test "x$cs_have_gui" = "xcheck" ; then

  if test "x$QT_SELECT" = "x" -o "x$QT_SELECT" = "x5" ; then

    AC_MSG_CHECKING([for PyQt5 version >= 5.0])

    prog="[import sys, string
try: import PyQt5
except ImportError: sys.exit(1)
from PyQt5.QtCore import *
if list(map(int, QT_VERSION_STR.split('.'))) < [5,0,0]: sys.exit(1)
if list(map(int, PYQT_VERSION_STR.split('.'))) < [5,0,0]: sys.exit(1)
sys.exit(0)]"

    ${PYTHON} -c "${prog}"
    retval=$?

    if test $retval -ne 0 ; then
      AC_MSG_RESULT([no])
    else
      AC_MSG_RESULT([yes])
      cs_have_gui=yes
      cs_have_gui_qt5=yes
    fi

  fi

fi

if test "x$cs_have_gui" = "xcheck" ; then

  cs_have_gui=no

  if test "x$QT_SELECT" = "x" -o "x$QT_SELECT" = "x4" ; then

    AC_MSG_CHECKING([for PyQt4 version >= 4.5])

    prog="[import sys, string
try: import PyQt4
except ImportError: sys.exit(1)
from PyQt4.QtCore import *
if list(map(int, QT_VERSION_STR.split('.'))) < [4,5,0]: sys.exit(1)
if list(map(int, PYQT_VERSION_STR.split('.'))) < [4,5,0]: sys.exit(1)
sys.exit(0)]"

    ${PYTHON} -c "${prog}"
    retval=$?

    if test $retval -ne 0 ; then
      AC_MSG_RESULT([no])
      AC_MSG_FAILURE([cannot find PyQt5 support (>= 5.0) or PyQt4 support (>= 4.5), Graphical User Interface cannot be installed])
    else
      AC_MSG_RESULT([yes])
      cs_have_gui=yes
      cs_have_gui_qt4=yes
    fi
  fi

fi

# Test for PyQt4 developer tools
if test "x$cs_have_gui_qt4" = "xyes" ; then

  # Try to find PyQt4 tools through Python interpreter
  prog="[import os, sys
from PyQt4.QtCore import QCoreApplication
app = QCoreApplication([])
path = app.applicationDirPath()
sys.stdout.write(path)]"

  PYPATH=`${PYTHON} -c "${prog}"`

  # On Cygwin hosts, the PATH does not seem to be correctly found
  if test "$host_os" = cygwin ; then
     PYPATH=/usr/lib/python$PYTHON_VERSION/site-packages/PyQt4
  fi

  AC_ARG_VAR([PYUIC4], [PyQt4 user interfaces compiler])
  AC_ARG_VAR([PYRCC4], [PyQt4 resources compiler])
  AC_ARG_VAR([LRELEASE], [Qt message files generator])

  # On MinGW (MSYS) or Cygwin (BASH) build systems, and for MinGW hosts,
  # pyuic4 is a batch file and cannot be run from the shell.
  # To bypass this issue, we launch pyuic4 through cmd.exe
  # with a "standard" PATH.
  # pyrcc4 does not need this trick as it is a standard executable.
  if test "$host_os" = mingw32 ; then
    if test "$build_os" = mingw32 ; then
      CMDPYUIC4="cmd //C $PYPATH/Lib/site-packages/PyQt4/pyuic4"
    elif test "$build_os" = cygwin ; then
      CMDPYUIC4="cmd /C $PYPATH/Lib/site-packages/PyQt4/pyuic4"
    fi
  fi

  if test "x$PYUIC4" = "x" ; then
    AC_PATH_PROG([PYUIC4], [pyuic4], [$CMDPYUIC4], [${PYPATH}:$PATH])
  fi
  if test "x$PYRCC4" = "x" ; then
    AC_PATH_PROG([PYRCC4], [pyrcc4], [],           [${PYPATH}:$PATH])
  fi

  if test "x$LRELEASE" = "x" ; then
    AC_PATH_PROGS([LRELEASE], [lrelease-qt4 lrelease])
  fi

  # lrelease not currently used, so not an issue if missing

  if test "x$PYUIC4" = "x" -o "x$PYRCC4" = "x" ; then
    AC_MSG_FAILURE([cannot find PyQt4 dev tools, Graphical User Interface cannot be installed])
    cs_have_gui=no
  fi

fi

# Test for PyQt5 developer tools
if test "x$cs_have_gui_qt5" = "xyes" ; then

  # Try to find PyQt5 tools through Python interpreter
  prog="[import os, sys
from PyQt5.QtCore import QCoreApplication
app = QCoreApplication([])
path = app.applicationDirPath()
sys.stdout.write(path)]"

  PYPATH=`${PYTHON} -c "${prog}"`

  ## On Cygwin hosts, the PATH does not seem to be correctly found
  #if test "$host_os" = cygwin ; then
  #   PYPATH=/usr/lib/python$PYTHON_VERSION/site-packages/PyQt5
  #fi

  AC_ARG_VAR([PYUIC5], [PyQt5 user interfaces compiler])
  AC_ARG_VAR([PYRCC5], [PyQt5 resources compiler])
  AC_ARG_VAR([LRELEASE], [Qt message files generator])

  ## On MinGW (MSYS) or Cygwin (BASH) build systems, and for MinGW hosts,
  ## pyuic5 is a batch file and cannot be run from the shell.
  ## To bypass this issue, we launch pyuic5 through cmd.exe
  ## with a "standard" PATH.
  ## pyrcc5 does not need this trick as it is a standard executable.
  #if test "$host_os" = mingw32 ; then
  #  if test "$build_os" = mingw32 ; then
  #    CMDPYUIC4="cmd //C $PYPATH/Lib/site-packages/PyQt4/pyuic4"
  #  elif test "$build_os" = cygwin ; then
  #    CMDPYUIC4="cmd /C $PYPATH/Lib/site-packages/PyQt4/pyuic4"
  #  fi
  #fi

  if test "x$PYUIC5" = "x" ; then
    AC_PATH_PROG([PYUIC5], [pyuic5], [$CMDPYUIC5], [${PYPATH}:$PATH])
  fi
  if test "x$PYRCC5" = "x" ; then
    AC_PATH_PROG([PYRCC5], [pyrcc5], [],           [${PYPATH}:$PATH])
  fi

  if test "x$LRELEASE" = "x" ; then
    AC_PATH_PROGS([LRELEASE], [lrelease-qt5 lrelease])
  fi

  # lrelease not currently used, so not an issue if missing

  if test "x$PYUIC5" = "x" -o "x$PYRCC5" = "x" ; then
    AC_MSG_FAILURE([cannot find PyQt5 dev tools, Graphical User Interface cannot be installed])
    cs_have_gui=no
  fi

fi

#------------------------------------------------------------------------------
# Checks for Qt tools.
#------------------------------------------------------------------------------

AM_CONDITIONAL(HAVE_GUI, [test "x$cs_have_gui" = "xyes"])
AC_SUBST(cs_have_gui)
AM_CONDITIONAL(HAVE_QT4, [test "x$cs_have_gui_qt4" = "xyes"])
AC_SUBST(cs_have_gui_qt4)
AM_CONDITIONAL(HAVE_QT5, [test "x$cs_have_gui_qt5" = "xyes"])
AC_SUBST(cs_have_gui_qt5)

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([sys/types.h sys/utsname.h sys/stat.h dirent.h stddef.h])
AC_CHECK_HEADERS([unistd.h fcntl.h sys/types.h sys/signal.h])
AC_CHECK_HEADERS([sys/procfs.h sys/sysinfo.h sys/resource.h])
AC_CHECK_HEADERS([float.h string.h sys/time.h])

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([mkdir chdir stat dup2 access])
AC_CHECK_FUNCS([snprintf])
AC_CHECK_FUNCS([getcwd sleep])
AC_CHECK_FUNCS([getpwuid geteuid])
AC_CHECK_FUNCS([uname])
AC_CHECK_FUNCS([clock_gettime clock_getcpuclockid])
AC_CHECK_FUNCS([getrusage gettimeofday sbrk sysinfo])
AC_CHECK_FUNCS([posix_memalign])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([strtok_r])

saved_LIBS="$LIBS"
LIBS="${LIBS} -lm"
AC_CHECK_FUNCS([pow modf tgamma erf])
LIBS="$saved_LIBS"

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

# Optional socket support
#------------------------

have_socket=no

AC_ARG_ENABLE(sockets,
  [  --disable-sockets       disable communications through IP sockets],
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

# Plugin modules support
#-----------------------

cs_have_plugins=no
cs_py_have_plugins=False
if test "x$enable_shared" = "xyes" ; then
  if test "x$cs_have_dlloader" = "xdlopen" ; then
    # Test for different plugins here
    if test x$cs_have_plugin_medcoupling = xyes ; then
      cs_have_plugins=yes
    elif test x$cs_have_plugin_catalyst = xyes ; then
      cs_have_plugins=yes
    elif test x$cs_have_ccm = xyes ; then
      cs_have_plugins=yes
    fi
  fi
fi

if test "x$cs_have_plugins" = "xyes" ; then
  cs_py_have_plugins=True
  AC_DEFINE([HAVE_PLUGINS], 1, [plugins support])
fi

AC_SUBST(cs_py_have_plugins)
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

NEPTUNE_CFD=
if test x$cs_have_neptune_cfd = xyes; then
  AC_CONFIG_SUBDIRS([neptune_cfd])
  NEPTUNE_CFD=neptune_cfd
  NEPTUNE_CFD_LIBS="-lneptune"
fi

AC_SUBST(HAVE_NEPTUNE_CFD)
AC_SUBST(NEPTUNE_CFD)

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------

# Main configuration

AC_CONFIG_FILES([Makefile po/Makefile.in
                 preprocessor/Makefile
                 preprocessor/appli/Makefile preprocessor/base/Makefile
                 preprocessor/pre-post/Makefile preprocessor/util/Makefile
                 salome/Makefile
                 src/Makefile      src/apps/Makefile src/gui/Makefile
                 src/bft/Makefile  src/fvm/Makefile  src/mei/Makefile
                 src/base/Makefile src/atmo/Makefile src/cdo/Makefile
                 src/cfbl/Makefile src/cogz/Makefile src/ctwr/Makefile
                 src/darc/Makefile src/elec/Makefile src/comb/Makefile
                 src/pprt/Makefile src/lagr/Makefile src/rayt/Makefile
                 src/turb/Makefile src/alge/Makefile src/mesh/Makefile
                 src/user/Makefile src/user_examples/Makefile
                 gui/Makefile gui/Base/Makefile gui/Pages/Makefile
                 gui/studymanager_gui/Makefile gui/trackcvg/Makefile
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
                 salome/cfd_study/src/CFDSTUDYOTURNS/Makefile
                 docs/Makefile docs/doxygen/Makefile docs/style/csvers.tex
                 docs/studymanager/Makefile docs/studymanager/graphics/Makefile
                 docs/developer/Makefile docs/developer/graphics/Makefile
                 docs/install/Makefile
                 docs/refcard/Makefile docs/refcard/graphics/Makefile
                 docs/theory/Makefile docs/theory/graphics/Makefile
                 docs/user/Makefile docs/user/graphics/Makefile
                 tests/Makefile
                 extras/salome/fsi_appli_config.xml
                 extras/windows/innosetup/setup_mingw.iss
                 extras/windows/innosetup/setup_cygwin32.iss
                 extras/windows/innosetup/setup_cygwin64.iss
                 bin/cs_config.py])

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
if test x$cs_have_neptune_cfd = xyes ; then
  echo " NEPTUNE_CFD module support"
fi
echo " use debugging code: "$debug""
echo " use malloc hooks: "$cs_have_malloc_hooks""
if test x$cs_have_frontend = xno ; then
  echo " use frontend: no"
else
  echo " use graphical user interface: "$cs_have_gui""
fi
if test x$cs_have_backend = xno ; then
  echo " use backend: no"
fi
echo " use long global numbers: "$cs_have_long_gnum""
echo " Zlib (gzipped file) support: "$cs_have_zlib""
echo " MPI (Message Passing Interface) support: "$cs_have_mpi""
if (test x$cs_have_mpi = xyes) ; then
  echo "   MPI I/O support: "$cs_have_mpi_io""
  echo "   MPI2 one-sided communication support: "$cs_have_mpi_one_sided""
  echo "   MPI3 neighborhood collectives support: "$cs_have_mpi_one_sided""
  echo "   MPI3 nonblocking barrier support: "$cs_have_mpi_ibarrier""
fi
echo " OpenMP support: "$cs_have_openmp""
if test x$cs_have_openmp = xyes ; then
  echo " OpenMP Fortran support: "$cs_have_openmp_f""
fi
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
echo " CCM support: "$cs_have_ccm""
echo " HDF (Hierarchical Data Format) support: "$cs_have_hdf5""
echo " CGNS (CFD General Notation System) support: "$cs_have_cgns""
echo " MED (Model for Exchange of Data) support: "$cs_have_med""
if (test x$cs_have_med = xyes) ; then
  echo "   MED MPI I/O support: "$cs_have_med_mpi""
fi
echo " MEDCoupling support: "$cs_have_medcoupling""
if test x$cs_have_medcoupling = xyes ; then
  echo "   MEDLoader support: "$cs_have_medcoupling_loader""
  echo "   ParaMEDMEM support: "$cs_have_paramedmem""
  echo "   MEDCoupling as plugin: "$cs_have_plugin_medcoupling""
fi
echo " Catalyst (ParaView co-processing) support: "$cs_have_catalyst""
if (test x$cs_have_catalyst = xyes) ; then
  echo "   Catalyst as plugin: "$cs_have_plugin_catalyst""
fi
echo " Melissa (in-situ statistics) support: "$cs_have_melissa""
if (test x$cs_have_melissa = xyes) ; then
  echo "   Melissa as plugin: "$cs_have_plugin_melissa""
fi
echo " Vofi (VoF initialization) support: "$cs_have_vofi""
echo " EOS support: "$cs_have_eos""
echo " freesteam support: "$cs_have_freesteam""
echo " CoolProp support: "$cs_have_coolprop""
echo " PETSc support: "$cs_have_petsc""
echo " SALOME GUI support: "$cs_have_salome_gui""
echo " SALOME Kernel support: "$cs_have_salome_kernel""
echo " Dynamic loader support: "$cs_have_dlloader""
if (test x$cs_env_modules != x) ; then
  echo " Auto load environment modules: "$cs_env_modules""
fi
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the Code_Saturne kernel"
echo


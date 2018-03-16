# Process this file with autoconf to produce a configure script.

#------------------------------------------------------------------------------
# This file is part of the "Parallel Location and Exchange" library,
# intended to provide mesh or particle-based code coupling services.
#
# Copyright (C) 2005-2018 EDF S.A.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([ple_licence_c_comment],
[/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2018  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
*/])

m4_define([ple_major_version], [2])
m4_define([ple_minor_version], [0])
m4_define([ple_release_version], [1])
m4_define([ple_version_extra], [-patch])
m4_define([ple_version_string],
          [ple_major_version.ple_minor_version.ple_release_version@&t@ple_version_extra])

AC_INIT([ple],[ple_version_string],[saturne-support@edf.fr])
AC_CONFIG_SRCDIR([src/ple_defs.h])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([ple_config_priv.h:ple_config_priv.h.in])

AM_GNU_GETTEXT_VERSION([0.18.1])

# Initialize automake with the following options:
# - foreign, so as to handle the possible absence of ChangeLog
# - tar-pax, so as to handle long lines (> 99 characters) in tar archives
# Warnings can be activated at bootstrap with 'autoreconf -vi --warnings=all'
AM_INIT_AUTOMAKE([foreign tar-pax])

AH_TOP([#ifndef __PLE_CONFIG_PRIV_H__
#define __PLE_CONFIG_PRIV_H__

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.
*/])

AH_BOTTOM([#endif /* __PLE_CONFIG_PRIV_H__ */])

#Define public configuration file

PLE_AC_CONFIG_PUBL_INIT([ple_config.h])
PLE_AC_CONFIG_PUBL_VERBATIM([ple_licence_c_comment])

# Enable maintainer mode by default for a developer checkout
AS_IF([test -d ${srcdir}/.svn -o -d ${srcdir}/.git],
      [AM_MAINTAINER_MODE([enable])],
      [AM_MAINTAINER_MODE([disable])])

# Enable silent rules with "./configure --enable-silent-rules" or "make V=0"
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

# Version information

PLE_AC_CONFIG_PUBL_DEFINE_STRING(PLE_VERSION, ple_version_string,
                                 [Complete version number])

PLE_AC_CONFIG_PUBL_DEFINE(PLE_MAJOR_VERSION, ple_major_version,
                          [Major version number])
PLE_AC_CONFIG_PUBL_DEFINE(PLE_MINOR_VERSION, ple_minor_version,
                          [Minor version number])
PLE_AC_CONFIG_PUBL_DEFINE(PLE_RELEASE_VERSION, ple_release_version,
                          [Release version number])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_LDFLAGS=$LDFLAGS
user_LIBS=$LIBS

AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S

AM_PROG_CC_C_O

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

  if test -f "$srcdir/config/ple_auto_flags.sh" ; then
    AC_MSG_NOTICE([sourcing config/ple_auto_flags.sh (test for known compilers)])
    . "$srcdir/config/ple_auto_flags.sh"
  else
    AC_MSG_WARN([config/ple_auto_flags.sh default configuration file not found])
  fi

  # Default flags
  CPPFLAGS="$cppflags_default $user_CPPFLAGS"
  CFLAGS="$cflags_default $user_CFLAGS"
  LDFLAGS="$ldflags_default $user_LDFLAGS"
  LIBS="$libs_default $user_LIBS"
  LDRPATH="$ldflags_rpath"

  if test "x$debug" = xyes; then
    # Debug flags
    CFLAGS_DBG="$cflags_default_dbg"
    # Add debug flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_dbg"
  else
    # Normal optimization flags
    CFLAGS_OPT="$cflags_default_opt"
    # Add optimization flags for linker
    LDFLAGS="$LDFLAGS $ldflags_default_opt"
  fi

  AC_SUBST(CFLAGS_DBG)
  AC_SUBST(CFLAGS_OPT)

  AC_SUBST(LDRPATH)

fi

LIBS="${LIBS} -lm"

# Preprocessor flags for debugging purposes

if test "x$debug" = "xyes"; then
  CPPFLAGS="${CPPFLAGS} -DDEBUG"
else
  CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi

#------------------------------------------------------------------------------

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

# Add some preprocessor flags on some systems to work around autotools issues
#----------------------------------------------------------------------------

if test "x$enable_shared" = "xyes" ; then
  if test "x$auto_flags" = "xyes" ; then
    case "$ple_ac_cc_version" in
      Cray* )
        CFLAGS="${CFLAGS} -fPIC"
        ;;
      *)
        ;;
    esac
  fi
fi

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([sys/types.h stddef.h])

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

PLE_AC_CHECK_SIZEOF([void *], [ple])
PLE_AC_CHECK_SIZEOF([short], [ple])
PLE_AC_CHECK_SIZEOF([int], [ple])
PLE_AC_CHECK_SIZEOF([long], [ple])
PLE_AC_CHECK_SIZEOF([long long], [ple])
PLE_AC_CHECK_SIZEOF([float], [ple])
PLE_AC_CHECK_SIZEOF([double], [ple])

#------------------------------------------------------------------------------
# Determine OpenMP support (not used in source code yet)
#------------------------------------------------------------------------------

# For consistency with the auto_flags mechanisme and with Code_Saturne,
# we do not use AC_OPENMP, our own method here, based on flags already set
# by default or by the user.

# From this point, no runtime tests should be run, so no issues should
# arise due to missing RPATH or LD_LIBRARY_PATH values in case of non-standard
# compiler install paths, such as may happen with gcc's libgomp.

ple_have_openmp=no

AC_ARG_ENABLE(openmp,
  [AS_HELP_STRING([--disable-openmp], [disable OpenMP support])],
  [
    case "${enableval}" in
      yes) ple_have_openmp=yes ;;
      no)  ple_have_openmp=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
    esac
  ],
  [ ple_have_openmp=no ]
)

if test "x$ple_have_openmp" = "xyes" ; then

  saved_CFLAGS="$CFLAGS"
  saved_FCFLAGS="$FCFLAGS"
  saved_LDFLAGS="$LDFLAGS"

  CFLAGS="${CFLAGS} ${cflags_default_omp}"
  FCFLAGS="${FCFLAGS} ${fcflags_default_omp}"
  LDFLAGS="${LDFLAGS} ${cflags_default_omp}"

  AC_MSG_CHECKING([for OpenMP (C)])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <omp.h>]],
                 [[ omp_get_num_threads();]])],
                 [ple_have_openmp=yes],
                 [ple_have_openmp=no])
  AC_MSG_RESULT($ple_have_openmp)
  if test "x$ple_have_openmp" = "xyes" ; then
    AC_DEFINE([HAVE_OPENMP], 1, [openmp support])
  else
    ple_have_openmp=no
    CFLAGS="$saved_CFLAGS"
    FCFLAGS="$saved_FCFLAGS"
    LDFLAGS="$saved_LDFLAGS"
  fi
fi
AC_SUBST(ple_have_openmp)

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([snprintf gettimeofday])

#------------------------------------------------------------------------------
# Checks for system services.
#------------------------------------------------------------------------------

#---------------------
# optional MPI support (use CC=mpicc with configure if necessary)
#---------------------

MPI_CPPFLAGS=""
MPI_LDFLAGS=""
MPI_LIBS=""

PLE_AC_TEST_MPI

# Whether we used mpicc or $CC, we now know if we have MPI support
PLE_AC_CONFIG_PUBL_SET(PLE_HAVE_MPI, $ple_have_mpi,
                       [PLE library configured with MPI support ?])

#------------------------------------------------------------------------------
# Internationalization
#------------------------------------------------------------------------------

AM_GNU_GETTEXT([external])

#------------------------------------------------------------------------------
# Documentation
#------------------------------------------------------------------------------

ple_have_doxygen=yes

AC_ARG_VAR([DOXYGEN], [source code documentation generator])

dnl where is doxygen ?
AC_PATH_PROG(DOXYGEN, [doxygen])
if test "x$DOXYGEN" = "x"; then
  AC_MSG_WARN(doxygen not found)
  ple_have_doxygen=no
fi

AC_ARG_VAR([DOT], [graphs generator])

dnl where is dot ?
AC_PATH_PROG(DOT, [dot])
if test "x$DOT" = "x"; then
  AC_MSG_WARN(dot not found; Doxygen documentation will not be built)
  ple_have_doxygen=no
fi

AM_CONDITIONAL(HAVE_DOXYGEN, [test $ple_have_doxygen = yes])
AC_SUBST(ple_have_doxygen)

#------------------------------------------------------------------------------
# Compiler options and system characteristics
#------------------------------------------------------------------------------

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------

# Determine CPPFLAGS, and LDFLAGS required for external code using PLE
#---------------------------------------------------------------------

for flag in [$CPPFLAGS] ; do
  case [$flag] in
    -I* | -isystem*)
    user_cppflags="${user_cppflags} $flag"
    ;;
  esac
done

for flag in [$LDFLAGS] ; do
  case [$flag] in
    -L* | -Wl,-L*)
    user_ldflags="${user_ldflags} $flag"
    ;;
  esac
done

# Output configuration info script
#---------------------------------

PLE_AC_CONFIG_INFO_INIT([ple-config])

PLE_AC_CONFIG_INFO_EXTRA([
# use debugging code:                         "$debug"
# MPI (Message Passing Interface) support:    "$ple_have_mpi"
])

PLE_AC_CONFIG_INFO_CC([${CC}], [${ple_ac_cc_version}], [${ple_ac_cc_version_full}])

PLE_AC_CONFIG_INFO_VERSION([ple_version_string])

PLE_AC_CONFIG_INFO_DIRS([$prefix], [$exec_prefix], [$includedir], [$libdir])

PLE_AC_CONFIG_INFO_FLAGS([$user_cppflags], [], [$user_ldflags], [$LIBS])

PLE_AC_CONFIG_INFO_FLAGS([$CPPFLAGS], [$CFLAGS], [$LDFLAGS], [$LIBS],
                         [build])

PLE_AC_CONFIG_INFO_FLAGS([$MPI_CPPFLAGS], [], [$MPI_LDFLAGS], [$MPI_LIBS],
                         [mpi])

PLE_AC_CONFIG_INFO_FINALIZE

#------------------------------------------------------------------------------

# If we use GNU gettext, we need to add po/Makefile.

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile
                 doc/Makefile doc/doxygen/Makefile
                 po/Makefile.in])

#------------------------------------------------------------------------------

PLE_AC_CONFIG_PUBL_FINALIZE

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
echo " MPI (Message Passing Interface) support: "$ple_have_mpi""
echo " OpenMP support: "$ple_have_openmp""
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the PLE package"
echo


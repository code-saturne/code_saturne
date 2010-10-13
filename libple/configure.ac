# Process this file with autoconf to produce a configure script.

#------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2010  EDF
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Initialization
#------------------------------------------------------------------------------

m4_define([ple_licence_c_comment],
[/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2010  EDF

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
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/])

m4_define([ple_major_version], [1])
m4_define([ple_minor_version], [0])
m4_define([ple_release_version], [0])
m4_define([ple_version_extra], [])
m4_define([ple_version_string],
          [ple_major_version.ple_minor_version.ple_release_version@&t@ple_version_extra])

AC_INIT([ple],[ple_version_string],[yvan.fournier@edf.fr])
AC_CONFIG_SRCDIR([src/ple_defs.h])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_SUBST(ac_aux_dir)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([ple_config_priv.h:ple_config_priv.h.in])

AM_INIT_AUTOMAKE

AH_TOP([#ifndef __PLE_CONFIG_PRIV_H__
#define __PLE_CONFIG_PRIV_H__

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.
*/])

AH_BOTTOM([#endif /* __PLE_CONFIG_PRIV_H__ */])

#Define public configuration file

PLE_AC_CONFIG_PUBL_INIT([ple_config.h])
PLE_AC_CONFIG_PUBL_VERBATIM([ple_licence_c_comment])

AM_MAINTAINER_MODE

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

AC_PROG_CC(cc gcc)
AC_PROG_INSTALL
AC_PROG_LN_S

AM_PROG_CC_C_O

# Debug or production compilation mode (debug by default) ?

AC_ARG_ENABLE(debug,
  [  --enable-debug          enable debugging (reduces optimization)],
  [
    case "${enableval}" in
      yes) debug=true ;;
      no)  debug=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
    esac
  ],
  [ debug=false ]
)
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

# Optionally deactivate automatic determination of flags on known systems

AC_ARG_ENABLE(auto-flags,
  [  --disable-auto-flags    do not define *FLAGS on known systems],
  [
    case "${enableval}" in
      yes) auto_flags=true ;;
      no)  auto_flags=false ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-auto-flags]) ;;
    esac
  ],
  [ auto_flags=true ]
)

# Default compiler options (may be modified
# by defining CFLAGS in the environment)
#------------------------------------------

if test "x$auto_flags" = "xtrue" ; then

  # Source associated recommended compiler options

  if test -f "$srcdir/config/ple_auto_flags.sh" ; then
    AC_MSG_NOTICE([sourcing config/ple_auto_flags.sh (test for known compilers)])
    . "$srcdir/config/ple_auto_flags.sh"
  else
    AC_MSG_WARN([config/ple_auto_flags.sh default configuration file not found])
  fi

  if test "x$user_CPPFLAGS" = "x"; then
    CPPFLAGS="$cppflags_default"
  else
    AC_MSG_NOTICE([default CPPFLAGS overridden by user definition])
  fi

  if test "x$user_CFLAGS" = "x"; then

    CFLAGS="$cflags_default"

    if test "x$debug" = xtrue; then
      CFLAGS="${CFLAGS} $cflags_default_dbg"
    else
      CFLAGS="${CFLAGS} $cflags_default_opt"
    fi

  else

    AC_MSG_NOTICE([default CFLAGS overridden by user definition])

  fi

  if test "x$user_LDFLAGS" = "x"; then

    LDFLAGS="$ldflags_default"

    if test "x$debug" = xtrue; then
      LDFLAGS="${LDFLAGS} $ldflags_default_dbg"
    else
      LDFLAGS="${LDFLAGS} $ldflags_default_opt"
    fi

  else

    AC_MSG_NOTICE([default LDFLAGS overridden by user definition])

  fi

fi

LIBS="${LIBS} -lm"

# Optional additional user flags
#-------------------------------

AC_ARG_VAR(CPP_ADD,  [Additional user C preprocessor flags])
AC_ARG_VAR(CC_ADD,   [Additional user C compiler flags])
AC_ARG_VAR(LD_ADD,   [Additional user linker flags, e.g. -L<libdir>])
AC_ARG_VAR(LIBS_ADD, [Additional user librairies, e.g. -l<lib>])

if test "x$CPP_ADD" != x ; then
  CPPFLAGS="${CPPFLAGS} $CPP_ADD"
fi

if test "x$CC_ADD" != x ; then
  CFLAGS="${CFLAGS} $CC_ADD"
fi

if test "x$LD_ADD" != x ; then
  LDFLAGS="${LDFLAGS} $LD_ADD"
fi

if test "x$LIBS_ADD" != x ; then
  LIBS="${LIBS} $LIBS_ADD"
fi

#------------------------------------------------------------------------------

# We may only turn on processing for libtool now that the basic compiler
# and linker flags are set (to avoid issues with linkers with different
# modes such as 32 and 64 bit which may be modified by compiler or
# linker flags).

# AC_DISABLE_SHARED may not be used inside a test, as it seems to
# take effect whether the result is true or not (maybe due to a bug),
# so we directly use autoconf's enable_shared variable.

if test "x$ple_disable_shared" = "xyes" ; then
  enable_shared=no
fi

LT_INIT

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

if test "${prefix}" = "NONE"; then
  localedir="/usr/local/share/locale"
else
  localedir="${prefix}/share/locale"
fi
AC_DEFINE_UNQUOTED(LOCALEDIR, "${localedir}", LOCALEDIR)

# Set of available languages now in po/LINGUAS.

# AM_GNU_GETTEXT

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

# If we use GNU gettext, we need intl/Makefile, po/Makefile, and m4/Makefile.

# AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile doc/Makefile
#                  intl/Makefile po/Makefile.in m4/Makefile])

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])

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
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the PLE package"
echo


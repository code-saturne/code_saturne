# Process this file with autoconf to produce a configure script.

#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009 EDF S.A., France
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

   Copyright (C) 2009 EDF S.A., France

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
m4_define([cs_version_extra], [-beta1])
m4_define([cs_version_string],
          [cs_major_version.cs_minor_version.cs_release_version@&t@cs_version_extra])

AC_INIT([ncs],[cs_version_string],[saturne-support@edf.fr])
AC_CONFIG_SRCDIR([src/base/cs_main.c])

# Use the config directory for libtool stuff ...
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_SUBST(ac_aux_dir)

AC_CANONICAL_TARGET
AC_CONFIG_HEADERS([cs_config.h])

AM_GNU_GETTEXT_VERSION([0.17])

AM_INIT_AUTOMAKE([-Wall])

AM_MAINTAINER_MODE

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------

# Check for C and Fortran compilers
#----------------------------------

user_CPPFLAGS=$CPPFLAGS
user_CFLAGS=$CFLAGS
user_FCFLAGS=$FCFLAGS
user_LDFLAGS=$LDFLAGS

AC_PROG_CC
AC_PROG_FC
AC_PROG_INSTALL
AC_PROG_LN_S

AC_FC_LIBRARY_LDFLAGS

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

# Optionally activate French translation of Fortran logs (compile-time)

AC_ARG_ENABLE(french,
  [  --enable-french         enable French translation of Fortran logs],
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
  LDRPATH="$ldflags_rpath"

  if test "x$debug" = xtrue; then
    CFLAGS_DBG="$cflags_default_dbg"
    FCFLAGS_DBG="$fcflags_default_dbg"
    LDFLAGS_DBG="$ldflags_default_dbg"
  else
    # Normal optimization flags
    CFLAGS_OPT="$cflags_default_opt"
    FCFLAGS_OPT="$fcflags_default_opt"
    LDFLAGS_OPT="$ldflags_default_opt"
    # Hot optimization flags
    CFLAGS_HOT="$cflags_default_hot"
    FCFLAGS_HOT="$fcflags_default_hot"
  fi

  AC_SUBST(CFLAGS_DBG)
  AC_SUBST(CFLAGS_OPT)
  AC_SUBST(CFLAGS_HOT)

  AC_SUBST(FCFLAGS_DBG)
  AC_SUBST(FCFLAGS_OPT)
  AC_SUBST(FCFLAGS_HOT)

  AC_SUBST(LDFLAGS_DBG)
  AC_SUBST(LDFLAGS_OPT)
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
  FCFLAGS="${FCFLAGS} -D_CS_LANG_FR"
else
  FCFLAGS="${FCFLAGS} -D_CS_LANG_FR"
fi

# We may only turn on processing for libtool now that the basic compiler
# and linker flags are set (to avoid issues with linkers with different
# modes such as 32 and 64 bit which may be modified by compiler or
# linker flags).

LT_INIT

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------

CS_AC_TEST_BFT(["1.0.0"])
CS_AC_TEST_FVM(["0.13.0"])
CS_AC_TEST_LIBXML2
CS_AC_TEST_BLAS
CS_AC_TEST_MPI
CS_AC_TEST_PREPRO
CS_AC_TEST_GUI
CS_AC_TEST_SYRCS
CS_AC_TEST_PYTHON

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS([sys/types.h sys/utsname.h sys/stat.h stddef.h])

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

AC_CHECK_SIZEOF([void *], [cs])
AC_CHECK_SIZEOF([short], [cs])
AC_CHECK_SIZEOF([int], [cs])
AC_CHECK_SIZEOF([long], [cs])
AC_CHECK_SIZEOF([float], [cs])
AC_CHECK_SIZEOF([double], [cs])

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------

AC_CHECK_FUNCS([snprintf])
AC_CHECK_FUNCS([getcwd sleep])
AC_CHECK_FUNCS([getpwuid geteuid])
AC_CHECK_FUNCS([uname])

#------------------------------------------------------------------------------
# Checks for system services.
#------------------------------------------------------------------------------

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

# Optional shared library support
#--------------------------------

have_dlloader=no

AC_ARG_ENABLE(dlloader,
  [  --disable-dlloader      disable dynamic shared library loading],
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
    LIBS="saved_LIBS"
    have_dlloader=no
  fi

fi

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

AM_GNU_GETTEXT([external])

#------------------------------------------------------------------------------
# Compiler options and system characteristics
#------------------------------------------------------------------------------

if test "$FCLIBS" != "x"; then
  LIBS="${LIBS} ${FCLIBS}"
fi

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
# MPI (Message Passing Interface) support:         "$cs_have_mpi"
# BLAS (Basic Linear Algebra Subprograms) support: "$have_blas"
# Libxml2 (XML Reader) support:                    "$have_libxml2"
# IP socket support (for SYRTHES 3 or CFD_Proxy): "$have_socket"
# Dynamic loader support (for YACS):               "$have_dlloader"
])

CS_AC_CONFIG_INFO_CC([${cs_ac_cc_version}], [${cs_ac_cc_version_full}])
CS_AC_CONFIG_INFO_FC([${cs_ac_fc_version}], [${cs_ac_fc_version_full}])

CS_AC_CONFIG_INFO_FINALIZE

AC_SUBST_FILE([cs_config_comments])
cs_config_comments=cs_config_comments_include

# Main configuration

AC_CONFIG_FILES([Makefile po/Makefile.in
                 src/base/Makefile src/atmo/Makefile src/cfbl/Makefile
                 src/cogz/Makefile src/cplv/Makefile src/ctwr/Makefile
                 src/elec/Makefile src/fuel/Makefile src/mati/Makefile
                 src/pprt/Makefile src/lagr/Makefile src/rayt/Makefile
                 bin/check_mesh bin/cs_config.py bin/cs_info bin/cs_create
                 bin/runcase bin/runcase_syrthes bin/SaturneGUI])

#------------------------------------------------------------------------------

AC_OUTPUT

rm -f cs_config_comments_include

echo
echo "Configuration options:"
echo " use debugging code: "$debug""
echo " MPI (Message Passing Interface) support: "$cs_have_mpi""
echo " BLAS (Basic Linear Algebra Subprograms) support: "$have_blas""
echo " Libxml2 (XML Reader) support: "$have_libxml2""
echo " IP socket support (for SYRTHES 3 or CFD_Proxy): "$have_socket""
echo " Dynamic loader support (for YACS): "$have_dlloader""
echo
echo "The package has been configured. Type:"
echo " make"
echo " make install"
echo
echo "To generate and install the Code_Saturne kernel"
echo


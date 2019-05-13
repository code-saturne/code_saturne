dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2019 EDF S.A.
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

# CS_AC_TEST_ZLIB
#----------------
# Checks for Zlib support
# modifies or sets cs_have_zlib, ZLIB_CPPFLAGS, ZLIB_LDFLAGS, and ZLIB_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_ZLIB],[

cs_have_zlib=no

AC_ARG_WITH(zlib,
            [AS_HELP_STRING([--with-zlib=PATH],
                            [specify prefix directory for Zlib])],
            [if test "x$withval" = "x"; then
               with_zlib=yes
             fi],
            [with_zlib=check])

AC_ARG_WITH(zlib-include,
            [AS_HELP_STRING([--with-zlib-include=PATH],
                            [specify directory for Zlib includes])],
            [if test "x$with_zlib" = "xcheck"; then
               with_zlib=yes
             fi
             ZLIB_CPPFLAGS="-I$with_zlib_include"],
            [if test "x$with_zlib" != "xno" -a "x$with_zlib" != "xyes" \
	          -a "x$with_zlib" != "xcheck"; then
               ZLIB_CPPFLAGS="-I$with_zlib/include"
             fi])

AC_ARG_WITH(zlib-lib,
            [AS_HELP_STRING([--with-zlib-lib=PATH],
                            [specify directory for Zlib library])],
            [if test "x$with_zlib" = "xcheck"; then
               with_zlib=yes
             fi
             ZLIB_LDFLAGS="-L$with_zlib_lib"],
            [if test "x$with_zlib" != "xno" -a "x$with_zlib" != "xyes" \
	          -a "x$with_zlib" != "xcheck"; then
               ZLIB_LDFLAGS="-L$with_zlib/lib"
             fi])

if test "x$with_zlib" != "xno" ; then

  saved_CPPFLAGS="$CPPFLAGS"
  saved_LDFLAGS="$LDFLAGS"
  saved_LIBS="$LIBS"

  CPPFLAGS="${CPPFLAGS} ${ZLIB_CPPFLAGS}"
  LDFLAGS="${LDFLAGS} ${ZLIB_LDFLAGS}"
  LIBS="-lz ${LIBS}"

  AC_MSG_CHECKING([for Zlib])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <zlib.h>]],
                 [[ gzFile *f; f = gzopen("filename", "r"); ]])],
                 [ AC_DEFINE([HAVE_ZLIB], 1, [gzipped file support])
                   cs_have_zlib=yes
                 ],
                 [if test "x$with_zlib" != "xcheck" ; then
                    AC_MSG_FAILURE([gzipped file support is requested, but test for Zlib failed!])
                  else
                    AC_MSG_WARN([no gzipped file support])
                  fi
                 ])
  AC_MSG_RESULT($cs_have_zlib)

  if test "x$cs_have_zlib" = "xno"; then
    LIBS="$saved_LIBS"
  fi

  # Additional test if zlib found to check for type sizes
  #------------------------------------------------------

  if test "x$cs_have_zlib" = "xyes"; then

    /bin/rm -f conftestval

    AC_MSG_CHECKING([size of z_off_t])
    AC_RUN_IFELSE([AC_LANG_SOURCE([
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
int main(int argc, char** argv)
{
  unsigned long i;
  int z_size_uint, z_size_ulong, z_size_voidpf, z_size_off_t;
  FILE *f = fopen("conftestval", "w");

  i = zlibCompileFlags(); /* See zlib.h for documentation */
  z_size_uint   = ((i >> 0) & 3) * 4;
  z_size_ulong  = ((i >> 2) & 3) * 4;
  z_size_voidpf = ((i >> 4) & 3) * 4;
  z_size_off_t  = ((i >> 6) & 3) * 4;

  if (   z_size_uint != sizeof(unsigned int)
      || z_size_ulong != sizeof(unsigned long)
      || z_size_voidpf != sizeof(void *)) {
    printf("\ncompile-time sizes for unsigned int, unsigned long, or void *\n"
           "for zlib are %d, %d and %d,\n"
           "while current values are %d, %d and %d.\n",
          z_size_uint, z_size_ulong, z_size_voidpf,
          sizeof(unsigned), sizeof(unsigned long), sizeof(void *));
    fprintf(f, "0\n");
  }
  else {
    fprintf(f, "%d\n", z_size_off_t);
  }

  fclose(f);
  exit(0);
}])],
[
if test -f conftestval ; then
  cs_ac_sizeof=`cat conftestval`
  if test "$cs_ac_sizeof" = "0" ; then
    AC_MSG_ERROR([zlib uint, ulong, or void sizes do not match current values.])
  else
    AC_MSG_RESULT([$cs_ac_sizeof])
    AC_DEFINE_UNQUOTED([SIZEOF_Z_OFF_T], $cs_ac_sizeof,
                       [The size of z_off_t, as returned by zlibCompileFlags.])
  fi
fi
],
[
AC_MSG_WARN([error running zlibCompileFlags configure test])
],
[
AC_MSG_WARN([unable to test for zlibCompileFlags when cross-compiling])
])

    unset cs_ac_sizeof

    /bin/rm -f conftest*]

  fi # "x$cs_have_zlib" = "xyes"

fi)dnl


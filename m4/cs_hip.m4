dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2026 EDF S.A.
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

# CS_AC_TEST_HIP
#----------------
# optional HIP support
# modifies or sets cs_have_hip, HIP_CPPFLAGS, HIP_LDFLAGS, and HIP_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_HIP], [

cs_have_hip=no
cs_have_hipblas=no
cs_have_rocsparse=no
cs_have_rccl=no

AC_ARG_ENABLE(hip,
  [AS_HELP_STRING([--enable-hip], [Enable hip offload])],
  [
    case "${enableval}" in
      yes) cs_have_hip=yes ;;
      no)  cs_have_hip=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-hip]) ;;
    esac
  ],
  [ cs_have_hip=no ]
)

AC_ARG_ENABLE(hip-amdclang,
  [AS_HELP_STRING([--enable-hip-amdclang], [Use amdclang compiler for HIP offload])],
  [
    case "${enableval}" in
      yes) cs_have_hip_amdclang=yes ;;
      no)  cs_have_hip_amdclang=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-hip-amdclang]) ;;
    esac
  ],
  [ cs_have_hip_amdclang=no ]
)

AC_LANG_PUSH(C++)

if test "x$cs_have_hip" != "xno" ; then

  # Check for hipcc compiler

  # Check first if we force amdclang++

  # if test "$CS_HIPCC" = ""; then
  if test "x$cs_have_hip_amdclang" != "xno" ; then
    AC_PATH_PROG(HIPCC, amdclang++, "no")
    AS_IF([test "x$HIPCC" = "xno"],
          [AC_MSG_ERROR([amdclang++ compiler not found!])])
  else
    AC_PATH_PROG(HIPCC, hipcc, "no")
    AS_IF([test "x$HIPCC" = "xno"],
          [AC_MSG_ERROR([HIPCC compiler not found!])])
  fi

  # Set flags, substituting "bin/hipcc" by "include".
  if test "x$cs_have_hip_amdclang" != "xno" ; then
    HIP_CPPFLAGS=" -I${HIPCC/'bin/amdclang++'/include}"

    cs_hip_lib_path="${HIPCC/'bin/amdclang++'/lib}"
  else
    HIP_CPPFLAGS=" -I${HIPCC/'bin/hipcc'/include}"

    cs_hip_lib_path="${HIPCC/'bin/hipcc'/lib}"
  fi
  AS_IF([echo $build_cpu | grep -q "_64"],
        [cs_hip_lib_path+="64"])
  HIP_LDFLAGS="-L${cs_hip_lib_path}"
  HIP_LIBS+=" -lamdhip64"

  # Try to detect available architectures.
  # As of 2024, we do not care to support HIP versions older than 11
  # (and even then,target machines should be at least Volta, though
  # developping/debugging on local machines using older hardware remains useful).

  if test "$HIP_ARCH_NUM" = ""; then
    HIP_ARCH_NUM="gfx90a"
  fi

  user_hipccflags="${HIPCCFLAGS}"

  if test "$HIP_ARCH_NUM" != ""; then
    #touch conftest.cpp
    for hip_arch in $HIP_ARCH_NUM; do
      #$HIPCC --dryrun -c conftest.cpp -o conftest.o --offload-arch=${hip_arch} >/dev/null 2>&1
      #if test $? -eq 0; then
      #  HIPCCFLAGS="${HIPCCFLAGS} --offload-arch=${hip_arch}"
      #fi
      HIPCCFLAGS="${HIPCCFLAGS} --offload-arch=${hip_arch}"
    done
    #rm -f conftest.cpp conftest.o
  fi

  HIPCCFLAGS="${HIPCCFLAGS} -x hip"
  #HIPCCFLAGS="${HIPCCFLAGS} -v"

  if test "$user_hipccflags" != ""; then
    HIPCCFLAGS="${HIPCCFLAGS} ${user_hipccflags}"
  fi

  AC_DEFINE([HAVE_HIP], 1, [HIP offload support])
  AC_ARG_VAR([HIPCCFLAGS], [Additional flags for HIP hipcc])

  AC_SUBST(cs_have_hip)
  AC_SUBST(HIPCC)
  AC_SUBST(HIPCCFLAGS)

fi

AM_CONDITIONAL([HAVE_HIP], [test "$cs_have_hip" = "yes"])

# Now check for libraries such as hipBLAS and rocsparse if HIP enabled.

cs_enable_hip_cpp=no

if test "x$cs_have_hip" != "xno" ; then

  AC_ARG_ENABLE(hip-cpp,
    [AS_HELP_STRING([--enable-hip-cpp], [Enable hip offload for .cpp files])],
    [
      case "${enableval}" in
        yes) cs_have_hip_cpp=yes ;;
        no)  cs_have_hip_cpp=no ;;
        *)   AC_MSG_ERROR([bad value ${enableval} for --enable-hip-cpp]) ;;
      esac
    ],
    [ cs_enable_hip_cpp=yes ]
  )

  AC_ARG_WITH(hipblas,
              [AS_HELP_STRING([--with-hipblas=PATH],
                              [specify prefix directory for hipblas])],
              [if test "x$withval" = "x"; then
                 with_hipblas=yes
               fi],
              [with_hipblas=check])

  AC_ARG_WITH(hipblas-include,
              [AS_HELP_STRING([--with-hipblas-include=PATH],
                              [specify directory for hipblas include files])],
              [if test "x$with_hipblas" = "xcheck"; then
                 with_hipblas=yes
               fi
               HIPBLAS_CPPFLAGS="-I$with_hipblas_include"],
              [if test "x$with_hipblas" != "xno" -a "x$with_hipblas" != "xyes" \
  	          -a "x$with_hipblas" != "xcheck"; then
                 if test "${HIPCC/'bin/hipcc'/include}" != "$with_hipblas/include" ; then
                   HIPBLAS_CPPFLAGS="-I$with_hipblas/include"
                 fi
               fi])

  AC_ARG_WITH(hipblas-lib,
              [AS_HELP_STRING([--with-hipblas-lib=PATH],
                              [specify directory for hipblas library])],
              [if test "x$with_hipblas" = "xcheck"; then
                 with_hipblas=yes
               fi
               HIPBLAS_LDFLAGS="-L$with_hipblas_lib"],
              [if test "x$with_hipblas" != "xno" -a "x$with_hipblas" != "xyes" \
	            -a "x$with_hipblas" != "xcheck"; then
                 if test "$cs_hip_lib_path" != "$with_hipblas/lib64" ; then
                   HIPBLAS_LDFLAGS="-L$with_hipblas/lib64"
                 fi
               fi])

  AC_ARG_WITH(rocsparse,
              [AS_HELP_STRING([--with-rocsparse=PATH],
                              [specify prefix directory for rocsparse])],
              [if test "x$withval" = "x"; then
                 with_rocsparse=yes
               fi],
              [with_rocsparse=check])

  AC_ARG_WITH(rocsparse-include,
              [AS_HELP_STRING([--with-rocsparse-include=PATH],
                              [specify directory for rocsparse include files])],
              [if test "x$with_rocsparse" = "xcheck"; then
                 with_rocsparse=yes
               fi
               ROCSPARSE_CPPFLAGS="-I$with_rocsparse_include"],
              [if test "x$with_rocsparse" != "xno" -a "x$with_rocsparse" != "xyes" \
  	          -a "x$with_rocsparse" != "xcheck"; then
                 if test "${HIPCC/'bin/hipcc'/include}" != "$with_rocsparse/include" ; then
                   ROCSPARSE_CPPFLAGS="-I$with_rocsparse/include"
                 fi
               fi])

  AC_ARG_WITH(rocsparse-lib,
              [AS_HELP_STRING([--with-rocsparse-lib=PATH],
                              [specify directory for rocsparse library])],
              [if test "x$with_rocsparse" = "xcheck"; then
                 with_rocsparse=yes
               fi
               ROCSPARSE_LDFLAGS="-L$with_rocsparse_lib"],
              [if test "x$with_rocsparse" != "xno" -a "x$with_rocsparse" != "xyes" \
	            -a "x$with_rocsparse" != "xcheck"; then
                 if test "$cs_hip_lib_path" != "$with_rocsparse/lib64" ; then
                   ROCSPARSE_LDFLAGS="-L$with_rocsparse/lib64"
                 fi
               fi])

  AC_ARG_WITH(rccl,
              [AS_HELP_STRING([--with-rccl=PATH],
                              [specify prefix directory for ROCm Collective Communications Library (RCCL)])],
              [if test "x$withval" = "x"; then
                 with_rccl=no
               fi],
              [with_rccl=no])

  AC_ARG_WITH(rccl-include,
              [AS_HELP_STRING([--with-rccl-include=PATH],
                              [specify directory for RCCL include files])],
              [if test "x$with_rccl" = "xcheck"; then
                 with_rccl=yes
               fi
               RCCL_CPPFLAGS="-I$with_rccl_include"],
              [if test "x$with_rccl" != "xno" -a "x$with_rccl" != "xyes" \
                  -a "x$with_rccl" != "xcheck"; then
                 if test "${HIPCC/'bin/hipcc'/include}" != "$with_rccl/include" ; then
                   RCCL_CPPFLAGS="-I$with_rccl/include"
                 fi
               fi])

  AC_ARG_WITH(rccl-lib,
              [AS_HELP_STRING([--with-rccl-lib=PATH],
                              [specify directory for RCCL library])],
              [if test "x$with_rccl" = "xcheck"; then
                 with_rccl=yes
               fi
               RCCL_LDFLAGS="-L$with_rccl_lib"],
              [if test "x$with_rccl" != "xno" -a "x$with_rccl" != "xyes" \
	            -a "x$with_rccl" != "xcheck"; then
                 if test "$cs_hip_lib_path" != "$with_rccl/lib" ; then
                   RCCL_LDFLAGS="-L$with_rccl/lib"
                 fi
               fi])

  # Check for hipBLAS

  if test "x$with_hipblas" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_HIP_CPPFLAGS="$HIP_CPPFLAGS"
    saved_HIP_LDFLAGS="$HIP_LDFLAGS"
    saved_HIP_LIBS="$HIP_LIBS"

    if test "x$HIPBLAS_CPPFLAGS" != "x" ; then
      HIP_CPPFLAGS="${HIP_CPPFLAGS} ${HIPBLAS_CPPFLAGS}"
    fi
    if test "x$HIPBLAS_LDFLAGS" != "x" ; then
      HIP_LDFLAGS="${HIP_LDFLAGS} ${HIPBLAS_LDFLAGS}"
    fi
    HIP_LIBS="-lhipblas ${HIP_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${HIP_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${HIP_LDFLAGS}"
    LIBS="${HIP_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for hipBLAS support])
    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <hipblas_v2.h>]],
[[hipblasHandle_t handle = NULL;
hipblasStatus_t status = hipblasCreate(&handle);]])
                   ],
                   [ AC_DEFINE([HAVE_HIPBLAS], 1, [hipBLAS support])
                     cs_have_hipblas=yes ],
                   [cs_have_hipblas=no])

    AC_MSG_RESULT($cs_have_hipblas)
    if test "x$cs_have_hipblas" = "xno" ; then
      if test "x$with_hipblas" != "xcheck" ; then
        AC_MSG_FAILURE([hipBLAS support is requested, but test for hipBLAS failed!])
      else
        HIP_CPPFLAGS="$saved_HIP_CPPFLAGS"
        HIP_LDFLAGS="$saved_HIP_LDFLAGS"
        HIP_LIBS="$saved_HIP_LIBS"
      fi
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi

  AC_SUBST(cs_have_hipblas)

  # Check for rocsparse

  if test "x$with_rocsparse" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_HIP_CPPFLAGS="$HIP_CPPFLAGS"
    saved_HIP_LDFLAGS="$HIP_LDFLAGS"
    saved_HIP_LIBS="$HIP_LIBS"

    if test "x$ROCSPARSE_CPPFLAGS" != "x" ; then
      HIP_CPPFLAGS="${HIP_CPPFLAGS} ${ROCSPARSE_CPPFLAGS}"
    fi
    if test "x$ROCSPARSE_LDFLAGS" != "x" ; then
      HIP_LDFLAGS="${HIP_LDFLAGS} ${ROCSPARSE_LDFLAGS}"
    fi
    HIP_LIBS="-lrocsparse ${HIP_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${HIP_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${HIP_LDFLAGS}"
    LIBS="${HIP_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for rocsparse support])

    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <rocsparse/rocsparse.h>]],
[[rocsparse_spmat_descr mat_a;
rocsparse_status status = rocsparse_destroy_spmat_descr(mat_a);]])
                   ],
                   [ AC_DEFINE([HAVE_ROCSPARSE], 1,
                               [rocsparse support ])
                     AC_MSG_RESULT([rocsparse found])
                     cs_have_rocsparse=yes ],
                   [cs_have_rocsparse=no])

    AC_MSG_RESULT($cs_have_rocsparse)
    if test "x$cs_have_rocsparse" = "xno" ; then
      if test "x$with_rocsparse" != "xcheck" ; then
        AC_MSG_FAILURE([rocsparse support is requested, but test for rocsparse failed!])
      else
        HIP_CPPFLAGS="$saved_HIP_CPPFLAGS"
        HIP_LDFLAGS="$saved_HIP_LDFLAGS"
        HIP_LIBS="$saved_HIP_LIBS"
      fi
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi

  AC_SUBST(cs_have_rocsparse)

  # Check for RCCL

  if test "x$with_rccl" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_HIP_CPPFLAGS="$HIP_CPPFLAGS"
    saved_HIP_LDFLAGS="$HIP_LDFLAGS"
    saved_HIP_LIBS="$HIP_LIBS"

    if test "x$RCCL_CPPFLAGS" != "x" ; then
      HIP_CPPFLAGS="${HIP_CPPFLAGS} ${RCCL_CPPFLAGS}"
    fi
    if test "x$RCCL_LDFLAGS" != "x" ; then
      HIP_LDFLAGS="${HIP_LDFLAGS} ${RCCL_LDFLAGS}"
    fi
    HIP_LIBS="-lrccl ${HIP_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${HIP_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${HIP_LDFLAGS}"
    LIBS="${HIP_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for RCCL support])

    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <rccl.h>]],
[[rcclComm_t* comm;
rcclUniqueId Id;
rcclCommInitRank(comm, 1, Id, 0);]])
                 ],
                 [ cs_have_rccl=yes ],
                 [ cs_have_rccl=no ])

    if test "$cs_have_rccl" = "yes"; then
      AC_DEFINE([HAVE_RCCL], 1, [RCCL support])
    fi

    AC_MSG_RESULT($cs_have_rccl)
    if test "x$cs_have_rccl" = "xno" ; then
      if test "x$with_rccl" != "xcheck" ; then
        AC_MSG_FAILURE([RCCL support is requested, but test for RCCL failed!])
      else
        HIP_CPPFLAGS="$saved_HIP_CPPFLAGS"
        HIP_LDFLAGS="$saved_HIP_LDFLAGS"
        HIP_LIBS="$saved_HIP_LIBS"
      fi
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi

  AC_SUBST(cs_have_rccl)

  # Finally set flags which can be extended by libraries and paths.

  AC_SUBST(HIP_CPPFLAGS)
  AC_SUBST(HIP_LDFLAGS)
  AC_SUBST(HIP_LIBS)

fi

AC_LANG_POP(C++)

AM_CONDITIONAL([HAVE_HIP_CPP], [test "$cs_enable_hip_cpp" = "yes"])

])dnl

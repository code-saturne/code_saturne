dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2024 EDF S.A.
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
cs_have_hipsparse=no
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

if test "x$cs_have_hip" != "xno" ; then

  # Check for hipcc compiler

  AC_PATH_PROG(HIPCC, hipcc, "no")
  AS_IF([test "x$HIPCC" = "xno"],
        [AC_MSG_ERROR([HIPCC compiler not found!])])

  # Set flags, substituting "bin/hipcc" by "include".
  HIP_CPPFLAGS=" -I${HIPCC/'bin/hipcc'/include}"

  cs_hip_lib_path="${HIPCC/'bin/hipcc'/lib}"
  AS_IF([echo $build_cpu | grep -q "_64"],
        [cs_hip_lib_path+="64"])
  HIP_LDFLAGS="-L${cs_hip_lib_path}"
  HIP_LIBS+=" -lamdhip64"

  # Try to detect available architectures.
  # As of 2024, we do not care to support HIP versions older than 11
  # (and even then,target machines should be at least Volta, though
  # developping/debugging on local machines using older hardware remains useful).

  if test "$HIP_ARCH_NUM" = ""; then
    # HIP_ARCH_NUM="60 61 62 70 72 75 80 86"
    HIP_ARCH_NUM="70 80"
  fi

  user_hipccflags="${HIPCCFLAGS}"

  if test "$HIP_ARCH_NUM" != ""; then
    touch conftest.cpp
    for hip_arch in $HIP_ARCH_NUM; do
      $HIPCC --dryrun -c conftest.cpp -o conftest.o -gencode arch=compute_${hip_arch},code=sm_${hip_arch} >/dev/null 2>&1
      if test $? -eq 0; then
        HIPCCFLAGS="${HIPCCFLAGS} -gencode arch=compute_${hip_arch},code=sm_${hip_arch}"
      fi
    done
    rm -f conftest.cpp conftest.o
  fi

  HIPCCFLAGS="${HIPCCFLAGS} -v"

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

# Now check for libraries such as hipBLAS and hipSPARSE if HIP enabled.

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

  AC_ARG_WITH(hipsparse,
              [AS_HELP_STRING([--with-hipsparse=PATH],
                              [specify prefix directory for hipSPARSE])],
              [if test "x$withval" = "x"; then
                 with_hipsparse=yes
               fi],
              [with_hipsparse=check])

  AC_ARG_WITH(hipsparse-include,
              [AS_HELP_STRING([--with-hipsparse-include=PATH],
                              [specify directory for hipSPARSE include files])],
              [if test "x$with_hipsparse" = "xcheck"; then
                 with_hipsparse=yes
               fi
               HIPSPARSE_CPPFLAGS="-I$with_hipsparse_include"],
              [if test "x$with_hipsparse" != "xno" -a "x$with_hipsparse" != "xyes" \
  	          -a "x$with_hipsparse" != "xcheck"; then
                 if test "${HIPCC/'bin/hipcc'/include}" != "$with_hipsparse/include" ; then
                   HIPSPARSE_CPPFLAGS="-I$with_hipsparse/include"
                 fi
               fi])

  AC_ARG_WITH(hipsparse-lib,
              [AS_HELP_STRING([--with-hipsparse-lib=PATH],
                              [specify directory for hipSPARSE library])],
              [if test "x$with_hipsparse" = "xcheck"; then
                 with_hipsparse=yes
               fi
               HIPSPARSE_LDFLAGS="-L$with_hipsparse_lib"],
              [if test "x$with_hipsparse" != "xno" -a "x$with_hipsparse" != "xyes" \
	            -a "x$with_hipsparse" != "xcheck"; then
                 if test "$cs_hip_lib_path" != "$with_hipsparse/lib64" ; then
                   HIPSPARSE_LDFLAGS="-L$with_hipsparse/lib64"
                 fi
               fi])

  AC_ARG_WITH(nccl,
              [AS_HELP_STRING([--with-nccl=PATH],
                              [specify prefix directory for NVIDIA Collective Communications Library (NCCL)])],
              [if test "x$withval" = "x"; then
                 with_nccl=no
               fi],
              [with_nccl=no])

  AC_ARG_WITH(nccl-include,
              [AS_HELP_STRING([--with-nccl-include=PATH],
                              [specify directory for NCCL include files])],
              [if test "x$with_nccl" = "xcheck"; then
                 with_nccl=yes
               fi
               NCCL_CPPFLAGS="-I$with_nccl_include"],
              [if test "x$with_nccl" != "xno" -a "x$with_nccl" != "xyes" \
                  -a "x$with_nccl" != "xcheck"; then
                 if test "${HIPCC/'bin/hipcc'/include}" != "$with_nccl/include" ; then
                   NCCL_CPPFLAGS="-I$with_nccl/include"
                 fi
               fi])

  AC_ARG_WITH(nccl-lib,
              [AS_HELP_STRING([--with-nccl-lib=PATH],
                              [specify directory for NCCL library])],
              [if test "x$with_nccl" = "xcheck"; then
                 with_nccl=yes
               fi
               NCCL_LDFLAGS="-L$with_nccl_lib"],
              [if test "x$with_nccl" != "xno" -a "x$with_nccl" != "xyes" \
	            -a "x$with_nccl" != "xcheck"; then
                 if test "$cs_hip_lib_path" != "$with_nccl/lib" ; then
                   NCCL_LDFLAGS="-L$with_nccl/lib"
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

  # Check for hipSPARSE

  if test "x$with_hipsparse" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_HIP_CPPFLAGS="$HIP_CPPFLAGS"
    saved_HIP_LDFLAGS="$HIP_LDFLAGS"
    saved_HIP_LIBS="$HIP_LIBS"

    if test "x$HIPSPARSE_CPPFLAGS" != "x" ; then
      HIP_CPPFLAGS="${HIP_CPPFLAGS} ${HIPSPARSE_CPPFLAGS}"
    fi
    if test "x$HIPSPARSE_LDFLAGS" != "x" ; then
      HIP_LDFLAGS="${HIP_LDFLAGS} ${HIPSPARSE_LDFLAGS}"
    fi
    HIP_LIBS="-lhipsparse ${HIP_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${HIP_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${HIP_LDFLAGS}"
    LIBS="${HIP_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for hipSPARSE support])

    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <hipsparse.h>]],
[[hipsparseSpMatDescr_t matA;
hipsparseStatus_t status = hipsparseDestroySpMat(matA);]])
                   ],
                   [ AC_DEFINE([HAVE_HIPSPARSE], 1,
                               [hipSPARSE support ])
                     AC_MSG_RESULT([hipSPARSE found])
                     cs_have_hipsparse=yes ],
                   [cs_have_hipsparse=no])

    AC_MSG_RESULT($cs_have_hipsparse)
    if test "x$cs_have_hipsparse" = "xno" ; then
      if test "x$with_hipsparse" != "xcheck" ; then
        AC_MSG_FAILURE([hipSPARSE support is requested, but test for hipSPARSE failed!])
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

  AC_SUBST(cs_have_hipsparse)

  # Check for NCCL

  if test "x$with_nccl" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_HIP_CPPFLAGS="$HIP_CPPFLAGS"
    saved_HIP_LDFLAGS="$HIP_LDFLAGS"
    saved_HIP_LIBS="$HIP_LIBS"

    if test "x$NCCL_CPPFLAGS" != "x" ; then
      HIP_CPPFLAGS="${HIP_CPPFLAGS} ${NCCL_CPPFLAGS}"
    fi
    if test "x$NCCL_LDFLAGS" != "x" ; then
      HIP_LDFLAGS="${HIP_LDFLAGS} ${NCCL_LDFLAGS}"
    fi
    HIP_LIBS="-lnccl ${HIP_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${HIP_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${HIP_LDFLAGS}"
    LIBS="${HIP_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for NCCL support])

    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <nccl.h>]],
[[ncclComm_t* comm;
ncclUniqueId Id;
ncclCommInitRank(comm, 1, Id, 0);]])
                 ],
                 [ cs_have_nccl=yes ],
                 [ cs_have_nccl=no ])

    if test "$cs_have_nccl" = "yes"; then
      AC_DEFINE([HAVE_NCCL], 1, [NCCL support])
    fi

    AC_MSG_RESULT($cs_have_nccl)
    if test "x$cs_have_nccl" = "xno" ; then
      if test "x$with_nccl" != "xcheck" ; then
        AC_MSG_FAILURE([NCCL support is requested, but test for NCCL failed!])
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

  AC_SUBST(cs_have_nccl)

  # Finally set flags which can be extended by libraries and paths.

  AC_SUBST(HIP_CPPFLAGS)
  AC_SUBST(HIP_LDFLAGS)
  AC_SUBST(HIP_LIBS)

fi

AM_CONDITIONAL([HAVE_HIP_CPP], [test "$cs_enable_hip_cpp" = "yes"])

])dnl

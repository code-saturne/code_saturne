dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of code_saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2022 EDF S.A.
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

# CS_AC_TEST_CUDA
#----------------
# optional CUDA support
# modifies or sets cs_have_cuda, CUDA_CPPFLAGS, CUDA_LDFLAGS, and CUDA_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CUDA], [

cs_have_cuda=no
cs_have_cublas=no
cs_have_cusparse=no

AC_ARG_ENABLE(cuda,
  [AS_HELP_STRING([--enable-cuda], [Enable cuda offload])],
  [
    case "${enableval}" in
      yes) cs_have_cuda=yes ;;
      no)  cs_have_cuda=no ;;
      *)   AC_MSG_ERROR([bad value ${enableval} for --enable-cuda]) ;;
    esac
  ],
  [ cs_have_cuda=no ]
)

if test "x$cs_have_cuda" != "xno" ; then

  # Check for nvcc compiler

  AC_PATH_PROG(NVCC, nvcc, "no")
  AS_IF([test "x$NVCC" = "xno"],
        [AC_MSG_ERROR([NVCC compiler not found!])])

  # Set flags, substituting "bin/nvcc" by "include".
  CUDA_CPPFLAGS=" -I${NVCC/'bin/nvcc'/include}"

  cs_cuda_lib_path="${NVCC/'bin/nvcc'/lib}"
  AS_IF([echo $build_cpu | grep -q "_64"],
        [cs_cuda_lib_path+="64"])
  CUDA_LDFLAGS="-L${cs_cuda_lib_path}"
  CUDA_LIBS+=" -lcudart"

  # Try to detect available architectures.
  # As of late 2021, we do not care to support CUDA versions older than 9
  # (and even then,target machines should be at least Volta,
  # though developping/debugging on local machines using CUDA 9 remains useful).

  if test "$CUDA_ARCH_NUM" = ""; then
    # CUDA_ARCH_NUM="60 61 62 70 72 75 80 86"
    CUDA_ARCH_NUM="60 70 80"
  fi

  NVCCFLAGS="-ccbin $CXX -DHAVE_CONFIG_H"  # wrap C++ compiler arount nvcc
  if test "$CUDA_ARCH_NUM" != ""; then
    touch conftest.cu
    for cu_arch in $CUDA_ARCH_NUM; do
      $NVCC --dryrun -c conftest.cu -o conftest.o -gencode arch=compute_${cu_arch},code=sm_${cu_arch} >/dev/null 2>&1
      if test $? -eq 0; then
        NVCCFLAGS="${NVCCFLAGS} -gencode arch=compute_${cu_arch},code=sm_${cu_arch}"
      fi
    done
    rm -f conftest.cu conftest.o
  fi

  NVCCFLAGS="${NVCCFLAGS} --maxrregcount=64 -Xptxas -v"

  AC_DEFINE([HAVE_CUDA], 1, [CUDA offload support])

  AC_SUBST(cs_have_cuda)
  AC_SUBST(NVCC)
  AC_SUBST(NVCCFLAGS)

fi

AM_CONDITIONAL([HAVE_CUDA], [test "$cs_have_cuda" = "yes"])

# Now check for libraries such as cuBLAS and cuSPARSE if CUDA enabled.

if test "x$cs_have_cuda" != "xno" ; then

  AC_ARG_WITH(cublas,
              [AS_HELP_STRING([--with-cublas=PATH],
                              [specify prefix directory for cublas])],
              [if test "x$withval" = "x"; then
                 with_cublas=yes
               fi],
              [with_cublas=check])

  AC_ARG_WITH(cublas-include,
              [AS_HELP_STRING([--with-cublas-include=PATH],
                              [specify directory for cublas include files])],
              [if test "x$with_cublas" = "xcheck"; then
                 with_cublas=yes
               fi
               CUBLAS_CPPFLAGS="-I$with_cublas_include"],
              [if test "x$with_cublas" != "xno" -a "x$with_cublas" != "xyes" \
  	          -a "x$with_cublas" != "xcheck"; then
                 if test "${NVCC/'bin/nvcc'/include}" != "$with_cublas/include" ; then
                   CUBLAS_CPPFLAGS="-I$with_cublas/include"
                 fi
               fi])

  AC_ARG_WITH(cublas-lib,
              [AS_HELP_STRING([--with-cublas-lib=PATH],
                              [specify directory for cublas library])],
              [if test "x$with_cublas" = "xcheck"; then
                 with_cublas=yes
               fi
               CUBLAS_LDFLAGS="-L$with_cublas_lib"],
              [if test "x$with_cublas" != "xno" -a "x$with_cublas" != "xyes" \
	            -a "x$with_cublas" != "xcheck"; then
                 if test "$cs_cuda_lib_path" != "$with_cublas/lib64" ; then
                   CUBLAS_LDFLAGS="-L$with_cublas/lib64"
                 fi
               fi])

  AC_ARG_WITH(cusparse,
              [AS_HELP_STRING([--with-cusparse=PATH],
                              [specify prefix directory for cuSPARSE])],
              [if test "x$withval" = "x"; then
                 with_cusparse=yes
               fi],
              [with_cusparse=check])

  AC_ARG_WITH(cusparse-include,
              [AS_HELP_STRING([--with-cusparse-include=PATH],
                              [specify directory for cuSPARSE include files])],
              [if test "x$with_cusparse" = "xcheck"; then
                 with_cusparse=yes
               fi
               CUSPARSE_CPPFLAGS="-I$with_cusparse_include"],
              [if test "x$with_cusparse" != "xno" -a "x$with_cusparse" != "xyes" \
  	          -a "x$with_cusparse" != "xcheck"; then
                 if test "${NVCC/'bin/nvcc'/include}" != "$with_cusparse/include" ; then
                   CUSPARSE_CPPFLAGS="-I$with_cusparse/include"
                 fi
               fi])

  AC_ARG_WITH(cusparse-lib,
              [AS_HELP_STRING([--with-cusparse-lib=PATH],
                              [specify directory for cuSPARSE library])],
              [if test "x$with_cusparse" = "xcheck"; then
                 with_cusparse=yes
               fi
               CUSPARSE_LDFLAGS="-L$with_cusparse_lib"],
              [if test "x$with_cusparse" != "xno" -a "x$with_cusparse" != "xyes" \
	            -a "x$with_cusparse" != "xcheck"; then
                 if test "$cs_cuda_lib_path" != "$with_cusparse/lib64" ; then
                   CUSPARSE_LDFLAGS="-L$with_cusparse/lib64"
                 fi
               fi])

  # Check for cuBLAS

  if test "x$with_cublas" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_CUDA_CPPFLAGS="$CUDA_CPPFLAGS"
    saved_CUDA_LDFLAGS="$CUDA_LDFLAGS"
    saved_CUDA_LIBS="$CUDA_LIBS"

    if test "x$CUBLAS_CPPFLAGS" != "x" ; then
      CUDA_CPPFLAGS="${CUDA_CPPFLAGS} ${CUBLAS_CPPFLAGS}"
    fi
    if test "x$CUBLAS_LDFLAGS" != "x" ; then
      CUDA_LDFLAGS="${CUDA_LDFLAGS} ${CUBLAS_LDFLAGS}"
    fi
    CUDA_LIBS="-lcublas ${CUDA_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${CUDA_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${CUDA_LDFLAGS}"
    LIBS="${CUDA_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for cuBLAS support])
    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <cublas_v2.h>]],
[[cublasHandle_t handle = NULL;
cublasStatus_t status = cublasCreate(&handle);]])
                   ],
                   [ AC_DEFINE([HAVE_CUBLAS], 1, [cuBLAS support])
                     cs_have_cublas=yes ],
                   [cs_have_cublas=no])

    AC_MSG_RESULT($cs_have_cublas)
    if test "x$cs_have_cublas" = "xno" ; then
      if test "x$with_cublas" != "xcheck" ; then
        AC_MSG_FAILURE([cuBLAS support is requested, but test for cuBLAS failed!])
      else
        CUDA_CPPFLAGS="$saved_CUDA_CPPFLAGS"
        CUDA_LDFLAGS="$saved_CUDA_LDFLAGS"
        CUDA_LIBS="$saved_CUDA_LIBS"
      fi
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi

  AC_SUBST(cs_have_cublas)

  # Check for cuSPARSE

  if test "x$with_cusparse" != "xno" ; then

    saved_CPPFLAGS="$CPPFLAGS"
    saved_LDFLAGS="$LDFLAGS"
    saved_LIBS="$LIBS"

    saved_CUDA_CPPFLAGS="$CUDA_CPPFLAGS"
    saved_CUDA_LDFLAGS="$CUDA_LDFLAGS"
    saved_CUDA_LIBS="$CUDA_LIBS"

    if test "x$CUSPARSE_CPPFLAGS" != "x" ; then
      CUDA_CPPFLAGS="${CUDA_CPPFLAGS} ${CUSPARSE_CPPFLAGS}"
    fi
    if test "x$CUSPARSE_LDFLAGS" != "x" ; then
      CUDA_LDFLAGS="${CUDA_LDFLAGS} ${CUSPARSE_LDFLAGS}"
    fi
    CUDA_LIBS="-lcusparse ${CUDA_LIBS}"

    CPPFLAGS="${CPPFLAGS} ${CUDA_CPPFLAGS}"
    LDFLAGS="${LDFLAGS} ${CUDA_LDFLAGS}"
    LIBS="${CUDA_LIBS} ${LIBS}"

    AC_MSG_CHECKING([for cuSPARSE support])

    AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <cusparse.h>]],
[[cusparseSpMatDescr_t matA;
cusparseStatus_t status = cusparseDestroySpMat(matA);]])
                   ],
                   [ AC_DEFINE([HAVE_CUSPARSE_GENERIC_API], 1,
                               [cuSPARSE generic API support ])
                     AC_MSG_RESULT([cuSPARSE generic API found])
                     cs_have_cusparse=yes ],
                   [cs_have_cusparse=fallback])

    if test "x$cs_have_cusparse" = "xfallback" ; then

      AC_LINK_IFELSE(
[AC_LANG_PROGRAM([[#include <cusparse.h>]],
[[cusparseHandle_t handle = NULL;
cusparseStatus_t status = cusparseCreate(&handle);]])
                   ],
                   [ cs_have_cusparse=yes ],
                   [ cs_have_cusparse=no ])

    fi

    if test "$cs_have_cusparse" = "yes"; then
      AC_DEFINE([HAVE_CUSPARSE], 1, [cuSPARSE support])
    fi

    AC_MSG_RESULT($cs_have_cusparse)
    if test "x$cs_have_cusparse" = "xno" ; then
      if test "x$with_cusparse" != "xcheck" ; then
        AC_MSG_FAILURE([cuSPARSE support is requested, but test for cuSPARSE failed!])
      else
        CUDA_CPPFLAGS="$saved_CUDA_CPPFLAGS"
        CUDA_LDFLAGS="$saved_CUDA_LDFLAGS"
        CUDA_LIBS="$saved_CUDA_LIBS"
      fi
    fi

    CPPFLAGS="$saved_CPPFLAGS"
    LDFLAGS="$saved_LDFLAGS"
    LIBS="$saved_LIBS"

  fi

  AC_SUBST(cs_have_cusparse)

  # Finally set flags which can be extended by libraries and paths.

  AC_SUBST(CUDA_CPPFLAGS)
  AC_SUBST(CUDA_LDFLAGS)
  AC_SUBST(CUDA_LIBS)

fi

])dnl

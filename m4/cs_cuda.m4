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

# CS_AC_TEST_CUDA
#----------------
# optional CUDA support
# modifies or sets cs_have_cuda, CUDA_CPPFLAGS, CUDA_LDFLAGS, and CUDA_LIBS
# depending on libraries found

AC_DEFUN([CS_AC_TEST_CUDA], [

cs_have_cuda=no

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

  CUDA_LDFLAGS=""
  CUDA_LIBS=" -L${NVCC/'bin/nvcc'/lib}"
  AS_IF([echo $build_cpu | grep -q "_64"],
        [CUDA_LIBS+="64"])
  CUDA_LIBS+=" -lcudart"

  # Try to detect available architectures

  if test "$CUDA_ARCH_NUM" = ""; then
    CUDA_ARCH_NUM="35 37 60 70"
  fi

  NVCC_FLAGS="-ccbin $CXX -DHAVE_CONFIG_H"  # wrap C++ compiler arount nvcc
  touch conftest.cu
  for cu_arch in "$CUDA_ARCH_NUM"; do
    $NVCC --dryrun -c conftest.cu -o conftest.o -gencode arch=compute_${cu_arch},code=sm_${cu_arch} >/dev/null 2>&1
    if test $? -eq 0; then
      NVCC_FLAGS="${NVCC_FLAGS} -gencode arch=compute_${cu_arch},code=sm_${cu_arch}"
    fi
  done
  rm -f conftest.cu conftest.o

  NVCC_FLAGS="${NVCC_FLAGS} --maxrregcount=64 -Xptxas -v"
  if test "x$enable_shared" = "xyes" ; then
    NVCC_FLAGS="${NVCC_FLAGS} --compiler-options -fPIC"
  fi

  AC_DEFINE([HAVE_CUDA], 1, [CUDA offload support])

  AC_SUBST(cs_have_cuda)
  AC_SUBST(CUDA_CPPFLAGS)
  AC_SUBST(CUDA_LDFLAGS)
  AC_SUBST(CUDA_LIBS)
  AC_SUBST(NVCC)
  AC_SUBST(NVCC_FLAGS)

fi

AM_CONDITIONAL([HAVE_CUDA], [test "$cs_have_cuda" = "yes"])

])dnl


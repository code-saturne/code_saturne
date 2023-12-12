## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
  LDFLAGS_PLE = -L"$(top_builddir)/libple/src" ${LDRPATH}"$(abs_top_builddir)/libple/src"
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
  LDFLAGS_PLE = $(PLE_LDFLAGS) $(PLERUNPATH)
endif
LIBS_PLE = $(PLE_LIBS)

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS   = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS  =

BUILT_SOURCES = \
cs_halo.c \
cs_range_set.c \
cs_sort.c \
cs_matrix.c \
cs_matrix_assembler.c \
cs_blas.c \
cs_random.c

cs_halo.c: Makefile $(top_srcdir)/src/base/cs_halo.c
	cat $(top_srcdir)/src/base/$@ >$@

cs_range_set.c: Makefile $(top_srcdir)/src/base/cs_range_set.c
	cat $(top_srcdir)/src/base/$@ >$@

cs_sort.c: Makefile $(top_srcdir)/src/base/cs_sort.c
	cat $(top_srcdir)/src/base/$@ >$@

cs_random.c: Makefile $(top_srcdir)/src/base/cs_random.c
	cat $(top_srcdir)/src/base/$@ >$@

cs_blas.c: Makefile $(top_srcdir)/src/alge/cs_blas.c
	cat $(top_srcdir)/src/alge/$@ >$@

cs_matrix.c: Makefile $(top_srcdir)/src/alge/cs_matrix.c
	cat $(top_srcdir)/src/alge/$@ >$@

cs_matrix_assembler.c: Makefile $(top_srcdir)/src/alge/cs_matrix_assembler.c
	cat $(top_srcdir)/src/alge/$@ >$@

check_PROGRAMS =

# BFT tests

check_PROGRAMS += \
bft_backtrace_test \
bft_error_test \
bft_mem_usage_test \
bft_mem_test \
bft_printf_test

CPPFLAGS_BFT_TESTS = -I$(top_srcdir)/src/bft
LDFLAGS_BFT_TESTS =
LDADD_BFT_TESTS = $(top_builddir)/src/bft/libbft.a

bft_backtrace_test_SOURCES  = bft_backtrace_test.c
bft_backtrace_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_backtrace_test_LDADD    = $(LDADD_BFT_TESTS)

bft_error_test_SOURCES  = bft_error_test.c
bft_error_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_error_test_LDADD    = $(LDADD_BFT_TESTS)

bft_mem_usage_test_SOURCES  = bft_mem_usage_test.c
bft_mem_usage_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_mem_usage_test_LDADD    = $(LDADD_BFT_TESTS)

bft_mem_test_SOURCES  = bft_mem_test.c
bft_mem_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_mem_test_LDADD    = $(LDADD_BFT_TESTS)

bft_printf_test_SOURCES  = bft_printf_test.c
bft_printf_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_printf_test_LDADD    = $(LDADD_BFT_TESTS)

# code_saturne tests

check_PROGRAMS += \
cs_all_to_all_test \
cs_blas_test \
cs_check_cdo \
cs_check_quadrature \
cs_check_sdm \
cs_core_test \
cs_file_test \
cs_geom_test \
cs_interface_test \
cs_map_test \
cs_matrix_test \
cs_moment_test \
cs_random_test \
cs_rank_neighbors_test \
fvm_selector_test \
fvm_selector_postfix_test \
cs_sizes_test \
cs_tree_test

if HAVE_ACCEL
check_PROGRAMS += cs_gpu_test
endif

LDFLAGS_CS_TESTS = $(CGNS_LDFLAGS) $(MED_LDFLAGS) $(HDF5_LDFLAGS) \
	$(LDFLAGS_PLE) $(MPI_LDFLAGS) $(CUDA_LDFLAGS)

LDADD_CS_TESTS = \
	$(top_builddir)/src/base/libcscore.a \
	$(top_builddir)/src/bft/libbft.a \
	$(PLE_LIBS) $(MPI_LIBS) $(CUDA_LIBS) -lm

cs_all_to_all_test_SOURCES  = cs_all_to_all_test.c
cs_all_to_all_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_all_to_all_test_LDADD    = $(LDADD_CS_TESTS)

cs_blas_test_SOURCES  = \
cs_blas_test.c \
cs_blas.c
cs_blas_test_CPPFLAGS  = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
 $(BLAS_CPPFLAGS) $(MPI_CPPFLAGS)
cs_blas_test_CFLAGS  = $(AM_CFLAGS) $(CFLAGS_HOT)
cs_blas_test_LDFLAGS  = $(BLAS_LDFLAGS) $(LDFLAGS_CS_TESTS)
cs_blas_test_LDADD    = $(BLASRUNPATH) $(BLAS_LIBS) $(LDADD_CS_TESTS)

cs_check_cdo$(EXEEXT):
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_check_cdo $(top_srcdir)/tests/cs_check_cdo.c

cs_check_quadrature$(EXEEXT):
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_check_quadrature $(top_srcdir)/tests/cs_check_quadrature.c

cs_check_sdm$(EXEEXT):
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_check_sdm $(top_srcdir)/tests/cs_check_sdm.c

cs_core_test_SOURCES  = cs_core_test.c
cs_core_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_core_test_LDADD    = $(LDADD_CS_TESTS)

cs_file_test_SOURCES  = cs_file_test.c
cs_file_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_file_test_LDADD    = $(LDADD_CS_TESTS)

cs_geom_test_SOURCES  = cs_geom_test.c
cs_geom_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_geom_test_LDADD    = $(LDADD_CS_TESTS)

if HAVE_ACCEL

cs_gpu_test_SOURCES = $(top_srcdir)/tests/cs_gpu_test.cxx
if HAVE_CUDA
  cs_gpu_test_SOURCES += $(top_srcdir)/tests/cs_gpu_cuda_test.cu
endif

cs_gpu_test$(EXEEXT): $(top_srcdir)/tests/cs_gpu_test.cxx $(top_srcdir)/tests/cs_gpu_cuda_test.cu
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_gpu_test \
	${cs_gpu_test_SOURCES}

endif

cs_interface_test_SOURCES  = cs_interface_test.c
cs_interface_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_interface_test_LDADD    = \
$(top_builddir)/src/fvm/libfvm.a \
$(LDADD_CS_TESTS)

cs_map_test_SOURCES  = cs_map_test.c
cs_map_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_map_test_LDADD    = $(LDADD_CS_TESTS)

cs_matrix_test$(EXEEXT): $(top_srcdir)/tests/cs_matrix_test.c
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_matrix_test $(top_srcdir)/tests/cs_matrix_test.c

cs_moment_test_SOURCES  = cs_moment_test.c
cs_moment_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_moment_test_LDADD    = $(LDADD_CS_TESTS)

cs_random_test_SOURCES  = \
cs_random_test.c \
cs_random.c
cs_random_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_random_test_LDADD    = $(LDADD_CS_TESTS)

cs_rank_neighbors_test_SOURCES  = cs_rank_neighbors_test.c
cs_rank_neighbors_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_rank_neighbors_test_LDADD    = $(LDADD_CS_TESTS)

fvm_selector_test$(EXEEXT):
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o fvm_selector_postfix_test \
	$(top_srcdir)/tests/fvm_selector_test.c

fvm_selector_postfix_test$(EXEEXT):
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o fvm_selector_postfix_test \
	$(top_srcdir)/tests/fvm_selector_postfix_test.c

cs_sizes_test_SOURCES  = cs_sizes_test.c
cs_sizes_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_sizes_test_LDADD    = $(LDADD_CS_TESTS)

cs_tree_test_SOURCES  = cs_tree_test.c
cs_tree_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_tree_test_LDADD    = $(LDADD_CS_TESTS)

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)

# Distribution
# Files not built through classical Automake rules need to be explicitely
# added to distribution.

EXTRA_DIST = \
cs_gpu_cuda_test.cu \
unittests.py \
$(top_srcdir)/tests/graphics

# Clean

clean-local:
	-rm -f $(BUILT_SOURCES)

.cu.$(OBJEXT):
	$(NVCC) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) -c -o $@ $<

## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
AM_CXXFLAGS = $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CFLAGS_OPT)

AM_LDFLAGS  = $(CXXFLAGS_SYCL)

BUILT_SOURCES = \
cs_halo.cpp \
cs_range_set.cpp \
cs_sort.cpp \
cs_matrix.cpp \
cs_matrix_assembler.cpp \
cs_blas.cpp \
cs_random.cpp

cs_halo.cpp: Makefile $(top_srcdir)/src/base/cs_halo.cpp
	cat $(top_srcdir)/src/base/$@ >$@

cs_range_set.cpp: Makefile $(top_srcdir)/src/base/cs_range_set.cpp
	cat $(top_srcdir)/src/base/$@ >$@

cs_sort.cpp: Makefile $(top_srcdir)/src/base/cs_sort.cpp
	cat $(top_srcdir)/src/base/$@ >$@

cs_random.cpp: Makefile $(top_srcdir)/src/base/cs_random.cpp
	cat $(top_srcdir)/src/base/$@ >$@

cs_blas.cpp: Makefile $(top_srcdir)/src/alge/cs_blas.cpp
	cat $(top_srcdir)/src/alge/$@ >$@

cs_matrix.cpp: Makefile $(top_srcdir)/src/alge/cs_matrix.cpp
	cat $(top_srcdir)/src/alge/$@ >$@

cs_matrix_assembler.cpp: Makefile $(top_srcdir)/src/alge/cs_matrix_assembler.cpp
	cat $(top_srcdir)/src/alge/$@ >$@

check_PROGRAMS =

# BFT tests

check_PROGRAMS += \
bft_backtrace_test \
bft_error_test \
bft_mem_usage_test \
bft_mem_test \
bft_printf_test

bft_backtrace_test_SOURCES  = bft_backtrace_test.c
bft_backtrace_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
bft_backtrace_test_LDADD    = $(LDADD_CS_TESTS)

bft_error_test_SOURCES  = bft_error_test.c
bft_error_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
bft_error_test_LDADD    = $(LDADD_CS_TESTS)

bft_mem_usage_test_SOURCES  = bft_mem_usage_test.c
bft_mem_usage_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
bft_mem_usage_test_LDADD    = $(LDADD_CS_TESTS)

bft_mem_test_SOURCES  = bft_mem_test.cpp
bft_mem_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
bft_mem_test_LDADD    = $(LDADD_CS_TESTS)

bft_printf_test_SOURCES  = bft_printf_test.c
bft_printf_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
bft_printf_test_LDADD    = $(LDADD_CS_TESTS)

# code_saturne tests

check_PROGRAMS += \
cs_all_to_all_test \
cs_blas_test \
cs_check_cdo \
cs_check_quadrature \
cs_check_sdm \
cs_core_test \
cs_dispatch_test \
cs_file_test \
cs_geom_test \
cs_interface_test \
cs_map_test \
cs_matrix_test \
cs_mesh_quantities_test \
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

cs_all_to_all_test_SOURCES  = cs_all_to_all_test.cpp
cs_all_to_all_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_all_to_all_test_LDADD    = $(LDADD_CS_TESTS)

cs_blas_test_SOURCES  = \
cs_blas_test.cpp \
cs_blas.cpp
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

cs_core_test_SOURCES  = cs_core_test.cpp
cs_core_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_core_test_LDADD    = $(LDADD_CS_TESTS)

cs_dispatch_test_SOURCES = cs_dispatch_test.cpp

cs_dispatch_test$(EXEEXT): cs_dispatch_test.cpp $(top_srcdir)/src/base/cs_dispatch.h
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_dispatch_test $(top_srcdir)/tests/cs_dispatch_test.cpp

cs_file_test_SOURCES  = cs_file_test.cpp
cs_file_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_file_test_LDADD    = $(LDADD_CS_TESTS)

cs_geom_test_SOURCES  = cs_geom_test.cpp
cs_geom_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_geom_test_LDADD    = $(LDADD_CS_TESTS)

if HAVE_ACCEL

cs_gpu_test_SOURCES = cs_gpu_test.cpp
if HAVE_CUDA
  cs_gpu_test_SOURCES += cs_gpu_cuda_test.cu
endif

cs_gpu_test$(EXEEXT): cs_gpu_test.cpp cs_gpu_cuda_test.cu
	@list=; for p in ${cs_gpu_test_SOURCES}; do \
	  list="$$list $(top_srcdir)/tests/$$p"; \
	done; \
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_gpu_test $$list

endif

cs_interface_test_SOURCES  = cs_interface_test.cpp
cs_interface_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_interface_test_LDADD    = \
$(top_builddir)/src/fvm/libfvm.a \
$(LDADD_CS_TESTS)

cs_map_test_SOURCES  = cs_map_test.cpp
cs_map_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_map_test_LDADD    = $(LDADD_CS_TESTS)

cs_matrix_test$(EXEEXT): $(top_srcdir)/tests/cs_matrix_test.c
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_matrix_test $(top_srcdir)/tests/cs_matrix_test.c

cs_mesh_quantities_test$(EXEEXT): $(top_srcdir)/tests/cs_mesh_quantities_test.c
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_mesh_quantities_test $(top_srcdir)/tests/cs_mesh_quantities_test.c

cs_moment_test_SOURCES  = cs_moment_test.c
cs_moment_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_moment_test_LDADD    = -lm

cs_random_test_SOURCES  = \
cs_random_test.cpp \
cs_random.cpp
cs_random_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_random_test_LDADD    = $(LDADD_CS_TESTS)

cs_rank_neighbors_test_SOURCES  = cs_rank_neighbors_test.cpp
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

cs_sizes_test_SOURCES  = cs_sizes_test.cpp
cs_sizes_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_sizes_test_LDADD    = $(LDADD_CS_TESTS)

cs_tree_test_SOURCES  = cs_tree_test.cpp
cs_tree_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_tree_test_LDADD    = $(LDADD_CS_TESTS)

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)

# Distribution
# Files not built through classical Automake rules need to be explicitely
# added to distribution.

EXTRA_DIST = \
unittests.py \
$(top_srcdir)/tests/graphics

# Clean

clean-local:
	-rm -f $(BUILT_SOURCES)
	-rm -rf .deps

.cu.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) -I$(top_srcdir)/src/base -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) -c -o $@ $<

## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
  LIBPLE_LA = $(top_builddir)/libple/src/libple.la
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
  LDFLAGS_PLE = $(PLE_LDFLAGS)
  LIBS_PLE = $(PLE_LIBS)
endif

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
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/gui \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS   = $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS  =

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
LDADD_BFT_TESTS = $(top_builddir)/src/bft/libbft.la

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

# MEI tests

check_PROGRAMS += \
mei_test

mei_test_SOURCES = mei_test_main.c mei_test_graph.c
mei_test_CPPFLAGS  = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_builddir)/src/mei \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)
mei_test_LDFLAGS =
mei_test_LDADD = $(top_builddir)/src/mei/libmei.la \
	$(top_builddir)/src/bft/libbft.la -lm

# Code_Saturne tests

check_PROGRAMS += \
cs_all_to_all_test \
cs_blas_test \
cs_check_cdo \
cs_check_quadrature \
cs_check_sdm \
cs_core_test \
cs_file_test \
cs_interface_test \
cs_map_test \
cs_matrix_test \
cs_moment_test \
cs_rank_neighbors_test \
fvm_selector_test \
fvm_selector_postfix_test \
cs_random_test \
cs_sizes_test \
cs_tree_test

LDFLAGS_CS_TESTS = $(CGNS_LDFLAGS) $(MED_LDFLAGS) $(HDF5_LDFLAGS) \
	$(PLE_LDFLAGS) $(MPI_LDFLAGS)
LDADD_CS_TESTS = \
	$(top_builddir)/src/base/libcscore.la \
	$(top_builddir)/src/bft/libbft.la \
	$(LIBPLE_LA) $(PLE_LIBS) $(MPI_LIBS) -lm

cs_all_to_all_test_SOURCES  = cs_all_to_all_test.c
cs_all_to_all_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_all_to_all_test_LDADD    = $(LDADD_CS_TESTS)

cs_blas_test_SOURCES  = \
cs_blas_test.c \
../src/alge/cs_blas.c
cs_blas_test_CPPFLAGS  = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
 $(BLAS_CPPFLAGS) $(MPI_CPPFLAGS)
cs_blas_test_CFLAGS  = $(AM_CFLAGS) $(CFLAGS_HOT)
cs_blas_test_LDFLAGS  = $(BLAS_LDFLAGS) $(LDFLAGS_CS_TESTS)
cs_blas_test_LDADD    = $(BLASRUNPATH) $(BLAS_LIBS) $(LDADD_CS_TESTS)
if HAVE_ESSL
cs_blas_test_LDADD    += $(FCLIBS)
endif

cs_check_cdo$(EXEEXT):
	PYTHONPATH=$(top_builddir)/bin:$(top_srcdir)/bin \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_check_cdo $(top_srcdir)/tests/cs_check_cdo.c

cs_check_quadrature$(EXEEXT):
	PYTHONPATH=$(top_builddir)/bin:$(top_srcdir)/bin \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_check_quadrature $(top_srcdir)/tests/cs_check_quadrature.c

cs_check_sdm$(EXEEXT):
	PYTHONPATH=$(top_builddir)/bin:$(top_srcdir)/bin \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	-o cs_check_sdm $(top_srcdir)/tests/cs_check_sdm.c

cs_core_test_SOURCES  = cs_core_test.c
cs_core_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_core_test_LDADD    = $(LDADD_CS_TESTS)

cs_file_test_SOURCES  = cs_file_test.c
cs_file_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_file_test_LDADD    = $(LDADD_CS_TESTS)

cs_interface_test_SOURCES  = cs_interface_test.c
cs_interface_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_interface_test_LDADD    = $(LDADD_CS_TESTS)

cs_map_test_SOURCES  = cs_map_test.c
cs_map_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_map_test_LDADD    = $(LDADD_CS_TESTS)

cs_matrix_test_SOURCES  = \
cs_matrix_test.c \
../src/base/cs_halo.c \
../src/base/cs_range_set.c \
../src/base/cs_sort.c \
../src/alge/cs_matrix.c \
../src/alge/cs_matrix_assembler.c
cs_matrix_test_CPPFLAGS  = \
-D_CS_UNIT_MATRIX_TEST \
$(AM_CPPFLAGS)
cs_matrix_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_matrix_test_LDADD    = $(LDADD_CS_TESTS)

cs_moment_test_SOURCES  = cs_moment_test.c
cs_moment_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_moment_test_LDADD    = $(LDADD_CS_TESTS)

cs_random_test_SOURCES  = \
cs_random_test.c \
../src/base/cs_random.c
cs_random_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_random_test_LDADD    = $(LDADD_CS_TESTS)

cs_rank_neighbors_test_SOURCES  = cs_rank_neighbors_test.c
cs_rank_neighbors_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_rank_neighbors_test_LDADD    = $(LDADD_CS_TESTS)

fvm_selector_test_SOURCES  = fvm_selector_test.c
fvm_selector_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
fvm_selector_test_LDADD    = $(LDADD_CS_TESTS)

fvm_selector_postfix_test_SOURCES  = fvm_selector_postfix_test.c
fvm_selector_postfix_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
fvm_selector_postfix_test_LDADD    = $(LDADD_CS_TESTS)

cs_sizes_test_SOURCES  = cs_sizes_test.c
cs_sizes_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_sizes_test_LDADD    = $(LDADD_CS_TESTS)

cs_tree_test_SOURCES  = cs_tree_test.c
cs_tree_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_tree_test_LDADD    = $(LDADD_CS_TESTS)

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)

# Distribution

EXTRA_DIST = \
unittests.py \
$(top_srcdir)/tests/graphics

dist-hook:
	rm -rf `find $(distdir)/tests/graphics -type d -name .svn`

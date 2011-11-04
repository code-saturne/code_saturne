## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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

check_PROGRAMS =

# BFT tests

check_PROGRAMS += \
bft_test \
bft_backtrace_test \
bft_error_test \
bft_mem_usage_test \
bft_mem_test \
bft_printf_test

CPPFLAGS_BFT_TESTS = -I$(top_srcdir)/src/bft
LDFLAGS_BFT_TESTS =
LDADD_BFT_TESTS = $(top_builddir)/src/bft/libbft.la

bft_test_SOURCES  = bft_test.c
bft_test_CPPFLAGS = $(CPPFLAGS_BFT_TESTS)
bft_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_test_LDADD    = $(LDADD_BFT_TESTS)

bft_backtrace_test_SOURCES  = bft_backtrace_test.c
bft_backtrace_test_CPPFLAGS = $(CPPFLAGS_BFT_TESTS)
bft_backtrace_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_backtrace_test_LDADD    = $(LDADD_BFT_TESTS)

bft_error_test_SOURCES  = bft_error_test.c
bft_error_test_CPPFLAGS = $(CPPFLAGS_BFT_TESTS)
bft_error_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_error_test_LDADD    = $(LDADD_BFT_TESTS)

bft_mem_usage_test_SOURCES  = bft_mem_usage_test.c
bft_mem_usage_test_CPPFLAGS = $(CPPFLAGS_BFT_TESTS)
bft_mem_usage_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_mem_usage_test_LDADD    = $(LDADD_BFT_TESTS)

bft_mem_test_SOURCES  = bft_mem_test.c
bft_mem_test_CPPFLAGS = $(CPPFLAGS_BFT_TESTS)
bft_mem_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_mem_test_LDADD    = $(LDADD_BFT_TESTS)

bft_printf_test_SOURCES  = bft_printf_test.c
bft_printf_test_CPPFLAGS = $(CPPFLAGS_BFT_TESTS)
bft_printf_test_LDFLAGS  = $(LDFLAGS_BFT_TESTS)
bft_printf_test_LDADD    = $(LDADD_BFT_TESTS)

# FVM tests

check_PROGRAMS += \
fvm_interface_test \
fvm_selector_test \
fvm_selector_postfix_test \
fvm_sizes_test

CPPFLAGS_FVM_TESTS = -I$(top_srcdir)/src/base -I$(top_srcdir)/src/bft \
       -I$(top_srcdir)/src/fvm $(MPI_CPPFLAGS)
LDFLAGS_FVM_TESTS = $(CGNS_LDFLAGS) $(MED_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
LDADD_FVM_TESTS = $(top_builddir)/src/fvm/libfvm.la \
        $(top_builddir)/src/bft/libbft.la \
        $(MPI_LIBS)

fvm_interface_test_SOURCES  = fvm_interface_test.c
fvm_interface_test_CPPFLAGS = $(CPPFLAGS_FVM_TESTS)
fvm_interface_test_LDFLAGS  = $(LDFLAGS_FVM_TESTS)
fvm_interface_test_LDADD    = $(LDADD_FVM_TESTS)

fvm_selector_test_SOURCES  = fvm_selector_test.c
fvm_selector_test_CPPFLAGS = $(CPPFLAGS_FVM_TESTS)
fvm_selector_test_LDFLAGS  = $(LDFLAGS_FVM_TESTS)
fvm_selector_test_LDADD    = $(LDADD_FVM_TESTS)

fvm_selector_postfix_test_SOURCES  = fvm_selector_postfix_test.c
fvm_selector_postfix_test_CPPFLAGS = $(CPPFLAGS_FVM_TESTS)
fvm_selector_postfix_test_LDFLAGS  = $(LDFLAGS_FVM_TESTS)
fvm_selector_postfix_test_LDADD    = $(LDADD_FVM_TESTS)

fvm_sizes_test_SOURCES  = fvm_sizes_test.c
fvm_sizes_test_CPPFLAGS = $(CPPFLAGS_FVM_TESTS)
fvm_sizes_test_LDFLAGS  = $(LDFLAGS_FVM_TESTS)
fvm_sizes_test_LDADD    = $(LDADD_FVM_TESTS)

# MEI tests

check_PROGRAMS += \
mei_test

mei_test_SOURCES = mei_test_main.c mei_test_graph.c
mei_test_CPPFLAGS = -I$(top_srcdir)/src/bft \
        -I$(top_srcdir)/src/mei -I$(top_builddir)/src/mei
mei_test_LDFLAGS =
mei_test_LDADD = $(top_builddir)/src/mei/libmei.la \
        $(top_builddir)/src/bft/libbft.la -lm

# Code_Saturne tests

check_PROGRAMS += \
cs_file_test

CPPFLAGS_CS_TESTS = -I$(top_srcdir)/src/base -I$(top_srcdir)/src/bft \
       -I$(top_srcdir)/src/fvm $(MPI_CPPFLAGS)
LDFLAGS_CS_TESTS = $(CGNS_LDFLAGS) $(MED_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
LDADD_CS_TESTS =  $(top_builddir)/src/base/libcscore.la \
        $(top_builddir)/src/fvm/libfvm.la \
        $(top_builddir)/src/bft/libbft.la \
        $(MPI_LIBS)

cs_file_test_SOURCES  = cs_file_test.c
cs_file_test_CPPFLAGS = $(CPPFLAGS_CS_TESTS)
cs_file_test_LDFLAGS  = $(LDFLAGS_CS_TESTS)
cs_file_test_LDADD    = $(LDADD_CS_TESTS)

# Graphical User Interface unit tests

EXTRA_DIST = \
unittests.py

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)


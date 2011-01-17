# Process this file with automake to produce Makefile.in

#------------------------------------------------------------------------------
#   This file is part of the "Mathematical Expression Interpreter" library.
#
#   Copyright (C) 2008-2010  EDF
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#------------------------------------------------------------------------------

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
fvm_file_test \
fvm_interface_test \
fvm_selector_test \
fvm_selector_postfix_test \
fvm_sizes_test

CPPFLAGS_FVM_TESTS = -I$(top_srcdir)/src/bft -I$(top_srcdir)/src/fvm \
        $(MPI_CPPFLAGS)
LDFLAGS_FVM_TESTS = $(CGNS_LDFLAGS) $(MED_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
LDADD_FVM_TESTS = $(top_builddir)/src/fvm/libfvm.la \
        $(top_builddir)/src/bft/libbft.la \
        $(CGNS_LIBS) $(MED_LIBS) $(HDF5_LIBS) $(MPI_LIBS)

fvm_file_test_SOURCES  = fvm_file_test.c
fvm_file_test_CPPFLAGS = $(CPPFLAGS_FVM_TESTS)
fvm_file_test_LDFLAGS  = $(LDFLAGS_FVM_TESTS)
fvm_file_test_LDADD    = $(LDADD_FVM_TESTS)

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
mei_test_CPPFLAGS = -I$(top_srcdir)/src/bft -I$(top_srcdir)/src/mei
mei_test_LDFLAGS =
mei_test_LDADD = $(top_builddir)/src/mei/libmei.la \
        $(top_builddir)/src/bft/libbft.la -lm

# CFD Proxy tests

check_PROGRAMS += \
cfd_proxy_unit_test

cfd_proxy_unit_test_SOURCES = cfd_proxy_unit_test.c
cfd_proxy_unit_test_CPPFLAGS = -I$(top_srcdir)/salome/cfd_proxy \
        $(MPI_CPPFLAGS) $(SALOME_KERNEL_CPPFLAGS)
cfd_proxy_unit_test_LDFLAGS = $(MPI_LDFLAGS) $(SALOME_KERNEL_LDFLAGS)
cfd_proxy_unit_test_LDADD = $(top_builddir)/salome/cfd_proxy/libcfdproxy.la \
        $(MPI_LIBS) -lm

# Uncomment for tests execution at "make check"
#TESTS=$(check_PROGRAMS)


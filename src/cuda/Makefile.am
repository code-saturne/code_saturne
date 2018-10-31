## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) IBM Corp. 2017, 2018
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

# New suffixes and targets

SUFFIXES = .cu .o

# Main part

AM_CUDAFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/cdo \
$(CXXFLAGS) $(CUDAFLAGS)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_cuda.h \
$(top_srcdir)/src/alge/cs_matrix.h

# Library source files

noinst_LTLIBRARIES = libcscuda.la

libcscuda_la_CUDAINCOMPATIBLEFLAGS = \
-W \
-Wall \
-Wshadow \
-Wpointer-arith \
-Wcast-qual \
-Wcast-align \
-Wwrite-strings \
-Wunused \
-Wfloat-equal \
-Wmisleading-indentation \
-Wduplicated-cond \
-fdiagnostics-color=auto

libcscuda_la_CUDAFLAGS = $(filter-out $(libcscuda_la_CUDAINCOMPATIBLEFLAGS),$(AM_CUDAFLAGS))
libcscuda_la_SOURCES = \
cs_cuda.cu

libcscuda_la_CUDAOBJECTS = $(addprefix libcscuda_la-,$(libcscuda_la_SOURCES:.cu=.lo))

libcscuda_la_LDFLAGS = -no-undefined

# We compile the cuda files with a special rule, so we pass them as dependencies to the 
# libtool module.
libcscuda_la_LIBADD = $(libcscuda_la_CUDAOBJECTS)
libcscuda_la_DEPENDENCIES = $(libcscuda_la_CUDAOBJECTS)

# If we con't support CUDA - compile the CUDA file as C++ code.
if HAVE_CUDA_OFFLOAD
libcscuda_la-%.lo: %.cu $(pkginclude_HEADERS)
	$(AM_V_CXX)$(LIBTOOL) $(AM_V_lt) \
	--tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile \
	$(CUDACC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libcscuda_la_CUDAFLAGS) \
	-c -o $@ \
	`test -f "$<" || echo "$(srcdir)/"`$<
else
libcscuda_la-%.lo: %.cu $(pkginclude_HEADERS)
	$(AM_V_CXX)$(LIBTOOL) $(AM_V_lt) \
	--tag=CXX $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=compile \
	$(CXX) -x c++ $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(libcscuda_la_CUDAFLAGS) \
	-c -o $@ \
	`test -f "$<" || echo "$(srcdir)/"`$<
endif
clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

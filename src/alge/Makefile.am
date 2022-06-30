## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Libtool modified for CUDA

LIBTOOL_CUDA = $(SHELL) $(top_builddir)/libtool_cuda

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/cdo \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS) $(BLAS_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_LDFLAGS =

# Conditionally compiled extensions

libcsalge_petsc_la_CPPFLAGS =\
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
$(PETSC_CPPFLAGS) $(MPI_CPPFLAGS)

if HAVE_MUMPS
  libcsalge_mumps_la_CPPFLAGS =\
  -I$(top_srcdir)/src/base \
  -I$(top_srcdir)/src/bft \
  -I$(top_srcdir)/src/fvm \
  -I$(top_srcdir)/src/mesh \
  $(MUMPS_CPPFLAGS) $(MPI_CPPFLAGS)
endif

libcsalge_amgx_la_CPPFLAGS =\
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
$(AMGX_CPPFLAGS) $(MPI_CPPFLAGS)

libcsalge_hypre_la_CPPFLAGS =\
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
$(HYPRE_CPPFLAGS) $(MPI_CPPFLAGS)

libcsalge_extension_la_CPPFLAGS =\
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
$(PETSC_CPPFLAGS)  $(MUMPS_CPPFLAGS) $(MPI_CPPFLAGS)

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_alge_headers.h \
cs_bad_cells_regularisation.h \
cs_balance.h \
cs_balance_by_zone.h \
cs_benchmark.h \
cs_benchmark_matrix.h \
cs_blas.h \
cs_cell_to_vertex.h \
cs_convection_diffusion.h \
cs_divergence.h \
cs_face_viscosity.h \
cs_gradient.h \
cs_gradient_priv.h \
cs_grid.h \
cs_matrix.h \
cs_matrix_assembler.h \
cs_matrix_assembler_priv.h \
cs_matrix_building.h \
cs_matrix_default.h \
cs_matrix_priv.h \
cs_matrix_spmv.h \
cs_matrix_tuning.h \
cs_matrix_util.h \
cs_multigrid.h \
cs_multigrid_smoother.h \
cs_param_sles.h \
cs_sles.h \
cs_sles_default.h \
cs_sles_it.h \
cs_sles_it_priv.h \
cs_sles_pc.h \
cs_sles_pc_priv.h

if HAVE_CUDA
pkginclude_HEADERS += \
cs_benchmark_cuda.h \
cs_blas_cuda.h \
cs_matrix_spmv_cuda.h \
cs_sles_it_cuda.h \
cs_sles_pc_cuda.h
endif

if HAVE_PETSC
pkginclude_HEADERS += \
cs_matrix_petsc.h \
cs_matrix_petsc_priv.h \
cs_sles_petsc.h
endif

if HAVE_MUMPS
pkginclude_HEADERS += cs_sles_mumps.h
endif

if HAVE_HYPRE
pkginclude_HEADERS += \
cs_matrix_hypre.h \
cs_matrix_hypre_priv.h \
cs_sles_hypre.h
endif

if HAVE_AMGX
pkginclude_HEADERS += cs_sles_amgx.h
endif

# Library source files

noinst_LTLIBRARIES = libcsalge.la

libcsalge_la_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_HOT)
libcsalge_la_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_HOT)
libcsalge_la_SOURCES = \
cs_bad_cells_regularisation.c \
cs_balance.c \
cs_balance_by_zone.c \
cs_benchmark.c \
cs_benchmark_matrix.c \
cs_blas.c \
cs_bw_time_diff.c \
cs_cell_to_vertex.c \
cs_convection_diffusion.c \
cs_divergence.c \
cs_face_viscosity.c \
cs_gradient.c \
cs_grid.c \
cs_matrix.c \
cs_matrix_assembler.c \
cs_matrix_building.c \
cs_matrix_default.c \
cs_matrix_spmv.c \
cs_matrix_tuning.c \
cs_matrix_util.c \
cs_multigrid.c \
cs_multigrid_smoother.c \
cs_sles.c \
cs_sles_default.c \
cs_sles_it.c \
cs_sles_it_priv.c \
cs_sles_pc.c
libcsalge_la_LDFLAGS = -no-undefined
libcsalge_la_LIBADD =

if HAVE_CUDA
libcsalge_la_SOURCES += \
cs_benchmark_cuda.cu \
cs_blas_cuda.cu \
cs_gradient_cuda.cu \
cs_matrix_spmv_cuda.cu \
cs_sles_it_cuda.cu \
cs_sles_pc_cuda.cu
endif

noinst_LTLIBRARIES += libcsalge_extension.la
libcsalge_la_LIBADD += libcsalge_extension.la
libcsalge_extension_la_CPPFLAGS = \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mesh \
$(PETSC_CPPFLAGS) $(MUMPS_CPPFLAGS) $(AMGX_CPPFLAGS) $(HYPRE_CPPFLAGS) $(MPI_CPPFLAGS)
libcsalge_extension_la_SOURCES = \
cs_param_sles.c

if HAVE_PETSC
noinst_LTLIBRARIES += libcsalge_petsc.la
libcsalge_la_LIBADD += libcsalge_petsc.la
libcsalge_petsc_la_SOURCES = \
cs_matrix_petsc.c \
cs_sles_petsc.c
endif

if HAVE_MUMPS
noinst_LTLIBRARIES += libcsalge_mumps.la
libcsalge_la_LIBADD += libcsalge_mumps.la
libcsalge_mumps_la_SOURCES = \
cs_sles_mumps.c
libcsalge_mumps_la_LIBADD = \
$(MUMPS_LDFLAGS) $(MUMPSRUNPATH) $(MUMPS_LIBS)
libcsalge_mumps_la_LDFLAGS = -no-undefined
endif

if HAVE_HYPRE
noinst_LTLIBRARIES += libcsalge_hypre.la
libcsalge_la_LIBADD += libcsalge_hypre.la
libcsalge_hypre_la_SOURCES = \
cs_matrix_hypre.c \
cs_sles_hypre.c
endif

if HAVE_AMGX
noinst_LTLIBRARIES += libcsalge_amgx.la
libcsalge_la_LIBADD += libcsalge_amgx.la
libcsalge_amgx_la_SOURCES = cs_sles_amgx.c
endif

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

.cu.lo:
	$(LIBTOOL_CUDA) --tag=CC --mode=compile $(NVCC) $(AM_CPPFLAGS) -I$(top_srcdir)/src/alge -I../.. $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) -c -o $@ $<

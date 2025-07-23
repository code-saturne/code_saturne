## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2025 EDF S.A.
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

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(CPPFLAGS) $(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS) $(BLAS_CPPFLAGS)

AM_CXXFLAGS = $(CXXFLAGS_SYCL) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Conditionally compiled extensions

libcsalge_extension_a_CPPFLAGS =\
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(CPPFLAGS_PLE) \
$(PETSC_CPPFLAGS) $(MUMPS_CPPFLAGS) $(HYPRE_CPPFLAGS) \
$(MPI_CPPFLAGS)

libcsalge_petsc_a_CPPFLAGS =\
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(PETSC_CPPFLAGS) $(MPI_CPPFLAGS)

if HAVE_MUMPS
  libcsalge_mumps_a_CPPFLAGS =\
  -I$(top_srcdir)/src -I$(top_builddir)/src \
  $(MUMPS_CPPFLAGS) $(MPI_CPPFLAGS)
endif

libcsalge_amgx_a_CPPFLAGS =\
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(AMGX_CPPFLAGS) $(MPI_CPPFLAGS)

libcsalge_cudss_a_CPPFLAGS =\
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(CUDSS_CPPFLAGS) $(MPI_CPPFLAGS)

libcsalge_hypre_a_CPPFLAGS =\
-I$(top_srcdir)/src -I$(top_builddir)/src \
$(HYPRE_CPPFLAGS) $(MPI_CPPFLAGS)

noinst_HEADERS = \
cs_convection_diffusion_priv.h

# Public header files (to be installed)

csalge_includedir=$(pkgincludedir)/alge
csalge_include_HEADERS = \
cs_alge_headers.h \
cs_bad_cells_regularisation.h \
cs_balance.h \
cs_balance_by_zone.h \
cs_benchmark.h \
cs_benchmark_matrix.h \
cs_blas.h \
cs_bw_time_diff.h \
cs_cell_to_vertex.h \
cs_convection_diffusion.h \
cs_divergence.h \
cs_face_viscosity.h \
cs_gradient.h \
cs_gradient_boundary.h \
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
cs_param_amg.h \
cs_param_hpddm.h \
cs_param_mumps.h \
cs_param_saddle.h \
cs_param_sles.h \
cs_param_sles_setup.h \
cs_saddle_solver.h \
cs_saddle_solver_setup.h \
cs_sles.h \
cs_sles_default.h \
cs_sles_it.h \
cs_sles_it_priv.h \
cs_sles_pc.h \
cs_sles_pc_priv.h \
cs_vertex_to_cell.h

if HAVE_CUDA
csalge_include_HEADERS += \
cs_benchmark_cuda.h \
cs_blas_cuda.h \
cs_matrix_spmv_cuda.h \
cs_sles_it_cuda.h \
cs_sles_pc_cuda.h
endif

if HAVE_PETSC
csalge_include_HEADERS += \
cs_matrix_petsc.h \
cs_matrix_petsc_priv.h \
cs_sles_petsc.h
endif

if HAVE_MUMPS
csalge_include_HEADERS += cs_sles_mumps.h
endif

if HAVE_HYPRE
csalge_include_HEADERS += \
cs_matrix_hypre.h \
cs_matrix_hypre_priv.h \
cs_sles_hypre.h
endif

if HAVE_AMGX
csalge_include_HEADERS += cs_sles_amgx.h
endif

if HAVE_CUDSS
csalge_include_HEADERS += cs_sles_cudss.h
endif

# Library source files

noinst_LIBRARIES = libcsalge.a

libcsalge_a_SOURCES = \
cs_bad_cells_regularisation.cpp \
cs_balance.cpp \
cs_balance_by_zone.cpp \
cs_benchmark.cpp \
cs_benchmark_matrix.cpp \
cs_blas.cpp \
cs_bw_time_diff.cpp \
cs_cell_to_vertex.cpp \
cs_convection_diffusion.cpp \
cs_convection_diffusion_steady.cpp \
cs_divergence.cpp \
cs_face_viscosity.cpp \
cs_gradient.cpp \
cs_gradient_boundary.cpp \
cs_grid.cpp \
cs_matrix.cpp \
cs_matrix_assembler.cpp \
cs_matrix_building.cpp \
cs_matrix_default.cpp \
cs_matrix_spmv.cpp \
cs_matrix_tuning.cpp \
cs_matrix_util.cpp \
cs_multigrid.cpp \
cs_multigrid_smoother.cpp \
cs_param_amg.cpp \
cs_param_hpddm.cpp \
cs_param_mumps.cpp \
cs_param_saddle.cpp \
cs_saddle_solver.cpp \
cs_sles.cpp \
cs_sles_default.cpp \
cs_sles_it.cpp \
cs_sles_it_priv.cpp \
cs_sles_pc.cpp \
cs_vertex_to_cell.cpp
libcsalge_a_LIBADD =

if HAVE_CUDA
libcsalge_a_SOURCES += \
cs_benchmark_cuda.cu \
cs_blas_cuda.cu \
cs_gradient_cuda.cu \
cs_matrix_spmv_cuda.cu \
cs_sles_it_cuda.cu \
cs_sles_pc_cuda.cu
endif

noinst_LIBRARIES += libcsalge_extension.a
libcsalge_extension_a_SOURCES = \
cs_param_sles.cpp \
cs_param_sles_setup.cpp \
cs_saddle_solver_setup.cpp
libcsalge_a_LIBADD += \
libcsalge_extension_a-cs_param_sles.$(OBJEXT) \
libcsalge_extension_a-cs_param_sles_setup.$(OBJEXT) \
libcsalge_extension_a-cs_saddle_solver_setup.$(OBJEXT)

if HAVE_PETSC
noinst_LIBRARIES += libcsalge_petsc.a
libcsalge_petsc_a_SOURCES = \
cs_matrix_petsc.cpp \
cs_sles_petsc.cpp
libcsalge_a_LIBADD += \
libcsalge_petsc_a-cs_matrix_petsc.$(OBJEXT) \
libcsalge_petsc_a-cs_sles_petsc.$(OBJEXT)

if HAVE_CUDA_CPP
cs_matrix_petsc = $(top_srcdir)/src/alge/cs_matrix_petsc.cpp
cs_sles_petsc = $(top_srcdir)/src/alge/cs_sles_petsc.cpp
libcsalge_petsc_a-cs_matrix_petsc.$(OBJEXT): $(cs_matrix_petsc)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(PETSC_CPPFLAGS) $(NVCCFLAGS) $(NVCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_matrix_petsc) -o libcsalge_petsc_a-cs_matrix_petsc.$(OBJEXT)
libcsalge_petsc_a-cs_sles_petsc.$(OBJEXT): $(cs_sles_petsc)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(PETSC_CPPFLAGS) $(NVCCFLAGS) $(NVCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_sles_petsc) -o libcsalge_petsc_a-cs_sles_petsc.$(OBJEXT)
endif
endif

if HAVE_MUMPS
noinst_LIBRARIES += libcsalge_mumps.a
libcsalge_mumps_a_SOURCES = \
cs_sles_mumps.cpp
libcsalge_a_LIBADD += libcsalge_mumps_a-cs_sles_mumps.$(OBJEXT)

if HAVE_CUDA_CPP
cs_sles_mumps = $(top_srcdir)/src/alge/cs_sles_mumps.cpp
libcsalge_mumps_a-cs_sles_mumps.$(OBJEXT): $(cs_sles_mumps)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(MUMPS_CPPFLAGS) $(NVCCFLAGS) $(NVCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_sles_mumps) -o libcsalge_mumps_a-cs_sles_mumps.$(OBJEXT)
endif
endif

if HAVE_HYPRE
noinst_LIBRARIES += libcsalge_hypre.a
libcsalge_hypre_a_SOURCES = \
cs_matrix_hypre.cpp \
cs_sles_hypre.cpp
libcsalge_a_LIBADD += \
libcsalge_hypre_a-cs_matrix_hypre.$(OBJEXT) \
libcsalge_hypre_a-cs_sles_hypre.$(OBJEXT)

if HAVE_CUDA_CPP
cs_matrix_hypre = $(top_srcdir)/src/alge/cs_matrix_hypre.cpp
cs_sles_hypre = $(top_srcdir)/src/alge/cs_sles_hypre.cpp
libcsalge_hypre_a-cs_matrix_hypre.$(OBJEXT): $(cs_matrix_hypre)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(HYPRE_CPPFLAGS) $(NVCCFLAGS) $(NVCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_matrix_hypre) -o libcsalge_hypre_a-cs_matrix_hypre.$(OBJEXT)
libcsalge_hypre_a-cs_sles_hypre.$(OBJEXT): $(cs_sles_hypre)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(HYPRE_CPPFLAGS) $(NVCCFLAGS) $(NVCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_sles_hypre) -o libcsalge_hypre_a-cs_sles_hypre.$(OBJEXT)
endif
endif

if HAVE_AMGX
noinst_LIBRARIES += libcsalge_amgx.a
libcsalge_amgx_a_SOURCES = cs_sles_amgx.cpp
libcsalge_a_LIBADD += libcsalge_amgx_a-cs_sles_amgx.$(OBJEXT)

if HAVE_CUDA_CPP
cs_sles_amgx = $(top_srcdir)/src/alge/cs_sles_amgx.cpp
libcsalge_amgx_a-cs_sles_amgx.$(OBJEXT): $(cs_sles_amgx)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(AMGX_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_sles_amgx) -o libcsalge_amgx_a-cs_sles_amgx.$(OBJEXT)
endif
endif

if HAVE_CUDSS
noinst_LIBRARIES += libcsalge_cudss.a
libcsalge_cudss_a_SOURCES = cs_sles_cudss.cpp
libcsalge_a_LIBADD += libcsalge_cudss_a-cs_sles_cudss.$(OBJEXT)

if HAVE_CUDA_CPP
cs_sles_cudss = $(top_srcdir)/src/alge/cs_sles_cudss.cpp
libcsalge_cudss_a-cs_sles_cudss.$(OBJEXT): $(cs_sles_cudss)
	$(NVCC) -c $(DEFS) $(AM_CPPFLAGS) $(CUDSS_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c  $(cs_sles_cudss) -o libcsalge_cudss_a-cs_sles_cudss.$(OBJEXT)
endif
endif

clean-local:
	-rm -rf .deps

if HAVE_CUDA_CPP
.cpp.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) $(NVCCFLAGS_CPP) -c -o $@ $<
endif

.cu.$(OBJEXT):
	$(NVCC) $(DEFS) $(AM_CPPFLAGS) $(NVCCFLAGS) $(NVCCFLAGS_DBG) $(NVCCFLAGS_OPT) -c -o $@ $<

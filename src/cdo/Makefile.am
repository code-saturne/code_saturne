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
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft  \
-I$(top_srcdir)/src/fvm  \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/gwf  \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) $(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)
AM_CXXFLAGS = $(CXXFLAGS_STD) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT) $(CXXFLAGS_SHARED)

AM_LDFLAGS =

# Conditionally compiled files relying on PETSc/MUMPS libraries

libcscdo_extension_a_CPPFLAGS =\
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft  \
-I$(top_srcdir)/src/fvm  \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/gwf  \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(PETSC_CPPFLAGS) $(MUMPS_CPPFLAGS) \
$(MPI_CPPFLAGS)

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_advection_field.h \
cs_basis_func.h \
cs_cdo_advection.h \
cs_cdo_assembly.h \
cs_cdo_bc.h \
cs_cdo_blas.h \
cs_cdo_connect.h \
cs_cdo_diffusion.h \
cs_cdo_field_interpolation.h \
cs_cdo_local.h \
cs_cdo_main.h \
cs_cdo_quantities.h \
cs_cdo_solve.h \
cs_cdo_system.h \
cs_cdo_toolbox.h \
cs_cdo_turbulence.h \
cs_cdocb_scaleq.h \
cs_cdocb_scaleq_sles.h \
cs_cdocb_priv.h \
cs_cdoeb_priv.h \
cs_cdoeb_vecteq.h \
cs_cdofb_ac.h \
cs_cdofb_monolithic.h \
cs_cdofb_monolithic_priv.h \
cs_cdofb_monolithic_sles.h \
cs_cdofb_navsto.h \
cs_cdofb_predco.h \
cs_cdofb_priv.h \
cs_cdofb_scaleq.h \
cs_cdofb_vecteq.h \
cs_cdovb_priv.h \
cs_cdovb_scaleq.h \
cs_cdovb_scalsys.h \
cs_cdovb_vecteq.h \
cs_cdovcb_scaleq.h \
cs_dbg.h \
cs_domain.h \
cs_domain_op.h \
cs_domain_setup.h \
cs_enforcement.h \
cs_equation.h \
cs_equation_bc.h \
cs_equation_builder.h \
cs_equation_param.h \
cs_equation_priv.h \
cs_equation_system.h \
cs_equation_system_param.h \
cs_equation_system_sles.h \
cs_evaluate.h \
cs_flag.h \
cs_hho_builder.h \
cs_hho_scaleq.h \
cs_hho_stokes.h \
cs_hho_vecteq.h \
cs_hodge.h \
cs_iter_algo.h \
cs_macfb_advection.h \
cs_macfb_builder.h \
cs_macfb_monolithic.h \
cs_macfb_monolithic_priv.h \
cs_macfb_monolithic_sles.h \
cs_macfb_navsto.h \
cs_macfb_priv.h \
cs_macfb_vecteq.h \
cs_maxwell.h \
cs_mesh_deform.h \
cs_navsto_coupling.h \
cs_navsto_param.h \
cs_navsto_system.h \
cs_param_cdo.h \
cs_property.h \
cs_quadrature.h \
cs_reco.h \
cs_reco_cw.h \
cs_saddle_system.h \
cs_scheme_geometry.h \
cs_sdm.h \
cs_solid_selection.h \
cs_solidification.h \
cs_source_term.h \
cs_static_condensation.h \
cs_thermal_system.h \
cs_walldistance.h \
cs_xdef.h \
cs_xdef_cw_eval.h \
cs_xdef_eval.h \
cs_cdo_headers.h

# Library source files

noinst_LIBRARIES = libcscdo.a
libcscdo_a_SOURCES = \
cs_advection_field.cxx \
cs_basis_func.cxx \
cs_cdo_advection.cxx \
cs_cdo_assembly.cxx \
cs_cdo_bc.cxx \
cs_cdo_blas.cxx \
cs_cdo_connect.cxx \
cs_cdo_diffusion.cxx \
cs_cdo_field_interpolation.cxx \
cs_cdo_local.cxx \
cs_cdo_main.cxx \
cs_cdo_quantities.cxx \
cs_cdo_solve.cxx \
cs_cdo_system.cxx \
cs_cdo_toolbox.cxx \
cs_cdo_turbulence.cxx \
cs_cdocb_priv.cxx \
cs_cdocb_scaleq.cxx \
cs_cdocb_scaleq_sles.cxx \
cs_cdoeb_vecteq.cxx \
cs_cdofb_ac.cxx \
cs_cdofb_monolithic.cxx \
cs_cdofb_monolithic_sles.cxx \
cs_cdofb_navsto.cxx \
cs_cdofb_predco.cxx \
cs_cdofb_priv.cxx \
cs_cdofb_scaleq.cxx \
cs_cdofb_vecteq.cxx \
cs_cdovb_scaleq.cxx \
cs_cdovb_scalsys.cxx \
cs_cdovb_vecteq.cxx \
cs_cdovcb_scaleq.cxx \
cs_dbg.cxx \
cs_domain.cxx \
cs_domain_op.cxx \
cs_domain_setup.cxx \
cs_enforcement.cxx \
cs_equation.cxx \
cs_equation_bc.cxx \
cs_equation_builder.cxx \
cs_equation_param.cxx \
cs_equation_system.cxx \
cs_equation_system_param.cxx \
cs_evaluate.cxx \
cs_flag.cxx \
cs_hho_builder.cxx \
cs_hho_scaleq.cxx \
cs_hho_stokes.cxx \
cs_hho_vecteq.cxx \
cs_hodge.cxx \
cs_iter_algo.cxx \
cs_macfb_advection.cxx \
cs_macfb_builder.cxx \
cs_macfb_monolithic.cxx \
cs_macfb_monolithic_sles.cxx \
cs_macfb_navsto.cxx \
cs_macfb_priv.cxx \
cs_macfb_vecteq.cxx \
cs_maxwell.cxx \
cs_mesh_deform.cxx \
cs_navsto_coupling.cxx \
cs_navsto_param.cxx \
cs_navsto_system.cxx \
cs_param_cdo.cxx \
cs_property.cxx \
cs_quadrature.cxx \
cs_reco.cxx \
cs_reco_cw.cxx \
cs_saddle_system.cxx \
cs_scheme_geometry.cxx \
cs_sdm.cxx \
cs_solid_selection.cxx \
cs_solidification.cxx \
cs_source_term.cxx \
cs_static_condensation.cxx \
cs_thermal_system.cxx \
cs_walldistance.cxx \
cs_xdef.cxx \
cs_xdef_cw_eval.cxx \
cs_xdef_eval.cxx

libcscdo_a_LIBADD =

noinst_LIBRARIES += libcscdo_extension.a
libcscdo_extension_a_SOURCES = \
cs_equation_system_sles.cxx
libcscdo_a_LIBADD += \
libcscdo_extension_a-cs_equation_system_sles.$(OBJEXT)

clean-local:
	-rm -rf .deps

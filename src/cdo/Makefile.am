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
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
$(CPPFLAGS_PLE) \
$(PETSC_CPPFLAGS) $(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_advection_field.h \
cs_basis_func.h \
cs_cdo_advection.h \
cs_cdo_bc.h \
cs_cdo_connect.h \
cs_cdo_diffusion.h \
cs_cdo_local.h \
cs_cdo_main.h \
cs_cdo_quantities.h \
cs_cdo_time.h \
cs_cdofb_scaleq.h \
cs_cdofb_stokes.h \
cs_cdofb_vecteq.h \
cs_cdovb_scaleq.h \
cs_cdovcb_scaleq.h \
cs_dbg.h \
cs_domain.h \
cs_domain_post.h \
cs_domain_setup.h \
cs_equation.h \
cs_equation_bc.h \
cs_equation_common.h \
cs_equation_param.h \
cs_evaluate.h \
cs_flag.h \
cs_gwf.h \
cs_gwf_soil.h \
cs_gwf_tracer.h \
cs_hho_builder.h \
cs_hho_scaleq.h \
cs_hho_stokes.h \
cs_hodge.h \
cs_mesh_deform.h \
cs_navsto_param.h \
cs_navsto_system.h \
cs_param.h \
cs_param_cdo.h \
cs_property.h \
cs_quadrature.h \
cs_reco.h \
cs_scheme_geometry.h \
cs_sdm.h \
cs_sla.h \
cs_source_term.h \
cs_walldistance.h \
cs_xdef.h \
cs_xdef_eval.h

# Library source files

noinst_LTLIBRARIES = libcscdo.la
libcscdo_la_SOURCES = \
cs_advection_field.c \
cs_basis_func.c \
cs_cdo_advection.c \
cs_cdo_bc.c \
cs_cdo_connect.c \
cs_cdo_diffusion.c \
cs_cdo_local.c \
cs_cdo_main.c \
cs_cdo_quantities.c \
cs_cdo_time.c \
cs_cdofb_scaleq.c \
cs_cdofb_stokes.c \
cs_cdofb_vecteq.c \
cs_cdovb_scaleq.c \
cs_cdovcb_scaleq.c \
cs_dbg.c \
cs_domain.c \
cs_domain_setup.c \
cs_domain_post.c \
cs_equation.c \
cs_equation_bc.c \
cs_equation_common.c \
cs_equation_param.c \
cs_evaluate.c \
cs_flag.c \
cs_gwf.c \
cs_gwf_soil.c \
cs_gwf_tracer.c \
cs_hho_builder.c \
cs_hho_scaleq.c \
cs_hho_stokes.c \
cs_hodge.c \
cs_mesh_deform.c \
cs_navsto_param.c \
cs_navsto_system.c \
cs_param.c \
cs_param_cdo.c \
cs_property.c \
cs_quadrature.c \
cs_reco.c \
cs_scheme_geometry.c \
cs_sdm.c \
cs_sla.c \
cs_source_term.c \
cs_walldistance.c \
cs_xdef.c \
cs_xdef_eval.c

libcscdo_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *.mod

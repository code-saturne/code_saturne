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
  CPPFLAGS_PLE = -I"$(top_srcdir)/libple/src" -I"$(top_builddir)/libple"
  LIBPLE_LA = $(top_builddir)/libple/src/libple.la
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
  LDFLAGS_PLE = $(PLE_LDFLAGS)
  LIBS_PLE = $(PLE_LIBS)
endif

# Relocatable installation

if HAVE_RELOCATABLE
  ORIGINRUNPATH = -R '\$$ORIGIN/../../lib'
  PLUGINRUNPATH_R = -R'\$$ORIGIN/../../lib/'$(PACKAGE)
else
  PLUGINRUNPATH_N = -R$(pkglibdir)
endif

if HAVE_PLUGINS
  PLUGINRUNPATH = $(PLUGINRUNPATH_R)$(PLUGINRUNPATH_N)
endif

# BLAS

if HAVE_MKL
  CPPFLAGS_BLAS =  $(BLAS_CPPFLAGS)
  LDADD_BLAS = $(BLAS_LDFLAGS) $(BLASRUNPATH) $(BLAS_LIBS)
endif

# Catalyst

if HAVE_CATALYST
if HAVE_PLUGIN_CATALYST
else
  LDADD_CATALYST = $(CATALYST_LDFLAGS) $(CATALYSTRUNPATH) $(CATALYST_LIBS)
endif
endif

# MEDCoupling

if HAVE_MEDCOUPLING
if HAVE_PLUGIN_MEDCOUPLING
else
  LDADD_MEDCOUPLING = $(MEDCOUPLING_LDFLAGS) $(MEDCOUPLINGRUNPATH) $(MEDCOUPLING_LIBS)
endif
endif

# Melissa

if HAVE_MELISSA
if HAVE_PLUGIN_MELISSA
else
  LDADD_MELISSA = $(MELISSA_LDFLAGS) $(MELISSARUNPATH) $(MELISSA_LIBS)
endif
endif

# CoolProp

if HAVE_COOLPROP
if HAVE_PLUGINS
else
  LDADD_COOLPROP = $(COOLPROP_LDFLAGS) $(COOLPROPRUNPATH) $(COOLPROP_LIBS)
endif
endif

# Main part

if HAVE_BACKEND

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/elec \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/darc \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS) $(LIBXML2_CPPFLAGS) $(CPPFLAGS_BLAS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

# Global LDFLAGS must be used so as to handle the possibility of a CXX link.
# Per-executable flags will override this global definition.
AM_LDFLAGS = $(ORIGINRUNPATH) $(PLUGINRUNPATH)

# Library source files

lib_LTLIBRARIES = libsaturne.la libcs_solver.la

libsaturne_la_SOURCES =
libsaturne_la_LIBADD = \
$(top_builddir)/src/fvm/libfvm_filters.la \
$(top_builddir)/src/base/libcsbase.la \
$(top_builddir)/src/base/libcscore.la \
$(top_builddir)/src/base/libcsrenumber.la \
$(top_builddir)/src/pprt/libcspprt.la \
$(top_builddir)/src/alge/libcsalge.la \
$(top_builddir)/src/mesh/libcsmesh.la \
$(top_builddir)/src/mesh/libcspartition.la \
$(top_builddir)/src/turb/libcsturb.la \
$(top_builddir)/src/atmo/libcsatmo.la \
$(top_builddir)/src/cdo/libcscdo.la \
$(top_builddir)/src/cfbl/libcscfbl.la \
$(top_builddir)/src/cogz/libcscogz.la \
$(top_builddir)/src/comb/libcscomb.la \
$(top_builddir)/src/ctwr/libcsctwr.la \
$(top_builddir)/src/darc/libcsdarc.la \
$(top_builddir)/src/elec/libcselec.la \
$(top_builddir)/src/gui/libcsgui.la \
$(top_builddir)/src/lagr/libcslagr.la \
$(top_builddir)/src/rayt/libcsrayt.la \
$(top_builddir)/src/user/libcsuser.la \
$(top_builddir)/src/libcsmodules.la \
$(top_builddir)/src/mei/libmei.la \
$(top_builddir)/src/bft/libbft.la \
$(LIBPLE_LA) \
$(LDFLAGS_PLE) $(LIBS_PLE) \
$(CCM_LDFLAGS) $(CCMRUNPATH) $(CCM_LIBS) \
$(CGNS_LDFLAGS) $(CGNSRUNPATH) $(CGNS_LIBS) \
$(MED_LDFLAGS) $(MED_LIBS) $(MEDRUNPATH) \
$(HDF5_LDFLAGS) $(HDF5_LIBS) $(HDF5RUNPATH) \
$(METIS_LDFLAGS) $(METIS_LIBS) $(METISRUNPATH) \
$(SCOTCH_LDFLAGS) $(SCOTCH_LIBS) $(SCOTCHRUNPATH) \
$(EOS_LDFLAGS) $(EOS_LIBS) $(EOSRUNPATH) \
$(FREESTEAM_LDFLAGS) $(FREESTEAM_LIBS) $(FREESTEAMRUNPATH) \
$(COOLPROP_LDFLAGS) $(COOLPROP_LIBS) $(COOLPROPRUNPATH) \
$(LDADD_CATALYST) $(LDADD_MEDCOUPLING) $(LDADD_MELISSA) \
$(LDADD_COOLPROP) \
$(MPI_LDFLAGS) $(MPI_LIBS) \
$(LIBXML2_LDFLAGS) $(LIBXML2_LIBS) \
$(LDADD_BLAS) \
$(LTLIBINTL) \
$(FCLIBS)
libsaturne_la_LDFLAGS = -no-undefined -version-info 5:2:0

libcs_solver_la_SOURCES = cs_solver.c
libcs_solver_la_LIBADD = libsaturne.la \
$(LIBPLE_LA) \
$(LDFLAGS_PLE) $(LIBS_PLE) \
$(CCM_LDFLAGS) $(CCMRUNPATH) $(CCM_LIBS) \
$(CGNS_LDFLAGS) $(CGNSRUNPATH) $(CGNS_LIBS) \
$(MED_LDFLAGS) $(MED_LIBS) $(MEDRUNPATH) \
$(HDF5_LDFLAGS) $(HDF5_LIBS) $(HDF5RUNPATH) \
$(METIS_LDFLAGS) $(METIS_LIBS) $(METISRUNPATH) \
$(SCOTCH_LDFLAGS) $(SCOTCH_LIBS) $(SCOTCHRUNPATH) \
$(EOS_LDFLAGS) $(EOS_LIBS) $(EOSRUNPATH) \
$(FREESTEAM_LDFLAGS) $(FREESTEAM_LIBS) $(FREESTEAMRUNPATH) \
$(COOLPROP_LDFLAGS) $(COOLPROP_LIBS) $(COOLPROPRUNPATH) \
$(LDADD_CATALYST) $(LDADD_MEDCOUPLING) $(LDADD_MELISSA) \
$(LDADD_COOLPROP) \
$(MPI_LDFLAGS) $(MPI_LIBS) \
$(LIBXML2_LDFLAGS) $(LIBXML2_LIBS) \
$(LDADD_BLAS) \
$(LTLIBINTL) \
$(FCLIBS)
libcs_solver_la_LDFLAGS = -no-undefined -version-info 5:2:0

endif

# Main program
#pkglibexecdir = $(libexecdir)/$(PACKAGE)
pkglibexec_PROGRAMS =

# Code_Saturne IO utility program

if HAVE_FRONTEND
pkglibexec_PROGRAMS += cs_io_dump
endif

# Code_Saturne syntax checker

if HAVE_FRONTEND
pkglibexec_PROGRAMS += cs_check_syntax
endif

if HAVE_BACKEND
cs_solver$(EXEEXT): libcs_solver.la
	PYTHONPATH=$(top_builddir)/bin:$(top_srcdir)/bin \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
        $(top_srcdir)/src/apps/cs_solver.c
endif

# Code_Saturne IO utility (minimal dependencies)

if HAVE_FRONTEND
cs_io_dump_CPPFLAGS = -DLOCALEDIR=\"'$(localedir)'\" -I$(top_srcdir)/src/base
cs_io_dump_SOURCES = cs_io_dump.c
cs_io_dump_LDADD = $(LTLIBINTL)
endif

# Code_Saturne syntax checker

if HAVE_FRONTEND

cs_check_syntax_CPPFLAGS = \
-DLOCALEDIR=\"'$(localedir)'\" \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/mei
cs_check_syntax_SOURCES = cs_check_syntax.c
cs_check_syntax_LDADD = \
$(top_builddir)/src/mei/libmei.la \
$(top_builddir)/src/bft/libbft.la \
-lm
cs_check_syntax_LDFLAGS = \
$(ORIGINRUNPATH)

endif

# cs_solver executable is built using the same Python script as when
# using user functions (rather than with libtool), so it must be
# handled separately

all-local: cs_solver$(EXEEXT)

install-exec-hook: libcs_solver.la
	PYTHONPATH=$(top_builddir)/bin:$(top_srcdir)/bin \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	 --mode=install --dest "$(DESTDIR)"

uninstall-local:
	rm $(pkglibexecdir)/cs_solver$(EXEEXT)

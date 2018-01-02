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

# Main part

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/alge \
$(FCMODINCLUDE)$(top_builddir)/src/base \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/comb \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCMODINCLUDE)$(top_builddir)/src/turb \
$(FCMODINCLUDE)$(top_builddir)/src/user \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Public header files (to be installed)

pkginclude_HEADERS =

# Library source files

noinst_LTLIBRARIES = libcscomb.la
libcscomb_la_SOURCES = \
cs_coal_bcond.f90 \
cplin1.f90 \
cplini.f90 \
cplph1.f90 \
cplphy.f90 \
cplpro.f90 \
cpltcl.f90 \
cpltss.f90 \
cplvar.f90 \
cpthp1.f90 \
cplver.f90 \
cppdf4.f90 \
cpltsv.f90 \
cplym1.f90 \
cppdfr.f90 \
cpteh1.f90 \
cs_coal_fp2st.f90 \
cs_coal_htconvers1.f90 \
cs_coal_htconvers2.f90 \
cs_coal_masstransfer.f90 \
cs_coal_noxst.f90 \
cs_coal_param.f90 \
cs_coal_physprop1.f90 \
cs_coal_physprop2.f90 \
cs_coal_physprop.f90 \
cs_coal_prop.f90 \
cs_coal_radst.f90 \
cs_coal_readata.f90 \
cs_coal_scast.f90 \
cs_coal_thfieldconv1.f90 \
cs_coal_thfieldconv2.f90 \
cs_coal_varini.f90 \
cs_coal_varpos.f90 \
cs_coal_verify.f90 \
cs_fuel_bcond.f90 \
cs_fuel_fp2st.f90 \
cs_fuel_htconvers1.f90 \
cs_fuel_htconvers2.f90 \
cs_fuel_masstransfer.f90 \
cs_fuel_noxst.f90 \
cs_fuel_param.f90 \
cs_fuel_physprop1.f90 \
cs_fuel_physprop2.f90 \
cs_fuel_physprop.f90 \
cs_fuel_prop.f90 \
cs_fuel_radst.f90 \
cs_fuel_readata.f90 \
cs_fuel_scast.f90 \
cs_fuel_thfieldconv1.f90 \
cs_fuel_thfieldconv2.f90 \
cs_fuel_varini.f90 \
cs_fuel_varpos.f90 \
cs_fuel_verify.f90 \
cs_gascomb.f90 \
cs_physical_properties_combustion_drift.f90 
libcscomb_la_LDFLAGS = -no-undefined

clean-local:
	-rm -f *__genmod.f90 *.mod

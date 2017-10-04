## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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

# Turn off the Automake warning about the makefile portability as we need
# GNUmake rules in order to account for the IDL files target.

AUTOMAKE_OPTIONS = -Wno-portability

include $(top_srcdir)/config/salome_module.mk

MODULE_NAME = CFDSTUDY

BUILT_SOURCES = CFDSTUDY_GenSK.cc
IDL_FILES = CFDSTUDY_Gen.idl

salomelib_LTLIBRARIES = libSalomeIDLCFDSTUDY.la
dist_salomeidl_DATA = $(IDL_FILES) 
libSalomeIDLCFDSTUDY_la_SOURCES =
nodist_libSalomeIDLCFDSTUDY_la_SOURCES = \
	CFDSTUDY_GenSK.cc \
	CFDSTUDY_GenDynSK.cc
nodist_salomeinclude_HEADERS = CFDSTUDY_Gen.hh

libSalomeIDLCFDSTUDY_la_CXXFLAGS = $(SALOME_INCLUDES) $(SALOME_GUI_CPPFLAGS) \
	$(CORBA_CXXFLAGS) \
	$(CORBA_INCLUDES)

libSalomeIDLCFDSTUDY_la_LIBADD  = \
	$(SALOME_GUI_LDFLAGS) -lSalomeIDLGUI \
	$(SALOME_KERNEL_LDFLAGS) -lSalomeIDLKernel \
	$(CORBA_LIBS)

libSalomeIDLCFDSTUDY_la_LDFLAGS = -bindir $(bindir) -no-undefined -version-info=0:0:0

SUFFIXES = .idl .hh SK.cc

%SK.cc %.hh : %.idl
	$(SALOMEENVCMD); \
	$(PYTHON) $(OMNIIDL) -bcxx $(OMNIORB_IDLCXXFLAGS) $(IDL_INCLUDES) $<
%_idl.py : %.idl
	$(SALOMEENVCMD); \
	$(PYTHON) $(OMNIIDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES) $<

CLEANFILES = *.hh *SK.cc *.py *.hxx *.cxx

EXTRA_DIST = $(IDL_FILES)
# potential problem on parallel make on the following - multiple outputs

install-data-local: $(IDL_FILES)
	$(INSTALL) -d  $(DESTDIR)$(salomepythondir)
	$(SALOMEENVCMD); \
	ls $^ | while read file; do \
	  $(OMNIIDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES) -C$(DESTDIR)$(salomepythondir) $$file; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(salomepythondir)/*

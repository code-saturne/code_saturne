#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2010-2011 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Kernel; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
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
	@CORBA_CXXFLAGS@ \
	@CORBA_INCLUDES@

libSalomeIDLCFDSTUDY_la_LIBADD  = \
	$(SALOME_GUI_LDFLAGS) -lSalomeIDLGUI \
	@CORBA_LIBS@

libSalomeIDLCFDSTUDY_la_LDFLAGS = -no-undefined -version-info=0:0:0

SUFFIXES = .idl .hh SK.cc

%SK.cc %.hh : %.idl
	$(OMNIORB_IDL) -bcxx $(OMNIORB_IDLCXXFLAGS) $(IDL_INCLUDES) $<
%_idl.py : %.idl
	$(OMNIORB_IDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES)  $<


CLEANFILES = *.hh *SK.cc *.py *.hxx *.cxx

EXTRA_DIST = $(IDL_FILES)
# potential problem on parallel make on the following - multiple outputs

install-data-local: $(IDL_FILES)
	$(INSTALL) -d  $(DESTDIR)$(salomepythondir)
	ls $^ | while read file; do \
	$(OMNIORB_IDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES) -C$(DESTDIR)$(salomepythondir) $$file ; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(salomepythondir)/*

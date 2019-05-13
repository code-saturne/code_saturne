## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

MODULE_NAME = FSI

BUILT_SOURCES = FSISK.cc
IDL_FILES = FSI.idl

salomelib_LTLIBRARIES = libSalomeIDLFSI.la
salomeidl_DATA = $(IDL_FILES)
libSalomeIDLFSI_la_SOURCES =
nodist_libSalomeIDLFSI_la_SOURCES = FSISK.cc
nodist_salomeinclude_HEADERS = FSI.hh
libSalomeIDLFSI_la_CXXFLAGS = $(SALOME_INCLUDES)
libSalomeIDLFSI_la_LIBADD = $(SALOME_IDL_LIBS)

%SK.cc %.hh : %.idl
	$(SALOMEENVCMD); \
	$(PYTHON) $(OMNIIDL) -bcxx $(OMNIORB_IDLCXXFLAGS) $(IDL_INCLUDES) $<
%_idl.py : %.idl
	$(SALOMEENVCMD); \
	$(PYTHON) $(OMNIIDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES) $<

CLEANFILES = *.hh *SK.cc *.py *.hxx *.cxx

EXTRA_DIST = $(IDL_FILES)

install-data-local: $(IDL_FILES)
	$(INSTALL) -d  $(DESTDIR)$(salomepythondir)
	$(SALOMEENVCMD); \
	ls $^ | while read file; do \
	  $(OMNIIDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES) -C$(DESTDIR)$(salomepythondir) $$file; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(salomepythondir)/*

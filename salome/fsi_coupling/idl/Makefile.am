
include $(top_srcdir)/adm_local/make_common_starter.am

BUILT_SOURCES = FSISK.cc  
IDL_FILES=FSI.idl 

lib_LTLIBRARIES = libSalomeIDLFSI.la
salomeidl_DATA = $(IDL_FILES) 
libSalomeIDLFSI_la_SOURCES      =
nodist_libSalomeIDLFSI_la_SOURCES = FSISK.cc 
nodist_salomeinclude_HEADERS= FSI.hh 
libSalomeIDLFSI_la_CXXFLAGS     = -I.  $(SALOME_INCLUDES)
libSalomeIDLFSI_la_LIBADD     = $(SALOME_IDL_LIBS)
##########################################################
%SK.cc %.hh : %.idl
	$(OMNIORB_IDL) -bcxx $(OMNIORB_IDLCXXFLAGS) $(IDL_INCLUDES) $<
%_idl.py : %.idl
	$(OMNIORB_IDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES)  $<
%PaCO.hxx %PaCO.cxx : %.idl %.xml
	$(OMNIORB_IDL) -I@KERNEL_ROOT_DIR@/idl/salome -p@PACOPATH@/lib/python -bpaco -Wb$(top_srcdir)/idl/$*.xml,$(srcdir):@PACOPATH@/idl:@KERNEL_ROOT_DIR@/idl/salome $(top_srcdir)/idl/$*.idl

CLEANFILES = *.hh *SK.cc *.py *.hxx *.cxx

EXTRA_DIST = $(IDL_FILES)

install-data-local: $(IDL_FILES)
	$(INSTALL) -d  $(DESTDIR)$(salomepythondir)
	ls $^ | while read file; do \
	$(OMNIORB_IDL) $(OMNIORB_IDLPYFLAGS) $(IDL_INCLUDES) -C$(DESTDIR)$(salomepythondir) $$file ; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(salomepythondir)/*



include $(top_srcdir)/adm_local/make_common_starter.am

AM_CFLAGS=$(SALOME_INCLUDES) -fexceptions

lib_LTLIBRARIES= libFSI_MILIEUEngine.la
salomeinclude_HEADERS= FSI_MILIEU.hxx

libFSI_MILIEUEngine_la_SOURCES      = FSI_MILIEU.cxx 
nodist_libFSI_MILIEUEngine_la_SOURCES =
libFSI_MILIEUEngine_la_CXXFLAGS = -I$(top_builddir)/idl  $(SALOME_INCLUDES) -I$(top_srcdir)/salome/libmilieu
libFSI_MILIEUEngine_la_FFLAGS = $(SALOME_INCLUDES) -fexceptions -I$(top_srcdir)/salome/libmilieu
libFSI_MILIEUEngine_la_LIBADD   = $(top_builddir)/salome/libmilieu/libmilieu.la -L$(top_builddir)/idl -lSalomeIDLFSI ${SALOME_LIBS} $(FLIBS)
libFSI_MILIEUEngine_la_LDFLAGS = 


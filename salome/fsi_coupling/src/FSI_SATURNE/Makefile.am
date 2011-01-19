
include $(top_srcdir)/adm_local/make_common_starter.am

AM_CFLAGS=$(SALOME_INCLUDES) -fexceptions

lib_LTLIBRARIES= libFSI_SATURNEEngine.la
salomeinclude_HEADERS= FSI_SATURNE.hxx

libFSI_SATURNEEngine_la_SOURCES      = FSI_SATURNE.cxx 
nodist_libFSI_SATURNEEngine_la_SOURCES =
libFSI_SATURNEEngine_la_CXXFLAGS = -I$(top_builddir)/idl  $(SALOME_INCLUDES) -I$(top_srcdir)/salome/cfd_proxy
libFSI_SATURNEEngine_la_FFLAGS = $(SALOME_INCLUDES) -fexceptions -I$(top_srcdir)/salome/cfd_proxy
libFSI_SATURNEEngine_la_LIBADD   = $(top_builddir)/salome/cfd_proxy/libcfdproxy.la -L$(top_builddir)/idl -lSalomeIDLFSI ${SALOME_LIBS} $(FLIBS)
libFSI_SATURNEEngine_la_LDFLAGS = 


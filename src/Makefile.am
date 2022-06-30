#-------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2022  EDF S.A.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = -DHAVE_CONFIG_H -I$(top_builddir) $(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS) $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = $(LDFLAGS) $(MPI_LDFLAGS)
AM_LIBS = $(MPI_LIBS) $(LIBS)

# Public header files (to be installed)

pleincludedir = $(includedir)
pleinclude_HEADERS = \
ple_coupling.h \
ple_defs.h \
ple_locator.h

noinst_HEADERS = \
ple_config_defs.h

EXTRA_DIST = \
ple_defs.c \
ple_coupling.c \
ple_locator.c

# Dependency rules
#
# - Not handled by AutoMake, since we do not define an execcutable or library
#   in the AutoMake way).
#
# - We could use a finer-grained logic, as only ple_coupling.c and ple_locator.c
#   depend on their respective headers, but since this is a small library
#   that optimization would only save a few seconds of compilation time.
#
# - So we simply use the common suffix rule below.

.c.o: $(pleinclude_HEADERS) $(noinst_HEADERS) ../$(top_builddir)/ple_config.h
	$(CC) $(AM_CPPFLAGS) $(AM_CFLAGS) -c $<

# Library source files

libple_sources = ple_defs.c ple_locator.c
libple_objs = \
ple_defs.o \
ple_locator.o

if HAVE_MPI
libple_sources += ple_coupling.c
libple_objs += \
ple_coupling.o
endif

ple_dep_headers = \
$(pleinclude_HEADERS) \
$(noinst_HEADERS) \
$(top_builddir)/ple_config.h \
$(top_builddir)/ple_config_priv.h

ple_defs.o: $(ple_dep_headers)
ple_locator.o: $(ple_dep_headers)
ple_coupling.o: $(ple_dep_headers)

# Library

libple.a: ${libple_objs}
	@echo "Building libple.a ..."
	rm -f $@
	${AR} $(ARFLAGS) $@ ${libple_objs}
	${RANLIB} $@

libple_soname = libple.${SHARED_LIB_EXT}.$(PLE_LIB_VERSION_STRING_MAJOR)
libple_soname_full = libple.${SHARED_LIB_EXT}.$(PLE_LIB_VERSION_STRING_FULL)
libple_soname_major = libple.${SHARED_LIB_EXT}.$(PLE_LIB_VERSION_STRING_MAJOR)

libple.${SHARED_LIB_EXT}: ${libple_objs}
	@echo "Building libple ..."
	rm -f $@
	rm -rf ${libple_soname_major}
	${PLE_LD} -o ${libple_soname_full} ${LDFLAGS_SONAME}${libple_soname} ${LDFLAGS_SHARED} ${LDFLAGS} ${libple_objs} $(AM_LDFLAGS) $(AM_LIBS)
	ln -s ${libple_soname_full} $@
	ln -s ${libple_soname_full} ${libple_soname_major}

# Local targets

all-local: libple.${LIB_EXT}

clean-local:
	@if test "${enable_shared}" = "yes" ; then \
	  rm -f "$(libple_soname_full) $(libple_soname_major)"; \
	fi;
	-rm -f libple.${LIB_EXT}
	-rm -f *.o

install-exec-local: libple.${LIB_EXT}
	@libraries=libple.${LIB_EXT}; \
	if test "${enable_shared}" = "yes" ; then \
	  libraries="$$libraries $(libple_soname_full) $(libple_soname_major)"; \
	fi; \
	list=; for p in $$libraries; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	test -z "$$list" -a -n "$(libdir)" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(libdir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(libdir)" || exit 1; \
	  echo " cp -fp $$list '$(DESTDIR)$(libdir)'"; \
	  cp -fpP $$list "$(DESTDIR)$(libdir)"; \
	}

uninstall-local:
	@libraries=libple.${LIB_EXT}; \
	if test "${enable_shared}" = "yes" ; then \
	  libraries="$$libraries $(libple_soname_full) $(libple_soname_major)"; \
	fi; \
	list=; for p in $$libraries; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	for p in $$list; do \
	  echo " rm -f '$(DESTDIR)$(libdir)/$$p'"; \
	  rm -f "$(DESTDIR)$(libdir)/$$p"; \
	done; \
	rmdir "$(DESTDIR)$(libdir)" || exit 0

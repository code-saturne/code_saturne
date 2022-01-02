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

if HAVE_MPI4PY

EXTRA_DIST = pyple_coupler.py.in

edit_python = sed \
	-e 's|@mpi4pydir[@]|$(MPI4PY_PATH)|g'

pyple_coupler.py: Makefile
	rm -f $@ $@.tmp
	srcidr=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_python) $${srcdir}$@.in >$@.tmp;
	mv $@.tmp $@

AM_CPPFLAGS = -DHAVE_CONFIG_H -I$(top_builddir) $(MPI_CPPFLAGS)
AM_CFLAGS = $(CFLAGS) $(CFLAGS_DBG) $(CFLAGS_OPT)
AM_LDFLAGS = $(LDFLAGS) $(MPI_LDFLAGS)
AM_LIBS = $(MPI_LIBS) $(LIBS) $(top_builddir)/src/libple.${LIB_EXT}

# Build rules

ple_coupling_wrapper.o: ple_coupling_wrapper.c $(top_builddir)/src/libple.${LIB_EXT}
	$(CC) $(PYTHON_CPPFLAGS) $(AM_CPPFLAGS) $(MPI4PY_CFLAGS) $(AM_CFLAGS) -I$(srcdir)/../src -c $(srcdir)/ple_coupling_wrapper.c

# Base Library

libpyplecoupling.a: ple_coupling_wrapper.o
	@echo "Building libple.a ..."
	rm -f $@
	${AR} $(ARFLAGS) $@ ple_coupling_wrapper.o
	${RANLIB} $@

libpyplecoupling_soname = libpyplecoupling.${SHARED_LIB_EXT}.$(PLE_LIB_VERSION_STRING_MAJOR)

libpyplecoupling.${SHARED_LIB_EXT}: ple_coupling_wrapper.o
	@echo "Building libpyplecoupling ..."
	@link_cmd="${PLE_LD} -o $@ ${LDFLAGS_SHARED} $(AM_LDFLAGS) ple_coupling_wrapper.o $(top_builddir)/src/libple.${LIB_EXT} $(AM_LIBS) $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS)" ; \
	if test -n "$(LDRPATH)"; then \
	  echo "$${link_cmd} $(LDRPATH)'$(abs_top_builddir)/src'" ; \
	  $${link_cmd} $(LDRPATH)'$(abs_top_builddir)/src' ; \
	else \
	  echo "$${link_cmd}" ; \
	  $${link_cmd} ; \
	fi;

pkgpython_PYTHON = __init__.py Coupling.py Init.py
nodist_pkgpython_PYTHON = pyple_coupler.py

DISTCLEANFILES = pyple_coupler.py

# Local targets

all-local: libpyplecoupling.${LIB_EXT} pyple_coupler.py

clean-local:
	-rm -f pyple_coupler.py
	-rm -f libpyplecoupling.${LIB_EXT}
	-rm -f *.o

install-exec-local: libpyplecoupling.${LIB_EXT} pyple_coupler.py
	@test -z -n "$(libdir)" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(pkgpythondir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(pkgpythondir)" || exit 1; \
	  if test "${enable_shared}" = "yes" -a test -n "$(LDRPATH)"; then \
	    echo "Installing (relinking) libpyplecoupling ..."; \
	    link_cmd="${PLE_LD} '$(DESTDIR)$(pkgpythondir)/libpyplecoupling.${LIB_EXT}' ${LDFLAGS_SHARED} $(AM_LDFLAGS) ple_coupling_wrapper.o $(top_builddir)/src/libple.${LIB_EXT} $(AM_LIBS) $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS) $(LDRPATH)'$(DESTDIR)$(libdir)'" ; \
	    echo "$${link_cmd}" ; \
	    $${link_cmd} ; \
	  else \
	    echo "cp -fpP libpyplecoupling.${LIB_EXT} '$(DESTDIR)$(pkgpythondir)'"; \
	    cp -fpP libpyplecoupling.${LIB_EXT} '$(DESTDIR)$(pkgpythondir)'; \
	  fi ; \
	}

uninstall-local:
	echo "lib a"
	-rm -rf "$(DESTDIR)$(pkgpythondir)/__pycache__";
	-rm -f "$(DESTDIR)$(pkgpythondir)/libpyplecoupling.${LIB_EXT}";
	echo "lib b"
	rmdir -p "$(DESTDIR)$(pkgpythondir)" || exit 0

endif


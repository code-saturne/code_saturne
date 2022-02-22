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

# The name of the module.
DOC_MODULE=saturne

docdir = $(datadir)/doc/@PACKAGE@
EXTRA_DIST = Doxyfile \
mainpage.dox \
ple_coupling_example_domains.svg \
ple_coupling_example_exchange.svg \
ple_coupling_example_init.svg

if HAVE_DOXYGEN

doxygen:
	srcdir="$(top_srcdir)" builddir="$(top_builddir)" "$(DOXYGEN)" $(top_srcdir)/doc/doxygen/Doxyfile

doxygen-install: install-html-local

html-local: doxygen

clean-local:
	rm -rf index.bt

endif

distclean-local:
	rm -rf *~
	rm -rf html latex man index.bt
	rm -rf Makefile

# Generate the documentation at dist stage so that the user do not need
# to have Doxygen (and mandatory tools).

if HAVE_DOXYGEN

dist-hook: html
	\cp -R html $(distdir)/

install-html-local:
	$(MKDIR_P) "$(DESTDIR)$(htmldir)/doxygen"
	\cp -R html/* $(DESTDIR)$(htmldir)/doxygen/

endif

uninstall-local:
	-rm -rf $(htmldir)/doxygen


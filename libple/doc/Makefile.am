# Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2010  EDF
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
DOC_MODULE=fvm

docdir = $(datadir)/doc/@PACKAGE@-@VERSION@
doc_DATA = configuration  internationalization
EXTRA_DIST = $(doc_DATA)

doxygen:
	srcdir="$(top_srcdir)" builddir="$(top_builddir)" doxygen $(top_srcdir)/doc/Doxyfile

doxygen-install: doxygen
	cp -r html $(docdir)

uninstall:
	rm -rf $(docdir)

clean:
	rm -rf index.bt

distclean:
	rm -rf *~
	rm -rf html latex man index.bt
	rm -rf Makefile

uninstall-local:
	rm -rf $(docdir)



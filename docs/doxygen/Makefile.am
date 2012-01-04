## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

EXTRA_DIST = Doxyfile-pre Doxyfile-src Doxyfile-ui \
mainpage.h images/icon.png examples/bft_error_example.c

# Environment variables for LaTeX

TEXINPUTS = $(top_srcdir)/docs/style/:

doxygen-pre:
	srcdir="$(top_srcdir)" builddir="$(top_builddir)" \
	doxygen $(top_srcdir)/docs/doxygen/Doxyfile-pre

doxygen-src:
	srcdir="$(top_srcdir)" builddir="$(top_builddir)" \
	doxygen $(top_srcdir)/docs/doxygen/Doxyfile-src

doxygen-ui:
	srcdir="$(top_srcdir)" builddir="$(top_builddir)" \
	doxygen $(top_srcdir)/docs/doxygen/Doxyfile-ui

doxygen: doxygen-pre doxygen-src doxygen-ui

html-local: doxygen

clean-local:
	for d in pre src ui ;\
	do \
	  if test -d $$d ;\
	  then \
	    (cd $$d && rm -rf html latex man index.bt) ;\
          fi ;\
	done

distclean-local:
	-rm -rf Makefile

install-html-local:
	$(MKDIR_P) "$(DESTDIR)$(htmldir)/doxygen"
	\cp -R pre/html $(DESTDIR)$(htmldir)/doxygen/pre
	\cp -R src/html $(DESTDIR)$(htmldir)/doxygen/src
	\cp -R ui/html  $(DESTDIR)$(htmldir)/doxygen/ui

uninstall-local:
	-rm -rf $(htmldir)/doxygen

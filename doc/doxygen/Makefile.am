## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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

EXTRA_DIST = Doxyfile mainpage.h images/icon.png examples/bft_error_example.c

doxygen:
	srcdir="$(top_srcdir)" builddir="$(top_builddir)" doxygen $(top_srcdir)/doc/doxygen/Doxyfile

html-local: doxygen

clean-local:
	-rm -rf html latex man index.bt

distclean-local:
	-rm -rf Makefile

install-html-local:
	$(MKDIR_P) "$(DESTDIR)$(htmldir)"
	\cp -R html $(DESTDIR)$(htmldir)/doxygen

uninstall-local:
	-rm -rf $(docdir)

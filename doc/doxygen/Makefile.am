#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne CFD tool.
#
#   Copyright (C) 2011 EDF S.A., France
#
#   The Code_Saturne CFD tool is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne CFD tool is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne CFD tool; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

EXTRA_DIST = Doxyfile mainpage.h images/icon.png

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

## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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

SUBDIRS = . model \
bft fvm gui \
base cdo pprt alge mesh turb \
atmo cfbl cogz comb ctwr elec gwf lagr rayt \
user user_examples meg apps

# General flags

AM_LDFLAGS = -no-undefined

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_headers.h

# Additional local targets

install-exec-local:
	@$(NORMAL_INSTALL)
	test -z "$(pkgincludedir)" || $(MKDIR_P) "$(DESTDIR)$(pkgincludedir)"
	test -h $(DESTDIR)$(pkgincludedir)/cs_defs.h || ln -s base/cs_defs.h $(DESTDIR)$(pkgincludedir)/cs_defs.h;

uninstall-local:
	@$(NORMAL_UNINSTALL)
	rm -f $(DESTDIR)$(pkgincludedir)/cs_defs.h;
	rmdir $(pkgincludedir) || exit 0

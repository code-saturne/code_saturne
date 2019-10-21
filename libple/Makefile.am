#-------------------------------------------------------------------------------
#   This file is part of the "Parallel Location and Exchange" library,
#   intended to provide mesh or particle-based code coupling services.
#
#   Copyright (C) 2005-2019  EDF S.A.
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

AUX_DIST = build-aux/compile \
           build-aux/config.guess \
           build-aux/config.sub \
           build-aux/depcomp \
           build-aux/install-sh \
           build-aux/ltmain.sh \
           build-aux/missing \
           m4/libtool.m4 \
           m4/lt~obsolete.m4 \
           m4/ltoptions.m4 \
           m4/ltsugar.m4 \
           m4/ltversion.m4

MAINTAINERCLEANFILES = aclocal.m4 configure $(AUX_DIST)

SUBDIRS = src pyple tests doc po

EXTRA_DIST = m4/ple_config_info.m4 \
             m4/ple_config_publ.m4 \
             m4/ple_mpi.m4 \
             m4/ple_mpi4py.m4 \
             m4/ple_types.m4 \
             config/ple_auto_flags.sh \
             sbin/backup \
             sbin/bootstrap \
             sbin/clean \
             sbin/rmb

# One adds the examples sub-directories with a dist-hook.
# The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += examples

bin_SCRIPTS = ple-config

# Install man pages

dist_man_MANS = \
doc/manpages/ple-config.1

configincludedir = $(includedir)
nodist_configinclude_HEADERS = ple_config.h

nodist_noinst_HEADERS = ple_config_priv.h

DISTCLEANFILES = ple-config ple_config.h ple_config_priv.h

ACLOCAL_AMFLAGS = -I m4

uninstall-local:
	rmdir $(configincludedir) || exit 0

# Update gettext dictionnaries

update-po:
	cd po && $(MAKE) $(AM_MAKEFLAGS) update-po

# Generate and install documentation

.PHONY: doc
doc:
	cd doc && $(MAKE) $(AM_MAKEFLAGS) doc

.PHONY: install-doc
install-doc:
	cd doc && $(MAKE) $(AM_MAKEFLAGS) install-doc

install-tests:
	@cd tests ; $(MAKE) install-tests

uninstall-tests:
	@cd tests ; $(MAKE) uninstall-tests

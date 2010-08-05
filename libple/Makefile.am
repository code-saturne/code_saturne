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

AUX_DIST = $(ac_aux_dir)/config.guess \
           $(ac_aux_dir)/config.sub \
           $(ac_aux_dir)/depcomp \
           $(ac_aux_dir)/install-sh \
           $(ac_aux_dir)/missing \
           $(ac_aux_dir)/mkinstalldirs \
           $(ac_aux_dir)/ltmain.sh

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure config-h.in stamp-h.in $(AUX_DIST)

# SUBDIRS = src tests doc intl m4 po

SUBDIRS = src tests doc

#EXTRA_SUBDIRS = include

EXTRA_DIST = $(ac_aux_dir)/ple_auto_flags.sh \
             $(ac_aux_dir)/ple_config_info.m4 \
             $(ac_aux_dir)/ple_config_publ.m4 \
             $(ac_aux_dir)/ple_mpi.m4 \
             $(ac_aux_dir)/ple_types.m4

bin_SCRIPTS = ple-config

# Install man pages

dist_man_MANS = \
doc/manpages/ple-config.1

configincludedir = $(includedir)
configinclude_HEADERS = ple_config.h

noinst_HEADERS = ple_config_priv.h

uninstall-local:
	rmdir $(configincludedir) || exit 0

DISTCLEANFILES = ple-config ple_config.h ple_config_priv.h

ACLOCAL_AMFLAGS = -I config

install-tests:
	@cd tests ; $(MAKE) install-tests

uninstall-tests:
	@cd tests ; $(MAKE) uninstall-tests

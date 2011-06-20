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
           $(ac_aux_dir)/ltmain.sh \
           $(ac_aux_dir)/libtool.m4 \
           $(ac_aux_dir)/lt~obsolete.m4 \
           $(ac_aux_dir)/ltoptions.m4 \
           $(ac_aux_dir)/ltsugar.m4

MAINTAINERCLEANFILES = aclocal.m4 configure $(AUX_DIST)

SUBDIRS = src tests doc

EXTRA_DIST = $(ac_aux_dir)/ple_auto_flags.sh \
             $(ac_aux_dir)/ple_config_info.m4 \
             $(ac_aux_dir)/ple_config_publ.m4 \
             $(ac_aux_dir)/ple_mpi.m4 \
             $(ac_aux_dir)/ple_types.m4 \
             sbin/backup \
             sbin/bootstrap \
             sbin/clean \
             sbin/rmb

# One adds the examples sub-directories and get rid of the .svn
# directories with a dist-hook. The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += examples

bin_SCRIPTS = ple-config

# Install man pages

dist_man_MANS = \
doc/manpages/ple-config.1

configincludedir = $(includedir)
nodist_configinclude_HEADERS = ple_config.h

nodist_noinst_HEADERS = ple_config_priv.h

uninstall-local:
	rmdir $(configincludedir) || exit 0

DISTCLEANFILES = ple-config ple_config.h ple_config_priv.h

ACLOCAL_AMFLAGS = -I m4

install-tests:
	@cd tests ; $(MAKE) install-tests

uninstall-tests:
	@cd tests ; $(MAKE) uninstall-tests

dist-hook:
	-rm -rf `find $(distdir)/examples -name .svn`

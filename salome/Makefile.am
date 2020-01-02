## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

ACLOCAL_AMFLAGS = -I m4

if HAVE_SALOME_GUI
CFDSTUDY = cfd_study
endif

SUBDIRS = $(CFDSTUDY)

AUX_DIST = \
build-aux/list_salome_modules.sh \
build-aux/compile \
build-aux/config.guess \
build-aux/config.sub \
build-aux/depcomp \
build-aux/install-sh \
build-aux/ltmain.sh \
build-aux/missing \
build-aux/py-compile \
m4/libtool.m4 \
m4/lt~obsolete.m4 \
m4/ltoptions.m4 \
m4/ltsugar.m4 \
m4/ltversion.m4

MAINTAINERCLEANFILES = \
$(AUX_DIST) \
aclocal.m4 \
configure \
build-aux/compile \
build-aux/config.guess \
build-aux/config.sub \
build-aux/depcomp \
build-aux/install-sh \
build-aux/ltmain.sh \
build-aux/missing \
build-aux/py-compile

EXTRA_DIST = \
$(AUX_DIST) \
m4/cs_omniorb.m4 \
m4/cs_salome.m4 \
README.md \
sbin/backup \
sbin/bootstrap \
sbin/clean \
sbin/rmb \
config/cs_config_test/py\
config/salome_module.mk

DISTCLEANFILES = \
bin/salome_cfd

if HAVE_RELOCATABLE
  PYTHON_SHEBANG = /usr/bin/env python3
else
  PYTHON_SHEBANG = $(PYTHON)
endif

# Replace paths in Python scripts.
# The following command and target is adapted from an example in
# autoconf. The autoconf documentation recommends changing scripts
# through Makefile targets rather than through configure, and
# cites its own automate and autoheader as examples.

edit_package = sed \
	-e 's|@salome_env_modules[@]|$(salome_env_modules)|g' \
	-e 's|@SALOMEENVCMD[@]|$(SALOMEENVCMD)|g' \
	-e 's|@PYTHON[@]|$(PYTHON_SHEBANG)|g' \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@pythondir[@]|$(pythondir)|g' \
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g' \
	-e 's|@docdir[@]|$(docdir)|g' \
	-e 's|@relocatable[@]|$(relocatable)|g'

# Install Python library

bin/salome_cfd: Makefile $(srcdir)/bin/salome_cfd.in
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_package) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

# Install binary scripts

nodist_bin_SCRIPTS = \
bin/salome_cfd

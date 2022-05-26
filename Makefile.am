## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

if HAVE_FRONTEND
PREPRO = preprocessor
DOCS = docs
endif

if HAVE_BACKEND
TSTS = tests
endif

if HAVE_INTERNAL_PLE
PLE = libple
endif

SUBDIRS = . $(PLE) $(PREPRO) src python $(DOCS) $(TSTS) $(NEPTUNE_CFD)

AUX_DIST =

MAINTAINERCLEANFILES = \
$(AUX_DIST) \
aclocal.m4 \
configure \
build-aux/ar-lib \
build-aux/compile \
build-aux/config.guess \
build-aux/config.rpath \
build-aux/config.sub \
build-aux/depcomp \
build-aux/install-sh \
build-aux/ltmain.sh \
build-aux/missing \
build-aux/py-compile \
build-aux/ylwrap

EXTRA_DIST = \
$(AUX_DIST) \
install_saturne.py \
INSTALL.md \
NEWS.md \
README.md \
QUALITY_ASSURANCE \
build-aux/catalyst \
build-aux/catalyst-5.6 \
build-aux/cs_build_tools.py \
build-aux/cs_compile_build.py \
build-aux/cs_config_test.py \
build-aux/cs_user_src_filter.py \
build-aux/cs_version.py \
build-aux/list_modules.sh \
build-aux/petsc-petscvariables.makefile \
build-aux/petsc-variables.makefile \
sbin/backup \
sbin/bootstrap \
sbin/clean \
sbin/rmb \
config/cs_auto_flags.sh \
config/code_saturne_build.cfg.in.in

# One adds several sub-directories with a dist-hook.
# The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += extras m4 patches

nodist_pkginclude_HEADERS = \
cs_config.h \
cs_version.h

DISTCLEANFILES = \
bin/code_saturne \
lib/code_saturne_build.cfg \
cs_config.h \
cs_version.h \
libtool_cuda

ACLOCAL_AMFLAGS = -I m4

if HAVE_RELOCATABLE
  PYTHON_SHEBANG = /usr/bin/env python
else
  PYTHON_SHEBANG = $(PYTHON)
endif

# Replace paths in Python scripts.
# The following command and targets are adapted from an example in
# autoconf. The autoconf documentation recommends changing scripts
# through Makefile targets rather than through configure, and
# cites its own automate and autoheader as examples.

edit_python = sed \
	-e 's|@PYTHON[@]|$(PYTHON_SHEBANG)|g' \
	-e 's|@libdir[@]|$(libdir)|g' \
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g' \
	-e 's|@relocatable[@]|$(relocatable)|g'

edit_cfg = sed \
	-e 's|@cfg_prefix[@]|$(prefix)|g' \
	-e 's|@cfg_exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@cfg_bindir[@]|$(bindir)|g' \
	-e 's|@cfg_includedir[@]|$(includedir)|g' \
	-e 's|@cfg_pkgincludedir[@]|$(pkgincludedir)|g' \
	-e 's|@cfg_libdir[@]|$(libdir)|g' \
	-e 's|@cfg_libexecdir[@]|$(libexecdir)|g' \
	-e 's|@cfg_pkglibexecdir[@]|$(pkglibexecdir)|g' \
	-e 's|@cfg_pythondir[@]|$(pythondir)|g' \
	-e 's|@cfg_pkgpythondir[@]|$(pkgpythondir)|g' \
	-e 's|@cfg_localedir[@]|$(localedir)|g' \
	-e 's|@cfg_datarootdir[@]|$(datarootdir)|g' \
	-e 's|@cfg_datadir[@]|$(datadir)|g' \
	-e 's|@cfg_pkgdatadir[@]|$(pkgdatadir)|g' \
	-e 's|@cfg_docdir[@]|$(docdir)|g' \
	-e 's|@cfg_pdfdir[@]|$(pdfdir)|g' \
	-e 's|@cfg_sysconfdir[@]|$(sysconfdir)|g' \
	-e 's|@PACKAGE_NAME[@]|$(PACKAGE_NAME)|g' \
	-e 's|@PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	-e 's|@PACKAGE_STRING[@]|$(PACKAGE_STRING)|g' \
	-e 's|@PACKAGE_BUGREPORT[@]|$(PACKAGE_BUGREPORT)|g' \
	-e 's|@PACKAGE_URL[@]|$(PACKAGE_URL)|g'

edit_shell = sed \
	-e 's|@SHELL[@]|$(SHELL)|g' \
	-e 's|@PYTHON[@]|$(PYTHON)|g' \
	-e 's|@cs_env_modules[@]|$(cs_env_modules)|g' \
	-e 's|@shell_env[@]|$(shell_env)|g'

if HAVE_SHELL_LAUNCHER
  LAUNCHER_EXT=sh
  edit_launcher=$(edit_shell)
else
  LAUNCHER_EXT=py
  edit_launcher=$(edit_python)
endif

bin/code_saturne: Makefile
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.$(LAUNCHER_EXT).in || srcdir=$(srcdir)/; \
	$(edit_launcher) $${srcdir}$@.$(LAUNCHER_EXT).in >$@.tmp
	mv $@.tmp $@
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace $@ $@

if HAVE_SHELL_LAUNCHER
bin/code_saturne.py: Makefile $(srcdir)/bin/code_saturne.py.in
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_python) $${srcdir}$@.in >$@.tmp;
	mv $@.tmp $@
endif

bin/code_saturne: $(srcdir)/bin/code_saturne.$(LAUNCHER_EXT).in

cs_version.h: Makefile $(top_srcdir)/src/base/cs_version.h.in
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace ${top_srcdir}/src/base/cs_version.h.in cs_version.h

lib/code_saturne_build.cfg: Makefile config/code_saturne_build.cfg.in
	rm -f $@ $@.tmp
	$(MKDIR_P) lib
	srcdir=''; \
	$(edit_cfg) config/code_saturne_build.cfg.in >$@.tmp
	mv $@.tmp $@
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace $@ $@

BUILT_SOURCES = \
cs_version.h \
lib/code_saturne_build.cfg

CLEANFILES = $(edit_scripts)
EXTRA_DIST += bin/code_saturne.py.in \
	      bin/code_saturne.sh.in \
              data/thch/dp_radiat_MFS.tar.xz

# Install binary scripts

nodist_bin_SCRIPTS = \
bin/code_saturne

if HAVE_SHELL_LAUNCHER
nodist_bin_SCRIPTS += bin/code_saturne.py
endif

# Install Python library

# Install user scripts

dist_pkgdata_DATA = \
COPYING

dist_pkgdata_SCRIPTS = \
bin/runcase_mpi_rank

# Install thermochemistry and other data files

nobase_dist_pkgdata_DATA = \
data/thch/dp_C3P_auto \
data/thch/dp_radiat_ADF50 \
data/thch/dp_radiat_ADF8 \
data/thch/dp_radiat_MFS.tar.xz \
data/thch/dp_radiat_Planck_CO2 \
data/thch/dp_radiat_Planck_H2O \
data/thch/dp_radiat_wave \
data/thch/JANAF \
data/user/cs_user_scripts.py \
data/user/dp_C3P \
data/user/dp_C3PSJ \
data/user/dp_C4P \
data/user/dp_ELE \
data/user/dp_FUE \
data/user/dp_transformers \
data/user/meteo \
data/setup.xml

# Install batch cards for cluster management

batchdir = ${pkgdatadir}/batch
batch_DATA = \
extras/batch/batch.CCC \
extras/batch/batch.LSF \
extras/batch/batch.OAR \
extras/batch/batch.PBS \
extras/batch/batch.SGE \
extras/batch/batch.SLURM \
extras/batch/batch_help.CCC \
extras/batch/batch_help.LSF \
extras/batch/batch_help.OAR \
extras/batch/batch_help.PBS \
extras/batch/batch_help.SGE \
extras/batch/batch_help.SLURM

# Install man pages

dist_man_MANS = \
docs/manpages/code_saturne.1

# Install bash completion script

bashcompletiondir = ${sysconfdir}/bash_completion.d
bashcompletion_DATA = \
extras/bash_completion/code_saturne

# Install global option file

sysconf_DATA = \
extras/code_saturne.cfg.template

# Generate and install documentation

doc_dirs = $(PLE) docs $(NEPTUNE_CFD)

doc:
	for p in $(doc_dirs); do \
	  cd $$p && $(MAKE) $(AM_MAKEFLAGS) doc ;\
	  cd $(abs_top_builddir) ;\
	done

install-doc:
	for p in $(doc_dirs); do \
	  cd $$p && $(MAKE) $(AM_MAKEFLAGS) install-doc ;\
	  cd $(abs_top_builddir) ;\
	done

# Update Copyright

update-name:
	lst=`find ${top_srcdir} -type d \( -name .git -o -name "build*" \) -prune -o \( -type f -a ! -name Makefile.in \) -print`; \
	for f in $$lst; do \
	  echo $$f && \
	  $(SED) -i "s/code_saturne documentation page/code_saturne documentation page/" $$f; \
	done

update-copyright:
	y=`date +%Y`; yold=$$(($$y-1)); \
	lst=`find ${top_srcdir} -type d \( -name .git -o -name "build*" \) -prune -o \( -type f -a ! -name Makefile.in \) -print`; \
	for f in $$lst; do \
	  echo $$f && \
	  $(SED) -i "s/Copyright (C) 1998-$${yold} EDF S.A./Copyright (C) 1998-$$y EDF S.A./" $$f; \
	done

dist-hook:
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --full > $(distdir)/build-aux/version

# Additional local targets

all-local:
	@$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace ${top_srcdir}/src/base/cs_version.h.in cs_version.h
	@if test -f bin/cs_package.py; then \
	  $(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace bin/cs_package.py bin/cs_package.py; \
	fi

install-data-local:
	@$(NORMAL_INSTALL)
	$(MKDIR_P) "$(DESTDIR)$(libdir)"
	cp lib/code_saturne_build.cfg "$(DESTDIR)$(libdir)/code_saturne_build.cfg"
	test -z "$(pkgdatadir)/data/thch" || $(MKDIR_P) "$(DESTDIR)$(pkgdatadir)/data/thch"
	tar xvJf "$(top_srcdir)/data/thch/dp_radiat_MFS.tar.xz" -C "$(DESTDIR)$(pkgdatadir)/data/thch"

uninstall-local:
	rm -rf "$(DESTDIR)$(pkgdatadir)/data/thch/dp_radiat_MFS"
	rm -f"$(DESTDIR)$(libdir)"
	@$(NORMAL_UNINSTALL)

install-shortcut:
	@$(NORMAL_INSTALL)
	@if test "${cs_have_gui}" != "yes"; then \
		echo "Cannot install shortcut without GUI support"; \
		exit 1 ; \
	fi; \
	DESKTOPPATH=`xdg-user-dir DESKTOP 2>/dev/null` ; \
	if test "$${DESKTOPPATH}" = "" ; then \
		if test -d "$${HOME}/Desktop" ; then \
		  DESKTOPPATH=$${HOME}/Desktop ; \
		elif test -d "$${HOME}/Bureau" ; then \
		  DESKTOPPATH=$${HOME}/Bureau ; \
		else \
			DESKTOPPATH="none" ; \
		fi ; \
	fi ; \
	if test "$${DESKTOPPATH}" = "none" ; then \
		echo "Could not find DESKTOPPATH" ; \
		exit 1 ; \
	fi ; \
	if test ! -d "$${DESKTOPPATH}" ; then \
		echo "'$${DESKTOPPATH}' is not a folder" ; \
		exit 1 ; \
	fi ; \
	if test "${shortcut-solver}" != "" ; then \
	  DESKTOPSOLVER="${shortcut-solver}" ; \
	else \
	  DESKTOPSOLVER=code_saturne ; \
	fi ; \
	cp "$(top_srcdir)/extras/code_saturne.desktop.template" "$${DESKTOPPATH}/code_saturne_$(VERSION).desktop"; \
	sed -i -e 's|CSPATH|$(prefix)|g' $${DESKTOPPATH}/code_saturne_$(VERSION).desktop; \
	sed -i -e 's|CSVERSION|$(VERSION)|g' $${DESKTOPPATH}/code_saturne_$(VERSION).desktop ; \
	sed -i -e 's|DSOLVER|'$${DESKTOPSOLVER}'|g' $${DESKTOPPATH}/code_saturne_$(VERSION).desktop ;

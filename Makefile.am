## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

if HAVE_GUI
GUI = gui
endif

if HAVE_FRONTEND
PREPRO = preprocessor
endif

if HAVE_INTERNAL_PLE
PLE = libple
endif

if HAVE_SALOME
SALOME = salome
endif

SUBDIRS = . $(PLE) $(PREPRO) src $(SALOME) po docs $(GUI) tests

AUX_DIST = \
build-aux/catalyst \
build-aux/cs_config_test.py \
build-aux/cs_version.py \
m4/codeset.m4 \
m4/fcntl-o.m4 \
m4/gettext.m4 \
m4/glibc21.m4 \
m4/glibc2.m4 \
m4/iconv.m4 \
m4/intdiv0.m4 \
m4/intldir.m4 \
m4/intl.m4 \
m4/intlmacosx.m4 \
m4/intmax.m4 \
m4/inttypes_h.m4 \
m4/inttypes-pri.m4 \
m4/lcmessage.m4 \
m4/lib-ld.m4 \
m4/lib-link.m4 \
m4/lib-prefix.m4 \
m4/libtool.m4 \
m4/lock.m4 \
m4/longlong.m4 \
m4/lt~obsolete.m4 \
m4/ltoptions.m4 \
m4/ltsugar.m4 \
m4/ltversion.m4 \
m4/nls.m4 \
m4/po.m4 \
m4/printf-posix.m4 \
m4/progtest.m4 \
m4/size_max.m4 \
m4/stdint_h.m4 \
m4/threadlib.m4 \
m4/uintmax_t.m4 \
m4/visibility.m4 \
m4/wchar_t.m4 \
m4/wint_t.m4 \
m4/xsize.m4

MAINTAINERCLEANFILES = \
$(AUX_DIST) \
aclocal.m4 \
configure \
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
QUALITY_ASSURANCE \
po/Makevars.template \
po/Rules-install \
sbin/backup \
sbin/bootstrap \
sbin/clean \
sbin/rmb \
config/cs_auto_flags.sh \
m4/acx_pthread.m4 \
m4/cs_blas.m4 \
m4/cs_catalyst.m4 \
m4/cs_ccm.m4 \
m4/cs_cgns.m4 \
m4/cs_docs.m4 \
m4/cs_freesteam.m4 \
m4/cs_fortran.m4 \
m4/cs_hdf5.m4 \
m4/cs_libxml2.m4 \
m4/cs_med.m4 \
m4/cs_metis.m4 \
m4/cs_modules.m4 \
m4/cs_mpi.m4 \
m4/cs_omniorb.m4 \
m4/cs_ple.m4 \
m4/cs_salome.m4 \
m4/cs_scotch.m4 \
m4/cs_zlib.m4

# One adds several sub-directories and get rid of the .svn
# directories with a dist-hook. The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += examples extras patches

nodist_pkginclude_HEADERS = cs_config.h

DISTCLEANFILES = \
bin/code_saturne \
bin/cs_package.py \
extras/salome/fsi_appli_config.xml \
cs_config.h

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
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g' \
	-e 's|@relocatable[@]|$(relocatable)|g'

edit_package = sed \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@bindir[@]|$(bindir)|g' \
	-e 's|@includedir[@]|$(includedir)|g' \
	-e 's|@pkgincludedir[@]|$(pkgincludedir)|g' \
	-e 's|@libdir[@]|$(libdir)|g' \
	-e 's|@libexecdir[@]|$(libexecdir)|g' \
	-e 's|@pkglibexecdir[@]|$(pkglibexecdir)|g' \
	-e 's|@pythondir[@]|$(pythondir)|g' \
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g' \
	-e 's|@localedir[@]|$(localedir)|g' \
	-e 's|@datarootdir[@]|$(datarootdir)|g' \
	-e 's|@datadir[@]|$(datadir)|g' \
	-e 's|@pkgdatadir[@]|$(pkgdatadir)|g' \
	-e 's|@docdir[@]|$(docdir)|g' \
	-e 's|@pdfdir[@]|$(pdfdir)|g' \
	-e 's|@sysconfdir[@]|$(sysconfdir)|g' \
	-e 's|@PACKAGE_NAME[@]|$(PACKAGE_NAME)|g' \
	-e 's|@PACKAGE_TARNAME[@]|$(PACKAGE_TARNAME)|g' \
	-e 's|@PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	-e 's|@PACKAGE_STRING[@]|$(PACKAGE_STRING)|g' \
	-e 's|@PACKAGE_BUGREPORT[@]|$(PACKAGE_BUGREPORT)|g' \
	-e 's|@PACKAGE_URL[@]|$(PACKAGE_URL)|g' \
	-e 's|@cs_src_revision[@]|$(cs_src_revision)|g' \
	-e 's|@cs_version[@]|$(cs_version)|g' \
	-e 's|@cs_version_full[@]|$(cs_version_full)|g' \
	-e 's|@cs_version_short[@]|$(cs_version_short)|g' \
	-e 's|@SALOME_KERNEL[@]|$(SALOME_KERNEL)|g' \
	-e 's|@SALOMEENVCMD[@]|$(SALOMEENVCMD)|g'

bin/code_saturne: Makefile
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_python) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

bin/cs_package.py: Makefile
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_package) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

bin/code_saturne: $(srcdir)/bin/code_saturne.in
bin/cs_package.py: $(srcdir)/bin/cs_package.py.in

CLEANFILES = $(edit_scripts)
EXTRA_DIST += bin/code_saturne.in bin/cs_package.py.in

# Install binary scripts

nodist_bin_SCRIPTS = \
bin/code_saturne

# Install Python library

nodist_pkgpython_PYTHON = \
bin/cs_config.py \
bin/cs_package.py

dist_pkgpython_PYTHON = \
bin/__init__.py \
bin/cs_autovnv.py \
bin/cs_bdiff.py \
bin/cs_bdump.py \
bin/cs_case.py \
bin/cs_case_domain.py \
bin/cs_case_coupling.py \
bin/cs_compile.py \
bin/cs_create.py \
bin/cs_exec_environment.py \
bin/cs_gui.py \
bin/cs_info.py \
bin/cs_run.py \
bin/cs_runcase.py \
bin/cs_salome.py \
bin/cs_script.py \
bin/cs_xml_reader.py

autovnvdir = $(pkgpythondir)/autovnv
autovnv_PYTHON = \
bin/autovnv/Command.py \
bin/autovnv/Drawing.py \
bin/autovnv/__init__.py \
bin/autovnv/Parser.py \
bin/autovnv/PlotVTK.py \
bin/autovnv/Study.py \
bin/autovnv/TexMaker.py

# Install user scripts

salomesysconfdir = ${sysconfdir}/salome

if HAVE_SALOME
nodist_salomesysconf_DATA = \
extras/salome/fsi_appli_config.xml
endif

dist_pkgdata_DATA = \
bin/cs_user_scripts.py \
COPYING

dist_pkgdata_SCRIPTS = \
bin/runcase_mpi_rank

# Install thermochemistry and other data files

nobase_dist_pkgdata_DATA = \
data/thch/dp_C3P \
data/thch/dp_C3PSJ \
data/thch/dp_C4P \
data/thch/dp_ELE \
data/thch/dp_FCP.xml \
data/thch/dp_FUE \
data/thch/dp_transfo \
data/thch/JANAF \
data/thch/meteo

# Install batch cards for cluster management

batchdir = ${pkgdatadir}/batch
batch_DATA = \
extras/batch/batch.CCC \
extras/batch/batch.LOADL \
extras/batch/batch.LOADL_BG \
extras/batch/batch.LSF \
extras/batch/batch.OAR \
extras/batch/batch.PBS \
extras/batch/batch.SGE \
extras/batch/batch.SLURM \
extras/batch/batch.SLURM_BGQ

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

# Install SALOME related files

salomedir = ${pkgdatadir}/salome
dist_salome_DATA = \
extras/salome/fsi_yacs_scheme.xml \
extras/salome/fsi.export

# Update gettext dictionnaries

update-po:
	cd po && $(MAKE) $(AM_MAKEFLAGS) update-po

# Generate and install documentation

doc:
	cd docs && $(MAKE) $(AM_MAKEFLAGS) doc

install-doc:
	cd docs && $(MAKE) $(AM_MAKEFLAGS) install-doc

# Generate ChangeLog

changelog:
	start=2018-01-01; \
	if test -d ${top_srcdir}/.svn ; then \
	  today=`date +%F`; \
	  svn2cl --reparagraph --include-rev --include-actions \
	  --authors=${top_srcdir}/extras/changelog/authors.xml \
	  --output=${top_srcdir}/ChangeLog \
	  ${top_srcdir} --revision {"$${today} 23:59"}:{"$${start} 00:00"} ; \
        elif test -d ${top_srcdir}/.git ; then \
	  git --no-pager --git-dir=${top_srcdir}/.git --work-tree=${top_srcdir} log \
	    --since={$${start}} --no-merges --date-order --date=iso \
	    --format="%n%cd %aN <%ce> %n%n%x09* [%H]%n%w(82,8,8)%s" --name-status > ${top_srcdir}/ChangeLog ; \
	else \
	  echo "ChangeLog must be generated from a Subversion or git working copy."  > ${top_srcdir}/ChangeLog ; \
	fi; \
	cat ${top_srcdir}/extras/changelog/ChangeLog.base >> ${top_srcdir}/ChangeLog;

# Update Copyright

update-copyright:
	y=`date +%Y`; yold=$$(($$y-1)); \
	lst=`find ${top_srcdir} -type d \( -name .svn -o -name .git -o -name "build*" \) -prune -o \( -type f -a ! -name Makefile.in \) -print`; \
	for f in $$lst; do \
	  echo $$f && \
	  $(SED) -i "s/Copyright (C) 1998-$${yold} EDF S.A./Copyright (C) 1998-$$y EDF S.A./" $$f; \
	done

dist-hook: changelog
	\cp -f ${top_srcdir}/ChangeLog ${distdir}
	-rm -rf `find $(distdir)/docs -name .svn`
	-rm -rf `find $(distdir)/examples  -name .svn`
	-rm -rf `find $(distdir)/extras  -name .svn`
	-rm -rf `find $(distdir)/gui  -name .svn`
	-rm -rf `find $(distdir)/patches -name .svn`
	-rm -rf `find $(distdir)/salome -name .svn`
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --full > $(distdir)/build-aux/version

all-local: bin/cs_package.py

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

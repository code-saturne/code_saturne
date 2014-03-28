#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009-2011 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Kernel; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

if HAVE_GUI
GUI = gui
endif

SUBDIRS = src libsyrcs po doc $(GUI)

MAINTAINERCLEANFILES = \
$(top_srcdir)/Makefile.in \
$(top_srcdir)/aclocal.m4 \
$(top_srcdir)/configure \
$(top_srcdir)/cs_config-h.in \
$(ac_aux_dir)/compile \
$(ac_aux_dir)/config.guess \
$(ac_aux_dir)/config.sub \
$(ac_aux_dir)/depcomp \
$(ac_aux_dir)/install-sh \
$(ac_aux_dir)/ltmain.sh \
$(ac_aux_dir)/missing \
$(ac_aux_dir)/py-compile

EXTRA_DIST = \
COMPATIBILITY \
QUALITY_ASSURANCE \
po/Makevars.template \
po/Rules-install \
bin/batch.CCC \
bin/batch.LOADL \
bin/batch.LOADL_BG \
bin/batch.LSF \
bin/batch.PBS \
bin/batch.SGE \
bin/batch.SLURM \
bin/batch.SLURM_BGQ \
sbin/backup \
sbin/bootstrap \
sbin/clean \
sbin/rmb \
$(ac_aux_dir)/cs_auto_flags.sh \
$(ac_aux_dir)/cs_batch.m4 \
$(ac_aux_dir)/cs_bft.m4 \
$(ac_aux_dir)/cs_blas.m4 \
$(ac_aux_dir)/cs_config_info.m4 \
$(ac_aux_dir)/cs_fortran.m4 \
$(ac_aux_dir)/cs_fvm.m4 \
$(ac_aux_dir)/cs_libxml2.m4 \
$(ac_aux_dir)/cs_mei.m4 \
$(ac_aux_dir)/cs_modules.m4 \
$(ac_aux_dir)/cs_mpi.m4 \
$(ac_aux_dir)/cs_syrthes.m4

# One adds two whole sub-directories and get rid of the .svn
# directories with a dist-hook. The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += examples extras patches

nodist_include_HEADERS = cs_config.h

DISTCLEANFILES = cs_config.h

ACLOCAL_AMFLAGS = -I config

# Replace paths in Python scripts.
# The following command and targets are adapted from an example in
# autoconf. The autoconf documentation recommends changing scripts
# through Makefile targets rather than through configure, and
# cites its own automate and autoheader as examples.

edit = sed \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@bindir[@]|$(bindir)|g' \
	-e 's|@includedir[@]|$(includedir)|g' \
	-e 's|@libdir[@]|$(libdir)|g' \
	-e 's|@datarootdir[@]|$(datarootdir)|g' \
	-e 's|@datadir[@]|$(datadir)|g' \
	-e 's|@pkgdatadir[@]|$(pkgdatadir)|g' \
	-e 's|@docdir[@]|$(docdir)|g' \
	-e 's|@pdfdir[@]|$(pdfdir)|g' \
	-e 's|@cs_python[@]|$(cs_python)|g' \
	-e 's|@pythondir[@]|$(pythondir)|g' \
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g' \
	-e 's|@mei_prefix[@]|$(mei_prefix)|g' \
	-e 's|@ecs_bindir[@]|$(ecs_bindir)|g' \
	-e 's|@syrthes_prefix[@]|$(syrthes_prefix)|g' \
	-e 's|@mpi_type[@]|$(mpi_type)|g' \
	-e 's|@mpi_bindir[@]|$(mpi_bindir)|g' \
	-e 's|@mpi_libdir[@]|$(mpi_libdir)|g' \
	-e 's|@cs_env_modules_cmd[@]|$(cs_env_modules_cmd)|g' \
	-e 's|@cs_env_modules[@]|$(cs_env_modules)|g' \
	-e 's|@cs_batch_template[@]|$(cs_batch_template)|g' \
	-e 's|@PYTHON_VERSION[@]|$(PYTHON_VERSION)|g' \
	-e 's|@PACKAGE_NAME[@]|$(PACKAGE_NAME)|g' \
	-e 's|@PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g'

edit_scripts = bin/code_saturne bin/cs_config.py bin/runcase.py 

$(edit_scripts): Makefile
	rm -f $@ $@.tmp
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

bin/code_saturne: $(srcdir)/bin/code_saturne.in 
bin/cs_config.py: bin/cs_config.py.in 
bin/runcase.py: $(srcdir)/bin/runcase.py.in 

CLEANFILES = $(edit_scripts)
EXTRA_DIST += bin/code_saturne.in bin/cs_config.py.in bin/runcase.py.in

# Install binary scripts

nodist_bin_SCRIPTS = \
bin/code_saturne

# Install Python library

nodist_pkgpython_PYTHON = \
bin/cs_config.py \
bin/cs_config_build.py

dist_pkgpython_PYTHON = \
bin/cs_case.py \
bin/cs_case_domain.py \
bin/cs_check_consistency.py \
bin/cs_check_mesh.py \
bin/cs_compile.py \
bin/cs_create.py \
bin/cs_exec_environment.py \
bin/cs_info.py \
bin/cs_plot_probes.py

# Install user scripts

nodist_pkgdata_DATA = \
bin/runcase \
bin/runcase.py \
bin/SaturneGUI

dist_pkgdata_DATA = \
bin/runcase.help

dist_pkgdata_SCRIPTS = \
bin/runcase_mpi_env \
bin/runcase_mpi_rank \
bin/runcase_syrthes

# Install thermochemistry and Matisse data files

nobase_dist_pkgdata_DATA = \
data/mati/emm.dat \
data/mati/emm.geom \
data/mati/vault.dat \
data/mati/vault.geom \
data/thch/dp_C3P \
data/thch/dp_C3PSJ \
data/thch/dp_C4P \
data/thch/dp_ELE \
data/thch/dp_FCP \
data/thch/dp_FUE \
data/thch/dp_transfo \
data/thch/JANAF \
data/thch/meteo

# Install man pages

dist_man_MANS = \
doc/manpages/code_saturne.1 \
doc/manpages/cs_solver.1

# Install bash completion script

bashcompletiondir = ${sysconfdir}/bash_completion.d
bashcompletion_DATA = \
extras/bash_completion/code_saturne

# Finish user files installation
# - tag with the package version number
# - remove the test preventing from using a source file
# Done with two sed passes (do not work otherwise)

install-data-hook:
	if test "@cs_batch_template@" != no; then \
	  if test -f $(srcdir)/bin/batch.@cs_batch_template@; then \
	    cat $(srcdir)/bin/batch.@cs_batch_template@ > $(DESTDIR)$(pkgdatadir)/batch.tmp; \
	  elif test -f @cs_batch_template@; then \
	    cat @cs_batch_template@ > $(DESTDIR)$(pkgdatadir)/batch.tmp; \
	  fi; \
	fi; \
	cd $(DESTDIR)$(pkgdatadir) && \
	if test ! -f batch.tmp; then \
	  echo "#" > batch.tmp; \
	fi; \
	for file in runcase runcase.py; do \
	  sed -e "/BATCHTEMPLATE/r batch.tmp" \
	      -e "/BATCHTEMPLATE/d" "$${file}" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done; \
	rm -f batch.tmp
	cd $(DESTDIR)$(pkgdatadir)/users && \
	user_files=`ls */*.f90 2>/dev/null`; \
	stubvers="^\!VERS\n"; \
	version1="\!                      Code_Saturne version $(VERSION)\n"; \
	version2="\!                      --------------------------"; \
	comment1="^\! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START"; \
	comment2="^\! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END"; \
	for file in $${user_files}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version1}$${version2}/" < "$${file}" | \
	  $(SED) -e "/$${comment1}/,/$${comment2}/d" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done

dist-hook:
	-rm -rf `find $(distdir)/examples  -name .svn`
	-rm -rf `find $(distdir)/extras  -name .svn`
	-rm -rf `find $(distdir)/patches -name .svn`

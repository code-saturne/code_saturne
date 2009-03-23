#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009 EDF S.A., France
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

SUBDIRS = src/atmo src/cfbl src/cogz src/cplv src/ctwr src/elec src/fuel \
          src/lagr src/mati src/pprt src/rayt src/base \
          po doc

EXTRA_DIST = \
COMPATIBILITY \
QUALITY_ASSURANCE \
sbin/backup \
sbin/bootstrap \
sbin/clean \
sbin/rmb \
$(ac_aux_dir)/cs_auto_flags.sh \
$(ac_aux_dir)/cs_bft.m4 \
$(ac_aux_dir)/cs_blas.m4 \
$(ac_aux_dir)/cs_fvm.m4 \
$(ac_aux_dir)/cs_libxml2.m4 \
$(ac_aux_dir)/cs_mpi.m4 \
$(ac_aux_dir)/cs_tools.m4

# One adds two whole sub-directories and get rid of the .svn
# directories with a dist-hook. The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += extras patches

nodist_include_HEADERS = cs_config.h

DISTCLEANFILES = cs_config.h

ACLOCAL_AMFLAGS = -I config

# Install binary scripts

dist_bin_SCRIPTS = \
bin/cs_check_consistency \
bin/cs_check_mesh \
bin/cs_compile \
bin/cs_create \
bin/cs_info \
bin/cs_plot_probes

# Install Python library

nodist_pkgpython_PYTHON = \
bin/cs_config.py

# Install user scripts

nodist_pkgdata_DATA = \
bin/runcase \
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
data/thch/JANAF

# Install man pages

dist_man_MANS = \
doc/manpages/cs_check_consistency.1 \
doc/manpages/cs_check_mesh.1 \
doc/manpages/cs_compile.1 \
doc/manpages/cs_create.1 \
doc/manpages/cs_info.1 \
doc/manpages/cs_plot_probes.1 \
doc/manpages/cs_solver.1

# Finish user files installation
# - tag with the package version number
# - remove the test preventing from using a source file
# Done with two sed passes (do not work otherwise)

install-data-hook:
	cd $(DESTDIR)$(pkgdatadir)/users && \
	user_files=`ls */*.f90 2>/dev/null`; \
	stubvers="^\!VERS\n"; \
	version1="\!VERS                  Code_Saturne version $(VERSION)\n"; \
	version2="\!                      --------------------------"; \
	comment1="^\! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT"; \
	comment2="^\! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN"; \
	for file in $${user_files}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version1}$${version2}/" < "$${file}" | \
	  $(SED) -e "/$${comment1}/,/$${comment2}/d" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done

dist-hook:
	-rm -rf `find $(distdir)/extras  -name .svn`
	-rm -rf `find $(distdir)/patches -name .svn`

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
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

AUX_DIST = $(ac_aux_dir)/compile \
           $(ac_aux_dir)/config.guess \
           $(ac_aux_dir)/config.sub \
           $(ac_aux_dir)/depcomp \
           $(ac_aux_dir)/install-sh \
           $(ac_aux_dir)/missing \
           $(ac_aux_dir)/mkinstalldirs \
           $(ac_aux_dir)/ltmain.sh

MAINTAINERCLEANFILES = \
$(top_srcdir)/Makefile.in \
$(top_srcdir)/aclocal.m4 \
$(top_srcdir)/configure \
$(top_srcdir)/cs_config-h.in \
$(AUX_DIST)

SUBDIRS = src/atmo src/cfbl src/cogz src/cplv src/ctwr src/elec src/fuel src/lagr src/mati src/pprt src/rayt src/base \
          po

#EXTRA_SUBDIRS = doc

EXTRA_DIST = $(ac_aux_dir)/config.rpath  \
             $(ac_aux_dir)/cs_auto_flags.sh \
             $(ac_aux_dir)/cs_bft.m4 \
             $(ac_aux_dir)/cs_blas.m4 \
             $(ac_aux_dir)/cs_fvm.m4 \
             $(ac_aux_dir)/cs_libxml2.m4 \
             $(ac_aux_dir)/cs_mpi.m4 \
             $(ac_aux_dir)/cs_tools.m4

include_HEADERS = cs_config.h

DISTCLEANFILES = cs_config.h

ACLOCAL_AMFLAGS = -I config

# Install binary scripts

dist_bin_SCRIPTS = \
bin/check_consistency \
bin/check_mesh \
bin/cs_create \
bin/cs_info \
bin/gracehst \
bin/mpi_rank.sh \
bin/runcase_mpi_env \
bin/runcase_syrthes

# Install user scripts

dist_data_DATA = \
bin/runcase \
bin/runcase.help \
bin/SaturneGUI

# Install thermochemistry and Matisse data files

nobase_dist_data_DATA = \
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

# Tag the user files with the package version number

install-data-hook:
	cd $(DESTDIR)$(datadir)/users && \
	stub="^\!VERS\n"; \
	version1="\!VERS                  Code_Saturne version $(VERSION)\n"; \
	version2="\!                      --------------------------"; \
	for file in */*.f90; do \
	  sed -i -s "N;s/$${stub}/$${version1}$${version2}/" "$${file}"; \
	done

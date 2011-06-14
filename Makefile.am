#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne CFD tool.
#
#   Copyright (C) 2009-2011 EDF S.A., France
#
#   The Code_Saturne CFD tool is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne CFD tool is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne CFD tool; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

if HAVE_GUI
GUI = gui
endif

if HAVE_FRONTEND
PREPRO = preprocessor
endif

if HAVE_INTERNAL_PLE
PLE = libple
PLE_CPPFLAGS_TMP = -I$(includedir)
else
PLE_CPPFLAGS_TMP = $(PLE_CPPFLAGS)
endif

if HAVE_SALOME_KERNEL
SALOME = salome
endif

SUBDIRS = $(PLE) $(PREPRO) src libsyrcs $(SALOME) po doc $(GUI) tests

MAINTAINERCLEANFILES = \
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
build-aux/ylwrap \
m4/gettext.m4 \
m4/iconv.m4 \
m4/lib-ld.m4 \
m4/lib-link.m4 \
m4/lib-prefix.m4 \
m4/libtool.m4 \
m4/lt~obsolete.m4 \
m4/ltoptions.m4 \
m4/ltsugar.m4 \
m4/ltversion.m4 \
m4/nls.m4 \
m4/po.m4 \
m4/progtest.m4

EXTRA_DIST = \
COMPATIBILITY \
QUALITY_ASSURANCE \
po/Makevars.template \
po/Rules-install \
sbin/backup \
sbin/bootstrap \
sbin/clean \
sbin/rmb \
config/cs_auto_flags.sh \
m4/cs_adf.m4 \
m4/cs_blas.m4 \
m4/cs_ccm.m4 \
m4/cs_cgns.m4 \
m4/cs_config_info.m4 \
m4/cs_fortran.m4 \
m4/cs_hdf5.m4 \
m4/cs_libxml2.m4 \
m4/cs_med.m4 \
m4/cs_metis.m4 \
m4/cs_modules.m4 \
m4/cs_mpi.m4 \
m4/cs_ple.m4 \
m4/cs_scotch.m4 \
m4/cs_syrthes.m4 \
m4/cs_zlib.m4

# One adds several sub-directories and get rid of the .svn
# directories with a dist-hook. The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += examples extras patches

nodist_pkginclude_HEADERS = cs_config.h

DISTCLEANFILES = cs_config.h

ACLOCAL_AMFLAGS = -I m4

# Replace paths in Python scripts.
# The following command and targets are adapted from an example in
# autoconf. The autoconf documentation recommends changing scripts
# through Makefile targets rather than through configure, and
# cites its own automate and autoheader as examples.

if HAVE_LINK_CXX
CS_LD=$(CXX)
else
CS_LD=$(CC)
endif

CSCPPFLAGS = $(CGNS_CPPFLAGS) $(MED_CPPFLAGS) $(HDF5_CPPFLAGS) \
$(MPI_CPPFLAGS) $(LIBXML2_CPPFLAGS) $(BLAS_CPPFLAGS)

CSCFLAGS = $(CFLAGS) $(CFLAGS_DBG) $(CFLAGS_OPT)
CSCXXFLAGS = $(CXXFLAGS) $(CXXFLAGS_DBG) $(CXXFLAGS_OPT)
CSFCFLAGS = $(FCFLAGS) $(FCFLAGS_DBG) $(FCFLAGS_OPT)

CSLDFLAGS = $(CGNS_LDFLAGS) $(MED_LDFLAGS) $(HDF5_LDFLAGS) \
$(MPI_LDFLAGS) $(LIBXML2_LDFLAGS) $(BLAS_LDFLAGS)

CSLIBS = $(CGNS_LIBS) $(MED_LIBS) $(HDF5_LIBS) \
$(MPI_LIBS) $(LIBXML2_LIBS) $(BLAS_LIBS)

edit_python = sed \
	-e 's|@PYTHON[@]|$(PYTHON)|g' \
	-e 's|@pkgpythondir[@]|$(pkgpythondir)|g'

edit_config = sed \
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
	-e 's|@datarootdir[@]|$(datarootdir)|g' \
	-e 's|@datadir[@]|$(datadir)|g' \
	-e 's|@pkgdatadir[@]|$(pkgdatadir)|g' \
	-e 's|@docdir[@]|$(docdir)|g' \
	-e 's|@pdfdir[@]|$(pdfdir)|g' \
	-e 's|@sysconfdir[@]|$(sysconfdir)|g' \
	-e 's|@syrthes_prefix[@]|$(syrthes_prefix)|g' \
	-e 's|@mpi_type[@]|$(mpi_type)|g' \
	-e 's|@mpi_bindir[@]|$(mpi_bindir)|g' \
	-e 's|@mpi_libdir[@]|$(mpi_libdir)|g' \
	-e 's|@cs_env_modules[@]|$(cs_env_modules)|g' \
	-e 's|@CC[@]|$(CC)|g' \
	-e 's|@CS_LD[@]|$(CS_LD)|g' \
	-e 's|@CXX[@]|$(CXX)|g' \
	-e 's|@FC[@]|$(FC)|g' \
	-e 's|@CS_CPPFLAGS[@]|$(CSCPPFLAGS)|g' \
	-e 's|@CS_CFLAGS[@]|$(CSCFLAGS)|g' \
	-e 's|@CS_CXXFLAGS[@]|$(CSCXXFLAGS)|g' \
	-e 's|@CS_FCFLAGS[@]|$(CSFCFLAGS)|g' \
	-e 's|@CS_LDFLAGS[@]|$(CSLDFLAGS)|g' \
	-e 's|@CS_LIBS[@]|$(CSLIBS)|g' \
	-e 's|@PLE_CPPFLAGS[@]|$(PLE_CPPFLAGS_TMP)|g' \
	-e 's|@PLE_LDFLAGS[@]|$(PLE_LDFLAGS)|g' \
	-e 's|@PLE_LIBS[@]|$(PLE_LIBS)|g' \
	-e 's|@LDFLAGS[@]|$(LDFLAGS)|g' \
	-e 's|@LIBS[@]|$(LIBS)|g' \
	-e 's|@FCLIBS[@]|$(FCLIBS)|g' \
	-e 's|@FCMODINCLUDE[@]|$(FCMODINCLUDE)|g' \
	-e 's|@LDRPATH[@]|$(LDRPATH)|g' \
	-e 's|@PYTHON[@]|$(PYTHON)|g' \
	-e 's|@PYUIC4[@]|$(PYUIC4)|g' \
	-e 's|@PYRCC4[@]|$(PYRCC4)|g' \
	-e 's|@USE_NLS[@]|$(USE_NLS)|g' \
	-e 's|@PACKAGE_NAME[@]|$(PACKAGE_NAME)|g' \
	-e 's|@PACKAGE_TARNAME[@]|$(PACKAGE_TARNAME)|g' \
	-e 's|@PACKAGE_VERSION[@]|$(PACKAGE_VERSION)|g' \
	-e 's|@PACKAGE_STRING[@]|$(PACKAGE_STRING)|g' \
	-e 's|@PACKAGE_BUGREPORT[@]|$(PACKAGE_BUGREPORT)|g' \
	-e 's|@PACKAGE_URL[@]|$(PACKAGE_URL)|g' \
	-e 's|@KERNEL_ROOT_DIR[@]|$(KERNEL_ROOT_DIR)|g' \
	-e 's|@ROOT_SALOME[@]|$(ROOT_SALOME)|g' \
	-e 's|@SALOME_VERSION[@]|$(SALOME_VERSION)|g'

bin/code_saturne bin/runcase_coupling: Makefile
	rm -f $@ $@.tmp
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_python) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

bin/cs_package.py bin/runcase_aster: Makefile
	rm -f $@ $@.tmp
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_config) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

bin/code_saturne: $(srcdir)/bin/code_saturne.in 
bin/cs_package.py: $(srcdir)/bin/cs_package.py.in 
bin/runcase_aster: $(srcdir)/bin/runcase_aster.in 
bin/runcase_coupling: $(srcdir)/bin/runcase_coupling.in 

CLEANFILES = $(edit_scripts)
EXTRA_DIST += bin/code_saturne.in bin/cs_package.py.in bin/runcase_aster.in bin/runcase_coupling.in

# Install binary scripts

if HAVE_SALOME_KERNEL
CS_SALOME_WRAPPER = bin/cs_salome_wrapper
endif

nodist_bin_SCRIPTS = \
bin/code_saturne $(CS_SALOME_WRAPPER)

# Install Python library

nodist_pkgpython_PYTHON = \
bin/cs_config.py \
bin/cs_package.py

dist_pkgpython_PYTHON = \
bin/__init__.py \
bin/cs_autovnv.py \
bin/cs_case.py \
bin/cs_case_domain.py \
bin/cs_case_coupling.py \
bin/cs_compile.py \
bin/cs_create.py \
bin/cs_exec_environment.py \
bin/cs_gui.py \
bin/cs_info.py \
bin/cs_run.py \
bin/cs_salome.py \
bin/cs_script.py \
bin/cs_xml_reader.py

autovnvdir = $(pkgpythondir)/autovnv
autovnv_PYTHON = \
bin/autovnv/Command.py \
bin/autovnv/Drawing.py \
bin/autovnv/__init__.py \
bin/autovnv/Parser.py \
bin/autovnv/Study.py \
bin/autovnv/TexMaker.py

# Install user scripts

nodist_pkgdata_DATA = \
bin/cs_user_scripts.py \
bin/runcase_aster \
bin/runcase_coupling \
bin/SaturneGUI

dist_pkgdata_SCRIPTS = \
bin/runcase_mpi_rank \
bin/runcase_syrthes

# Install thermochemistry and other data files

nobase_dist_pkgdata_DATA = \
data/thch/dp_C3P \
data/thch/dp_C3PSJ \
data/thch/dp_C4P \
data/thch/dp_ELE \
data/thch/dp_FCP \
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
extras/batch/batch.PBS \
extras/batch/batch.SGE \
extras/batch/batch.SLURM

# Install man pages

dist_man_MANS = \
doc/manpages/code_saturne.1

# Install bash completion script

bashcompletiondir = ${sysconfdir}/bash_completion.d
bashcompletion_DATA = \
extras/bash_completion/code_saturne

# Install global option file

sysconf_DATA = \
extras/code_saturne.cfg

# Install SALOME related files

salomedir = ${pkgdatadir}/salome
salome_DATA = \
extras/salome/fsi_yacs_scheme.xml

nodist_salome_DATA = \
extras/salome/fsi_appli_config.xml

# Finish user files installation
# - tag with the package version number
# - remove the test preventing from using a source file
# Done with two sed passes (do not work otherwise)

install-data-hook:
	cd $(DESTDIR)$(pkgdatadir)/users && \
	f90_user_files=`ls */*.f90 2>/dev/null`; \
	stubvers="^\!VERS\n"; \
	version1="\!                      Code_Saturne version $(VERSION)\n"; \
	version2="\!                      --------------------------"; \
	comment1="^\! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START"; \
	comment2="^\! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END"; \
	for file in $${f90_user_files}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version1}$${version2}/" < "$${file}" | \
	  $(SED) -e "/$${comment1}/,/$${comment2}/d" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done ;\
	c_user_files=`ls */*.c 2>/dev/null`; \
	stubvers="\/\* VERS \*\/"; \
	version="\/\* Code_Saturne version $(VERSION) \*\/"; \
	comment="\/\* REMOVE_LINE_FOR_USE_OF_SUBROUTINE \*\/"; \
	for file in $${c_user_files}; do \
	  $(SED) -e "N;s/$${stubvers}/$${version}/" < "$${file}" | \
	  $(SED) -e "/$${comment}/d" > "$${file}.tmp" && \
	  mv "$${file}.tmp" "$${file}"; \
	done

dist-hook:
	-rm -rf `find $(distdir)/gui -name .svn`
	-rm -rf `find $(distdir)/examples  -name .svn`
	-rm -rf `find $(distdir)/extras  -name .svn`
	-rm -rf `find $(distdir)/patches -name .svn`

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

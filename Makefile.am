## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
DOCS = docs
endif

if HAVE_BACKEND
TSTS = tests
endif

if HAVE_INTERNAL_PLE
PLE = libple
endif

if HAVE_SALOME
SALOME = salome
endif

SUBDIRS = . $(PLE) $(PREPRO) src $(SALOME) po $(DOCS) $(GUI) $(TSTS) $(NEPTUNE_CFD)

AUX_DIST = \
build-aux/catalyst \
build-aux/catalyst-5.6 \
build-aux/cs_compile_build.py \
build-aux/cs_config_test.py \
build-aux/cs_version.py \
build-aux/list_modules.sh \
build-aux/petsc-petscvariables.makefile \
build-aux/petsc-variables.makefile \
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
NEWS.md \
README.md \
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
m4/cs_coolprop.m4 \
m4/cs_docs.m4 \
m4/cs_eos.m4 \
m4/cs_fortran.m4 \
m4/cs_freesteam.m4 \
m4/cs_hdf5.m4 \
m4/cs_med.m4 \
m4/cs_medcoupling.m4 \
m4/cs_melissa.m4 \
m4/cs_metis.m4 \
m4/cs_modules.m4 \
m4/cs_mpi.m4 \
m4/cs_omniorb.m4 \
m4/cs_petsc.m4 \
m4/cs_ple.m4 \
m4/cs_salome.m4 \
m4/cs_scotch.m4 \
m4/cs_zlib.m4

# One adds several sub-directories with a dist-hook.
# The reason for this trick is that
# the number of extra files or patches can grow and change a lot.
EXTRA_DIST += examples extras patches

nodist_pkginclude_HEADERS = \
cs_config.h \
cs_version.h

DISTCLEANFILES = \
bin/code_saturne \
bin/cs_package.py \
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
	-e 's|@SALOME_KERNEL[@]|$(SALOME_KERNEL)|g' \
	-e 's|@SALOMEENVCMD[@]|$(SALOMEENVCMD)|g'

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
	${top_srcdir}/build-aux/cs_version.py --replace $@ $@

if HAVE_SHELL_LAUNCHER
bin/code_saturne.py: Makefile $(srcdir)/bin/code_saturne.py.in
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_python) $${srcdir}$@.in >$@.tmp;
	mv $@.tmp $@
endif

bin/cs_package.py: Makefile $(srcdir)/bin/cs_package.py.in
	rm -f $@ $@.tmp
	$(MKDIR_P) bin
	srcdir=''; \
	test -f ./$@.in || srcdir=$(srcdir)/; \
	$(edit_package) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace $@ $@

bin/code_saturne: $(srcdir)/bin/code_saturne.$(LAUNCHER_EXT).in

cs_version.h: Makefile $(top_srcdir)/src/base/cs_version.h.in
	$(PYTHON) ${top_srcdir}/build-aux/cs_version.py --replace ${top_srcdir}/src/base/cs_version.h.in cs_version.h

BUILT_SOURCES = \
bin/cs_package.py \
cs_version.h

CLEANFILES = $(edit_scripts)
EXTRA_DIST += bin/code_saturne.py.in bin/cs_package.py.in \
	      bin/code_saturne.sh.in \
              data/thch/dp_radiat_MFS.tar.xz

# Install binary scripts

nodist_bin_SCRIPTS = \
bin/code_saturne

if HAVE_SHELL_LAUNCHER
nodist_bin_SCRIPTS += bin/code_saturne.py
endif

# Install Python library

nodist_pkgpython_PYTHON = \
bin/cs_config.py \
bin/cs_package.py

# Install Python library

dist_pkgpython_PYTHON = \
bin/__init__.py \
bin/cs_studymanager.py \
bin/cs_batch.py \
bin/cs_bdiff.py \
bin/cs_bdump.py \
bin/cs_case.py \
bin/cs_case_domain.py \
bin/cs_case_coupling.py \
bin/cs_compile.py \
bin/cs_control.py \
bin/cs_create.py \
bin/cs_debug_wrapper.py \
bin/cs_exec_environment.py \
bin/cs_studymanager_gui.py \
bin/cs_trackcvg.py \
bin/cs_gui.py \
bin/cs_info.py \
bin/cs_run.py \
bin/cs_runcase.py \
bin/cs_salome.py \
bin/cs_script.py \
bin/cs_submit.py \
bin/cs_mei_to_c.py \
bin/cs_update.py \
bin/cs_xml_reader.py

studymanagerdir = $(pkgpythondir)/studymanager
studymanager_PYTHON = \
bin/studymanager/cs_studymanager_run.py \
bin/studymanager/cs_studymanager_drawing.py \
bin/studymanager/__init__.py \
bin/studymanager/cs_studymanager_parser.py \
bin/studymanager/cs_studymanager_study.py \
bin/studymanager/cs_studymanager_texmaker.py \
bin/studymanager/cs_studymanager_pathes_model.py \
bin/studymanager/cs_studymanager_xml_init.py

modeldir = $(pkgpythondir)/model
model_PYTHON = \
bin/model/Common.py \
bin/model/AtmosphericFlowsModel.py \
bin/model/BalanceModelNeptune.py \
bin/model/BalanceModel.py \
bin/model/BatchRunningModel.py \
bin/model/BodyForcesModel.py \
bin/model/BoundaryConditionsModelNeptune.py \
bin/model/BoundaryNeptune.py \
bin/model/Boundary.py \
bin/model/CathareCouplingModel.py \
bin/model/CoalCombustionModel.py \
bin/model/CompressibleModel.py \
bin/model/ConjugateHeatTransferModel.py \
bin/model/CoriolisSourceTermsModel.py \
bin/model/DefineUserScalarsModel.py \
bin/model/DropletCondensationEvaporationModel.py \
bin/model/ElectricalModel.py \
bin/model/FansModel.py \
bin/model/FluidCharacteristicsModel.py \
bin/model/FluidStructureInteractionModel.py \
bin/model/GasCombustionModel.py \
bin/model/GlobalNumericalParametersModel.py \
bin/model/GroundwaterLawModel.py \
bin/model/GroundwaterModel.py \
bin/model/HgnModel.py \
bin/model/HeadLossesModel.py \
bin/model/IdentityAndPathesModel.py \
bin/model/ImmersedBoundariesModel.py \
bin/model/InitializationModel.py \
bin/model/__init__.py \
bin/model/InterfacialAreaModel.py \
bin/model/InterfacialEnthalpyModel.py \
bin/model/InterfacialForcesModel.py \
bin/model/LagrangianBoundariesModel.py \
bin/model/LagrangianModel.py \
bin/model/LagrangianOutputModel.py \
bin/model/LagrangianStatisticsModel.py \
bin/model/LocalizationModelNeptune.py \
bin/model/LocalizationModel.py \
bin/model/MainFieldsInitializationModel.py \
bin/model/MainFieldsSourceTermsModel.py \
bin/model/MainFieldsModel.py \
bin/model/MobileMeshModel.py \
bin/model/NonCondensableModel.py \
bin/model/NotebookModel.py \
bin/model/NucleateBoilingModel.py \
bin/model/NumericalParamEquationModelNeptune.py \
bin/model/NumericalParamEquationModel.py \
bin/model/NumericalParamGlobalModel.py \
bin/model/OpenTurnsModel.py \
bin/model/OutputControlModel.py \
bin/model/OutputFieldsModel.py \
bin/model/OutputSurfacicFieldsModel.py \
bin/model/OutputSurfacicVariablesModel.py \
bin/model/OutputVolumicVariablesModel.py \
bin/model/PerformanceTuningModel.py \
bin/model/PorosityModel.py \
bin/model/ProfilesModel.py \
bin/model/ScriptRunningModel.py \
bin/model/SolidModel.py \
bin/model/SolutionDomainModel.py \
bin/model/SourceTermsModel.py \
bin/model/SpeciesModel.py \
bin/model/StartRestartModel.py \
bin/model/ThermalRadiationModel.py \
bin/model/ThermalScalarModel.py \
bin/model/ThermodynamicsModel.py \
bin/model/TimeAveragesModel.py \
bin/model/TimeStepModelNeptune.py \
bin/model/TimeStepModel.py \
bin/model/TurboMachineryModel.py \
bin/model/TurbulenceModel.py \
bin/model/TurbulenceNeptuneModel.py \
bin/model/UsersControlModel.py \
bin/model/XMLengine.py \
bin/model/XMLinitializeNeptune.py \
bin/model/XMLinitialize.py \
bin/model/XMLmodel.py \
bin/model/XMLvariables.py

# Install user scripts

salomesysconfdir = ${sysconfdir}/salome

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
data/thch/dp_FUE \
data/thch/dp_radiat_ADF50 \
data/thch/dp_radiat_ADF8 \
data/thch/dp_radiat_Planck_CO2 \
data/thch/dp_radiat_Planck_H2O \
data/thch/dp_radiat_wave \
data/thch/dp_transformers \
data/thch/JANAF \
data/thch/meteo \
data/setup.xml

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

# Update gettext dictionnaries

update-po:
	cd po && $(MAKE) $(AM_MAKEFLAGS) update-po

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
	test -z "$(pkgdatadir)/data/thch" || $(MKDIR_P) "$(DESTDIR)$(pkgdatadir)/data/thch"
	tar xvJf "$(top_srcdir)/data/thch/dp_radiat_MFS.tar.xz" -C "$(DESTDIR)$(pkgdatadir)/data/thch"

uninstall-local:
	rm -rf "$(DESTDIR)$(pkgdatadir)/data/thch/dp_radiat_MFS"
	@$(NORMAL_UNINSTALL)

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

SUBDIRS = . bft mei

if HAVE_BACKEND
SUBDIRS += \
fvm gui \
base cdo pprt alge mesh turb darc \
atmo cfbl cogz comb ctwr elec lagr rayt \
user user_examples meg apps
endif

# General flags

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS = -no-undefined

# Library source files

noinst_LTLIBRARIES = libcsmodules.la

if HAVE_BACKEND

libcsmodules_la_SOURCES = \
base/paramx.f90 \
base/mesh.f90 \
base/field.f90 \
base/field_operator.f90 \
base/alaste.f90 \
base/albase.f90 \
base/alstru.f90 \
base/cavitation.f90 \
base/cdomod.f90 \
base/cplsat.f90 \
base/cstphy.f90 \
base/cstnum.f90 \
base/dimens.f90 \
base/entsor.f90 \
base/ihmpre.f90 \
base/cs_nz_condensation.f90 \
base/cs_nz_tagmr.f90 \
base/cs_tagmr.f90 \
base/cs_tagms.f90 \
base/numvar.f90 \
base/optcal.f90 \
base/parall.f90 \
base/period.f90 \
base/pointe.f90 \
base/cs_c_bindings.f90 \
base/cs_f_interfaces.f90 \
base/post.f90 \
base/rotation.f90 \
atmo/atincl.f90 \
atmo/atsoil.f90 \
atmo/atchem.f90 \
atmo/atimbr.f90 \
atmo/siream.f90 \
atmo/spefun.f90 \
cfbl/cfpoin.f90 \
cfbl/cs_cf_bindings.f90 \
comb/cs_coal_incl.f90 \
comb/cs_fuel_incl.f90 \
comb/cpincl.f90 \
ctwr/ctincl.f90 \
lagr/lagran.f90 \
pprt/ppcpfu.f90 \
pprt/ppincl.f90 \
pprt/ppppar.f90 \
pprt/ppthch.f90 \
cogz/coincl.f90 \
rayt/radiat.f90 \
base/vof.f90 \
turb/vorinc.f90 \
base/ptrglo.f90 \
darc/darcy_module.f90 \
base/turbomachinery.f90

# Module dependencies (handled manually, as Automake does
# not yet handle Fortran module dependencies)

base/albase.lo: base/cplsat.lo base/optcal.lo base/cs_c_bindings.lo
base/alstru.lo: base/paramx.lo
base/cavitation.lo: base/cstphy.lo base/field.lo base/numvar.lo base/optcal.lo\
 base/pointe.lo base/mesh.lo base/paramx.lo base/entsor.lo base/parall.lo \
 base/rotation.lo base/turbomachinery.lo base/vof.lo
base/vof.lo: base/cstphy.lo base/field.lo base/numvar.lo base/optcal.lo \
 base/pointe.lo base/mesh.lo base/paramx.lo base/entsor.lo base/parall.lo \
 base/rotation.lo base/turbomachinery.lo
base/cplsat.lo: base/paramx.lo
base/cstphy.lo: base/paramx.lo
base/cs_c_bindings.lo: base/paramx.lo base/entsor.lo base/mesh.lo base/period.lo \
 base/field.lo base/field_operator.lo
base/cs_f_interfaces.lo: base/cs_c_bindings.lo base/paramx.lo base/dimens.lo \
 base/mesh.lo base/optcal.lo base/field.lo lagr/lagran.lo
base/cs_nz_condensation.lo: base/pointe.lo
base/cs_nz_tagmr.lo: base/optcal.lo base/pointe.lo base/parall.lo \
 base/cs_nz_condensation.lo base/cs_tagmr.lo
base/cs_tagmr.lo: base/pointe.lo
base/cs_tagms.lo: base/mesh.lo
base/entsor.lo: base/paramx.lo base/field.lo
base/field_operator.lo: base/paramx.lo base/dimens.lo base/mesh.lo \
 base/optcal.lo base/field.lo
base/numvar.lo: base/paramx.lo
base/optcal.lo: base/paramx.lo base/field.lo base/numvar.lo base/cstphy.lo
base/pointe.lo: base/paramx.lo base/numvar.lo base/mesh.lo base/parall.lo \
 base/period.lo base/optcal.lo base/entsor.lo lagr/lagran.lo rayt/radiat.lo \
 pprt/ppincl.lo base/ihmpre.lo base/albase.lo base/cstnum.lo base/cs_c_bindings.lo
base/ptrglo.lo: base/paramx.lo base/dimens.lo base/mesh.lo base/numvar.lo \
 base/optcal.lo base/parall.lo base/period.lo base/mesh.lo base/field.lo
base/rotation.lo: base/paramx.lo
base/turbomachinery.lo:  base/mesh.lo base/cstphy.lo base/cplsat.lo
atmo/atincl.lo: pprt/ppppar.lo pprt/ppincl.lo base/mesh.lo atmo/atsoil.lo \
 base/cs_c_bindings.lo
atmo/atchem.lo: pprt/ppppar.lo base/mesh.lo base/entsor.lo base/numvar.lo
atmo/atimbr.lo: base/cstphy.lo base/entsor.lo atmo/atincl.lo
atmo/siream.lo: pprt/ppppar.lo
atmo/spefun.lo: base/entsor.lo
cogz/coincl.lo: base/paramx.lo pprt/ppppar.lo pprt/ppincl.lo
comb/cpincl.lo: pprt/ppppar.lo pprt/ppthch.lo
comb/cs_coal_incl.lo: pprt/ppppar.lo pprt/ppthch.lo
comb/cs_fuel_incl.lo: pprt/ppppar.lo pprt/ppthch.lo
ctwr/ctincl.lo: base/cstphy.lo
lagr/lagran.lo: base/cstnum.lo base/optcal.lo pprt/ppincl.lo rayt/radiat.lo \
 comb/cpincl.lo pprt/ppppar.lo base/cstphy.lo
pprt/ppcpfu.lo: pprt/ppthch.lo pprt/ppppar.lo
pprt/ppincl.lo: pprt/ppppar.lo pprt/ppthch.lo
pprt/ppthch.lo: base/cstphy.lo
rayt/radiat.lo: pprt/ppppar.lo pprt/ppincl.lo base/optcal.lo pprt/ppcpfu.lo base/numvar.lo
darc/darcy_module.lo: base/paramx.lo
cfbl/cfpoin.lo: pprt/ppincl.lo

# Additional local targets

clean-local:
	-rm -f *__genmod.f90 *.mod

install-exec-local: $(libcsmodules_la_SOURCES)
	@$(NORMAL_INSTALL)
	test -z "$(pkgincludedir)" || $(MKDIR_P) "$(DESTDIR)$(pkgincludedir)"
	@list='$(libcsmodules_la_SOURCES)'; test -n "$(pkgincludedir)" || list=; \
	modlist=`for p in $$list; do echo $$p; done | $(SED) -e 's|.f90|.$(FCMODEXT)|g'`; \
	files=`for p in $$modlist; do echo $$p; done | $(SED) -e 's|^.*/||'`; \
	test -n "$$files" || exit 0; \
	for p in $$files; do \
	  if test ! -f $$p; then \
	    p=`echo $$p | sed -e 's/.$(FCMODEXT)//' -e 's/\(^.*\)/\U\1/'`.$(FCMODEXT); \
	  fi; \
	  echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(pkgincludedir)'"; \
	  $(INSTALL_HEADER) $$p "$(DESTDIR)$(pkgincludedir)" || exit $$?; \
	done

uninstall-local:
	@$(NORMAL_UNINSTALL)
	@list='$(libcsmodules_la_SOURCES)'; test -n "$(pkgincludedir)" || list=; \
	modlist=`for p in $$list; do echo $$p; done | $(SED) -e 's|.f90|.$(FCMODEXT)|g'`; \
	lcfiles=`for p in $$modlist; do echo $$p; done | $(SED) -e 's|^.*/||'`; \
	test -d "$(DESTDIR)$(pkgincludedir)" || exit 0; \
	test -n "$$lcfiles" || exit 0; \
	echo " ( cd '$(DESTDIR)$(pkgincludedir)' )"; \
	files=""; \
	for p in $$lcfiles; do \
	  if test ! -f $$p; then \
	    p=`echo $$p | sed -e 's/.$(FCMODEXT)//' -e 's/\(^.*\)/\U\1/'`.$(FCMODEXT); \
	  fi; \
	  files="$$files $$p"; \
	done; \
	echo " ( cd '$(DESTDIR)$(pkgincludedir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(pkgincludedir)" && rm -f $$files
	rmdir $(pkgincludedir) || exit 0

endif

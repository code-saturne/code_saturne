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

SUBDIRS = . bft

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
$(FCFLAGS_DBG) $(FCFLAGS_OPT) $(FCFLAGS_SHARED)

AM_LDFLAGS = -no-undefined

# Library source files

noinst_LIBRARIES = libcsmodules.a

if HAVE_BACKEND

libcsmodules_a_SOURCES = \
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
atmo/sshaerosol.f90 \
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
base/ptrglo.f90 \
darc/darcy_module.f90 \
base/turbomachinery.f90

# Module dependencies (handled manually, as Automake does
# not yet handle Fortran module dependencies)

base/albase.$(OBJEXT): base/cplsat.$(OBJEXT) base/optcal.$(OBJEXT) \
 base/cs_c_bindings.$(OBJEXT)

base/alstru.$(OBJEXT): base/paramx.$(OBJEXT)

base/cavitation.$(OBJEXT): base/cstphy.$(OBJEXT) base/field.$(OBJEXT) \
 base/numvar.$(OBJEXT) base/optcal.$(OBJEXT) base/pointe.$(OBJEXT) \
 base/mesh.$(OBJEXT) base/paramx.$(OBJEXT) base/entsor.$(OBJEXT) \
 base/parall.$(OBJEXT) base/rotation.$(OBJEXT) \
 base/turbomachinery.$(OBJEXT) base/vof.$(OBJEXT)

base/vof.$(OBJEXT): base/cstphy.$(OBJEXT) base/field.$(OBJEXT) \
 base/numvar.$(OBJEXT) base/optcal.$(OBJEXT) base/pointe.$(OBJEXT) \
 base/mesh.$(OBJEXT) base/paramx.$(OBJEXT) base/entsor.$(OBJEXT) \
 base/parall.$(OBJEXT) base/rotation.$(OBJEXT) base/turbomachinery.$(OBJEXT)

base/cplsat.$(OBJEXT): base/paramx.$(OBJEXT)

base/cstphy.$(OBJEXT): base/paramx.$(OBJEXT)

base/cs_c_bindings.$(OBJEXT): base/paramx.$(OBJEXT) base/entsor.$(OBJEXT) \
 base/mesh.$(OBJEXT) base/period.$(OBJEXT) base/field.$(OBJEXT) \
 base/cplsat.$(OBJEXT)

base/cs_f_interfaces.$(OBJEXT): base/cs_c_bindings.$(OBJEXT) \
 base/paramx.$(OBJEXT) base/dimens.$(OBJEXT) base/mesh.$(OBJEXT) \
 base/optcal.$(OBJEXT) base/field.$(OBJEXT) lagr/lagran.$(OBJEXT)

base/cs_nz_tagmr.$(OBJEXT): base/optcal.$(OBJEXT) base/pointe.$(OBJEXT) \
 base/parall.$(OBJEXT) base/cs_nz_condensation.$(OBJEXT)

base/cs_tagms.$(OBJEXT): base/mesh.$(OBJEXT)

base/entsor.$(OBJEXT): base/paramx.$(OBJEXT) base/field.$(OBJEXT)

base/field_operator.$(OBJEXT): base/paramx.$(OBJEXT) base/dimens.$(OBJEXT) \
 base/mesh.$(OBJEXT) base/optcal.$(OBJEXT) base/field.$(OBJEXT)

base/numvar.$(OBJEXT): base/paramx.$(OBJEXT)

base/optcal.$(OBJEXT): base/paramx.$(OBJEXT) base/field.$(OBJEXT) \
 base/numvar.$(OBJEXT) base/cstphy.$(OBJEXT) base/cs_c_bindings.$(OBJEXT) \
 base/cdomod.$(OBJEXT)

base/pointe.$(OBJEXT): base/paramx.$(OBJEXT) base/numvar.$(OBJEXT) \
 base/mesh.$(OBJEXT) base/parall.$(OBJEXT) base/period.$(OBJEXT) \
 base/optcal.$(OBJEXT) base/entsor.$(OBJEXT) lagr/lagran.$(OBJEXT) \
 rayt/radiat.$(OBJEXT) pprt/ppincl.$(OBJEXT) base/albase.$(OBJEXT) \
 base/cstnum.$(OBJEXT) base/cs_c_bindings.$(OBJEXT)

base/ptrglo.$(OBJEXT): base/paramx.$(OBJEXT) base/dimens.$(OBJEXT) \
 base/mesh.$(OBJEXT) base/numvar.$(OBJEXT) base/optcal.$(OBJEXT) \
 base/parall.$(OBJEXT) base/period.$(OBJEXT) base/mesh.$(OBJEXT) \
 base/field.$(OBJEXT)

base/rotation.$(OBJEXT): base/paramx.$(OBJEXT)

base/turbomachinery.$(OBJEXT): base/mesh.$(OBJEXT) base/cstphy.$(OBJEXT) \
 base/cplsat.$(OBJEXT)

atmo/atincl.$(OBJEXT): pprt/ppppar.$(OBJEXT) pprt/ppincl.$(OBJEXT) \
 base/mesh.$(OBJEXT) atmo/atsoil.$(OBJEXT) base/cs_c_bindings.$(OBJEXT) \
 atmo/atchem.$(OBJEXT) atmo/sshaerosol.$(OBJEXT)

atmo/atchem.$(OBJEXT): pprt/ppppar.$(OBJEXT) base/mesh.$(OBJEXT) \
 base/entsor.$(OBJEXT) base/numvar.$(OBJEXT) \
 base/cs_c_bindings.$(OBJEXT)

atmo/atimbr.$(OBJEXT): base/cstphy.$(OBJEXT) base/entsor.$(OBJEXT) \
 atmo/atincl.$(OBJEXT) ctwr/ctincl.$(OBJEXT)

atmo/sshaerosol.$(OBJEXT): pprt/ppppar.$(OBJEXT) base/cs_c_bindings.$(OBJEXT) \
 atmo/atchem.$(OBJEXT) base/mesh.$(OBJEXT)

atmo/spefun.$(OBJEXT): base/entsor.$(OBJEXT)

cogz/coincl.$(OBJEXT): base/paramx.$(OBJEXT) pprt/ppppar.$(OBJEXT) \
 pprt/ppincl.$(OBJEXT) rayt/radiat.$(OBJEXT)

comb/cpincl.$(OBJEXT): pprt/ppppar.$(OBJEXT) pprt/ppthch.$(OBJEXT)

comb/cs_coal_incl.$(OBJEXT): pprt/ppppar.$(OBJEXT) pprt/ppthch.$(OBJEXT)

comb/cs_fuel_incl.$(OBJEXT): pprt/ppppar.$(OBJEXT) pprt/ppthch.$(OBJEXT)

ctwr/ctincl.$(OBJEXT): base/cstphy.$(OBJEXT)

lagr/lagran.$(OBJEXT): base/cstnum.$(OBJEXT) base/optcal.$(OBJEXT) \
 pprt/ppincl.$(OBJEXT) rayt/radiat.$(OBJEXT) comb/cpincl.$(OBJEXT) \
 pprt/ppppar.$(OBJEXT) base/cstphy.$(OBJEXT)

pprt/ppcpfu.$(OBJEXT): pprt/ppthch.$(OBJEXT) pprt/ppppar.$(OBJEXT)

pprt/ppincl.$(OBJEXT): pprt/ppppar.$(OBJEXT) pprt/ppthch.$(OBJEXT)

pprt/ppthch.$(OBJEXT): base/cstphy.$(OBJEXT)

rayt/radiat.$(OBJEXT): pprt/ppppar.$(OBJEXT) pprt/ppincl.$(OBJEXT) \
 base/optcal.$(OBJEXT) pprt/ppcpfu.$(OBJEXT) base/numvar.$(OBJEXT)

darc/darcy_module.$(OBJEXT): base/paramx.$(OBJEXT)

cfbl/cfpoin.$(OBJEXT): pprt/ppincl.$(OBJEXT)

# Additional local targets

clean-local:
	-rm -f *__genmod.f90 *.mod

install-exec-local: $(libcsmodules_a_SOURCES)
	@$(NORMAL_INSTALL)
	test -z "$(pkgincludedir)" || $(MKDIR_P) "$(DESTDIR)$(pkgincludedir)"
	@list='$(libcsmodules_a_SOURCES)'; test -n "$(pkgincludedir)" || list=; \
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
	@list='$(libcsmodules_a_SOURCES)'; test -n "$(pkgincludedir)" || list=; \
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

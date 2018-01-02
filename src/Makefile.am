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

SUBDIRS = \
. bft fvm mei gui \
base cdo pprt alge mesh turb darc \
atmo cfbl cogz comb ctwr elec lagr rayt \
user user_examples apps

# General flags

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS = -no-undefined

# Library source files

noinst_LTLIBRARIES = libcsmodules.la

libcsmodules_la_SOURCES = \
$(top_srcdir)/src/base/paramx.f90 \
$(top_srcdir)/src/base/mesh.f90 \
$(top_srcdir)/src/base/field.f90 \
$(top_srcdir)/src/base/field_operator.f90 \
$(top_srcdir)/src/base/alaste.f90 \
$(top_srcdir)/src/base/albase.f90 \
$(top_srcdir)/src/base/alstru.f90 \
$(top_srcdir)/src/base/cavitation.f90 \
$(top_srcdir)/src/base/cplsat.f90 \
$(top_srcdir)/src/base/cstphy.f90 \
$(top_srcdir)/src/base/cstnum.f90 \
$(top_srcdir)/src/base/dimens.f90 \
$(top_srcdir)/src/base/entsor.f90 \
$(top_srcdir)/src/base/ihmpre.f90 \
$(top_srcdir)/src/base/cs_nz_condensation.f90 \
$(top_srcdir)/src/base/cs_nz_tagmr.f90 \
$(top_srcdir)/src/base/cs_tagmr.f90 \
$(top_srcdir)/src/base/cs_tagms.f90 \
$(top_srcdir)/src/base/numvar.f90 \
$(top_srcdir)/src/base/optcal.f90 \
$(top_srcdir)/src/base/parall.f90 \
$(top_srcdir)/src/base/period.f90 \
$(top_srcdir)/src/base/pointe.f90 \
$(top_srcdir)/src/base/cs_c_bindings.f90 \
$(top_srcdir)/src/base/cs_f_interfaces.f90 \
$(top_srcdir)/src/base/post.f90 \
$(top_srcdir)/src/base/rotation.f90 \
$(top_srcdir)/src/atmo/atincl.f90 \
$(top_srcdir)/src/atmo/atsoil.f90 \
$(top_srcdir)/src/atmo/atchem.f90 \
$(top_srcdir)/src/atmo/atimbr.f90 \
$(top_srcdir)/src/atmo/siream.f90 \
$(top_srcdir)/src/atmo/spefun.f90 \
$(top_srcdir)/src/cfbl/cfpoin.f90 \
$(top_srcdir)/src/cfbl/cs_cf_bindings.f90 \
$(top_srcdir)/src/comb/cs_coal_incl.f90 \
$(top_srcdir)/src/comb/cs_fuel_incl.f90 \
$(top_srcdir)/src/comb/cpincl.f90 \
$(top_srcdir)/src/ctwr/ctincl.f90 \
$(top_srcdir)/src/lagr/lagran.f90 \
$(top_srcdir)/src/pprt/ppcpfu.f90 \
$(top_srcdir)/src/pprt/ppincl.f90 \
$(top_srcdir)/src/pprt/ppppar.f90 \
$(top_srcdir)/src/pprt/ppthch.f90 \
$(top_srcdir)/src/cogz/coincl.f90 \
$(top_srcdir)/src/rayt/radiat.f90 \
$(top_srcdir)/src/base/vof.f90 \
$(top_srcdir)/src/turb/vorinc.f90 \
$(top_srcdir)/src/base/ptrglo.f90 \
$(top_srcdir)/src/darc/darcy_module.f90 \
$(top_srcdir)/src/base/turbomachinery.f90

# Module dependencies (handled manually, as Automake does
# not yet handle Fortran module dependencies)

$(top_srcdir)/src/base/albase.f90: cplsat.lo optcal.lo
$(top_srcdir)/src/base/alstru.f90: paramx.lo
$(top_srcdir)/src/base/cavitation.f90: cstphy.lo field.lo numvar.lo optcal.lo pointe.lo mesh.lo paramx.lo entsor.lo parall.lo rotation.lo turbomachinery.lo vof.lo
$(top_srcdir)/src/base/vof.f90: cstphy.lo field.lo numvar.lo optcal.lo pointe.lo mesh.lo paramx.lo entsor.lo parall.lo rotation.lo turbomachinery.lo
$(top_srcdir)/src/base/cplsat.f90: paramx.lo
$(top_srcdir)/src/base/cstphy.f90: paramx.lo
$(top_srcdir)/src/base/cs_c_bindings.f90: paramx.lo mesh.lo period.lo field.lo field_operator.lo
$(top_srcdir)/src/base/cs_f_interfaces.f90: cs_c_bindings.lo paramx.lo dimens.lo mesh.lo optcal.lo field.lo lagran.lo
$(top_srcdir)/src/base/cs_nz_condensation.f90: pointe.lo
$(top_srcdir)/src/base/cs_nz_tagmr.f90: optcal.lo pointe.lo parall.lo cs_nz_condensation.lo cs_tagmr.lo
$(top_srcdir)/src/base/cs_tagmr.f90: pointe.lo
$(top_srcdir)/src/base/cs_tagms.f90: mesh.lo
$(top_srcdir)/src/base/entsor.f90: paramx.lo field.lo
$(top_srcdir)/src/base/field_operator.f90: paramx.lo dimens.lo mesh.lo optcal.lo field.lo
$(top_srcdir)/src/base/numvar.f90: paramx.lo
$(top_srcdir)/src/base/optcal.f90: paramx.lo field.lo numvar.lo cstphy.lo
$(top_srcdir)/src/base/pointe.f90: paramx.lo numvar.lo mesh.lo parall.lo period.lo optcal.lo entsor.lo lagran.lo radiat.lo ppincl.lo ihmpre.lo albase.lo cstnum.lo cs_c_bindings.lo
$(top_srcdir)/src/base/ptrglo.f90: paramx.lo dimens.lo mesh.lo numvar.lo optcal.lo parall.lo period.lo mesh.lo field.lo
$(top_srcdir)/src/base/rotation.f90:  paramx.lo
$(top_srcdir)/src/base/turbomachinery.f90:  mesh.lo cstphy.lo cplsat.lo
$(top_srcdir)/src/atmo/atincl.f90: ppppar.lo ppincl.lo mesh.lo atsoil.lo
$(top_srcdir)/src/atmo/atchem.f90: ppppar.lo mesh.lo entsor.lo numvar.lo
$(top_srcdir)/src/atmo/atimbr.f90: cstphy.lo entsor.lo atincl.lo
$(top_srcdir)/src/atmo/siream.f90: ppppar.lo
$(top_srcdir)/src/atmo/spefun.f90: entsor.lo
$(top_srcdir)/src/cogz/coincl.f90: paramx.lo ppppar.lo ppincl.lo
$(top_srcdir)/src/comb/cpincl.f90: ppppar.lo ppthch.lo
$(top_srcdir)/src/comb/cs_coal_incl.f90: ppppar.lo ppthch.lo
$(top_srcdir)/src/comb/cs_fuel_incl.f90: ppppar.lo ppthch.lo
$(top_srcdir)/src/ctwr/ctincl.f90: cstphy.lo
$(top_srcdir)/src/lagr/lagran.f90: cstnum.lo optcal.lo ppincl.lo radiat.lo cpincl.lo ppppar.lo cstphy.lo
$(top_srcdir)/src/pprt/ppcpfu.f90: ppthch.lo ppppar.lo
$(top_srcdir)/src/pprt/ppincl.f90: ppppar.lo ppthch.lo
$(top_srcdir)/src/pprt/ppthch.f90: cstphy.lo
$(top_srcdir)/src/rayt/radiat.f90: ppppar.lo ppincl.lo optcal.lo ppcpfu.lo numvar.lo
$(top_srcdir)/src/darc/darcy_module.f90: paramx.lo

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

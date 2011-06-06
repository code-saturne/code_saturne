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

SUBDIRS = . bft fvm mei gui base pprt atmo cfbl cogz cplv ctwr elec fuel lagr \
          rayt apps

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
$(top_srcdir)/src/base/alaste.f90 \
$(top_srcdir)/src/base/albase.f90 \
$(top_srcdir)/src/base/alstru.f90 \
$(top_srcdir)/src/base/cplsat.f90 \
$(top_srcdir)/src/base/cstphy.f90 \
$(top_srcdir)/src/base/cstnum.f90 \
$(top_srcdir)/src/base/dimens.f90 \
$(top_srcdir)/src/base/entsor.f90 \
$(top_srcdir)/src/base/ihmpre.f90 \
$(top_srcdir)/src/base/mltgrd.f90 \
$(top_srcdir)/src/base/numvar.f90 \
$(top_srcdir)/src/base/optcal.f90 \
$(top_srcdir)/src/base/parall.f90 \
$(top_srcdir)/src/base/period.f90 \
$(top_srcdir)/src/base/pointe.f90 \
$(top_srcdir)/src/base/vorinc.f90 \
$(top_srcdir)/src/atmo/atincl.f90 \
$(top_srcdir)/src/cfbl/cfpoin.f90 \
$(top_srcdir)/src/cogz/coincl.f90 \
$(top_srcdir)/src/cplv/cpincl.f90 \
$(top_srcdir)/src/ctwr/ctincl.f90 \
$(top_srcdir)/src/elec/elincl.f90 \
$(top_srcdir)/src/fuel/fuincl.f90 \
$(top_srcdir)/src/lagr/lagdim.f90 \
$(top_srcdir)/src/lagr/lagpar.f90 \
$(top_srcdir)/src/lagr/lagran.f90 \
$(top_srcdir)/src/pprt/ppcpfu.f90 \
$(top_srcdir)/src/pprt/ppincl.f90 \
$(top_srcdir)/src/pprt/ppppar.f90 \
$(top_srcdir)/src/pprt/ppthch.f90 \
$(top_srcdir)/src/rayt/radiat.f90

# Module dependencies (handled manually, as Automake does
# not yet handle Fortran module dependencies)

$(top_srcdir)/src/base/albase.f90: cplsat.lo
$(top_srcdir)/src/base/alstru.f90: paramx.lo
$(top_srcdir)/src/base/cplsat.f90: paramx.lo
$(top_srcdir)/src/base/cstphy.f90: paramx.lo
$(top_srcdir)/src/base/entsor.f90: paramx.lo
$(top_srcdir)/src/base/mltgrd.f90: paramx.lo
$(top_srcdir)/src/base/numvar.f90: paramx.lo
$(top_srcdir)/src/base/optcal.f90: paramx.lo
$(top_srcdir)/src/base/pointe.f90: paramx.lo parall.lo period.lo optcal.lo entsor.lo lagran.lo radiat.lo ppincl.lo ihmpre.lo albase.lo
$(top_srcdir)/src/atmo/atincl.f90: ppppar.lo
$(top_srcdir)/src/cogz/coincl.f90: paramx.lo ppppar.lo
$(top_srcdir)/src/cplv/cpincl.f90: ppppar.lo ppthch.lo
$(top_srcdir)/src/elec/elincl.f90: paramx.lo ppthch.lo
$(top_srcdir)/src/fuel/fuincl.f90: ppppar.lo ppthch.lo
$(top_srcdir)/src/lagr/lagran.f90: lagpar.lo
$(top_srcdir)/src/pprt/ppincl.f90: ppthch.lo
$(top_srcdir)/src/pprt/ppthch.f90: cstphy.lo
$(top_srcdir)/src/rayt/radiat.f90: ppppar.lo

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
	  echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(pkgincludedir)'"; \
	  $(INSTALL_HEADER) $$p "$(DESTDIR)$(pkgincludedir)" || exit $$?; \
	done

uninstall-local:
	@$(NORMAL_UNINSTALL)
	@list='$(libcsmodules_la_SOURCES)'; test -n "$(pkgincludedir)" || list=; \
	modlist=`for p in $$list; do echo $$p; done | $(SED) -e 's|.f90|.$(FCMODEXT)|g'`; \
	files=`for p in $$modlist; do echo $$p; done | $(SED) -e 's|^.*/||'`; \
	test -d "$(DESTDIR)$(pkgincludedir)" || exit 0; \
	test -n "$$files" || exit 0; \
	echo " ( cd '$(DESTDIR)$(pkgincludedir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(pkgincludedir)" && rm -f $$files
	rmdir $(pkgincludedir) || exit 0


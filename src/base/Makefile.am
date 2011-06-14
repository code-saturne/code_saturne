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

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I$(top_srcdir)/libple/src -I$(top_builddir)/libple
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/gui \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

AM_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCMODINCLUDE)$(top_builddir)/src/atmo \
$(FCMODINCLUDE)$(top_builddir)/src/cfbl \
$(FCMODINCLUDE)$(top_builddir)/src/cogz \
$(FCMODINCLUDE)$(top_builddir)/src/cplv \
$(FCMODINCLUDE)$(top_builddir)/src/ctwr \
$(FCMODINCLUDE)$(top_builddir)/src/elec \
$(FCMODINCLUDE)$(top_builddir)/src/fuel \
$(FCMODINCLUDE)$(top_builddir)/src/lagr \
$(FCMODINCLUDE)$(top_builddir)/src/pprt \
$(FCMODINCLUDE)$(top_builddir)/src/rayt \
$(FCFLAGS_DBG) $(FCFLAGS_OPT)

AM_LDFLAGS =

# Add user_modules file to distribution as it is not automatically accounted
# for by automake

EXTRA_DIST = $(top_srcdir)/users/base/user_modules.f90

# Public header files (to be installed)

pkginclude_HEADERS = \
cs_ale.h \
cs_ast_coupling.h \
cs_base.h \
cs_base_fortran.h \
cs_boundary_conditions.h \
cs_calcium.h \
cs_coupling.h \
cs_ext_neighborhood.h \
cs_halo.h \
cs_io.h \
cs_map.h \
cs_numbering.h \
cs_opts.h \
cs_parall.h \
cs_perio.h \
cs_post.h \
cs_preprocessor_data.h \
cs_prototypes.h \
cs_proxy_comm.h \
cs_renumber.h \
cs_resource.h \
cs_restart.h \
cs_sat_coupling.h \
cs_search.h \
cs_selector.h \
cs_sort.h \
cs_syr3_comm.h \
cs_syr3_coupling.h \
cs_syr3_messages.h \
cs_syr4_coupling.h \
cs_syr_coupling.h \
cs_time_plot.h \
cs_tpar1d.h \
cs_ventil.h

# Library source files

noinst_LTLIBRARIES = libcscore.la \
                     libcsbase.la

###

libcscore_la_SOURCES = \
cs_base.c \
cs_io.c \
cs_map.c
libcscore_la_LDFLAGS = -no-undefined

###

libcsbase_la_SOURCES = \
alelap.f90 \
alemaj.f90 \
altycl.f90 \
calhyd.f90 \
calmom.f90 \
caltri.f90 \
catsma.f90 \
causta.f90 \
clca66.f90 \
clpsca.f90 \
clptrg.f90 \
clptur.f90 \
clsyvt.f90 \
condli.f90 \
cou1di.f90 \
cou1do.f90 \
coupbi.f90 \
coupbo.f90 \
codits.f90 \
covofi.f90 \
cregeo.f90 \
cs_ale.c \
cs_ast_coupling.c \
cs_base_fortran.c \
cs_boundary_conditions.c \
cs_calcium.c \
csc2cl.f90 \
csc2ts.f90 \
csccel.f90 \
cscfbr.f90 \
cscini.f90 \
csclli.f90 \
cscloc.f90 \
cscpce.f90 \
cscpfb.f90 \
cs_coupling.c \
cs_ext_neighborhood.c \
csflsh.f90 \
cs_halo.c \
csinit.f90 \
cs_numbering.c \
csopli.f90 \
cs_opts.c \
cs_parall.c \
cs_perio.c \
cs_post.c \
cs_post_f2c.f90 \
cs_preprocessor_data.c \
cs_proxy_comm.c \
cs_renumber.c \
cs_resource.c \
cs_restart.c \
cs_sat_coupling.c \
cs_search.c \
cs_selector.c \
cs_selector_f2c.f90 \
cs_sort.c \
cs_syr3_comm.c \
cs_syr3_coupling.c \
cs_syr3_messages.c \
cs_syr4_coupling.c \
cs_syr_coupling.c \
cs_time_plot.c \
cs_tpar1d.c \
cs_ventil.c \
csprnt.f90 \
distpr.f90 \
distyp.f90 \
dttvar.f90 \
dvvpst.f90 \
ecrava.f90 \
ecrhis.f90 \
ecrlis.f90 \
findpt.f90 \
haltyp.f90 \
hturbp.f90 \
idrbla.f90 \
impini.f90 \
iniini.f90 \
initi1.f90 \
initi2.f90 \
iniusi.f90 \
iniva0.f90 \
inivar.f90 \
iprbla.f90 \
lecamo.f90 \
lecamp.f90 \
lecamx.f90 \
majgeo.f90 \
memfin.f90 \
memtri.f90 \
modini.f90 \
modpar.f90 \
mmtycl.f90 \
navsto.f90 \
newmrk.f90 \
numvec.f90 \
ordita.f90 \
pergra.f90 \
perinr.f90 \
perinu.f90 \
phyvar.f90 \
precli.f90 \
preduv.f90 \
prmoy0.f90 \
recvmc.f90 \
resolp.f90 \
scalai.f90 \
schtmp.f90 \
stdtcl.f90 \
strdep.f90 \
strhis.f90 \
strini.f90 \
strpre.f90 \
tdesi1.f90 \
testel.f90 \
tridim.f90 \
tsepdc.f90 \
tstvec.f90 \
typecl.f90 \
undscr.f90 \
varpos.f90 \
vericl.f90 \
verini.f90 \
verlon.f90 \
vert1d.f90 \
vissec.f90 \
zufall.f90


baseusers = \
$(top_srcdir)/users/base/usalcl.f90 \
$(top_srcdir)/users/base/usalin.f90 \
$(top_srcdir)/users/base/usaste.f90 \
$(top_srcdir)/users/base/usclim.f90 \
$(top_srcdir)/users/base/ushist.f90 \
$(top_srcdir)/users/base/usini1.f90 \
$(top_srcdir)/users/base/usiniv.f90 \
$(top_srcdir)/users/base/uskpdc.f90 \
$(top_srcdir)/users/base/usmpst.f90 \
$(top_srcdir)/users/base/usphyv.f90 \
$(top_srcdir)/users/base/usproj.f90 \
$(top_srcdir)/users/base/uspt1d.f90 \
$(top_srcdir)/users/base/ussmag.f90 \
$(top_srcdir)/users/base/usstru.f90 \
$(top_srcdir)/users/base/usthht.f90 \
$(top_srcdir)/users/base/ustmgr.f90 \
$(top_srcdir)/users/base/ustske.f90 \
$(top_srcdir)/users/base/ustskw.f90 \
$(top_srcdir)/users/base/ustsma.f90 \
$(top_srcdir)/users/base/ustsns.f90 \
$(top_srcdir)/users/base/ustsri.f90 \
$(top_srcdir)/users/base/ustssa.f90 \
$(top_srcdir)/users/base/ustssc.f90 \
$(top_srcdir)/users/base/ustsv2.f90 \
$(top_srcdir)/users/base/usvima.f90 \
$(top_srcdir)/users/base/usvist.f90 \
$(top_srcdir)/users/base/usvort.f90 \
$(top_srcdir)/users/base/usvpst.f90 \
$(top_srcdir)/users/base/cs_user_join.c \
$(top_srcdir)/users/base/cs_user_mesh_input.c \
$(top_srcdir)/users/base/cs_user_mesh_modify.c \
$(top_srcdir)/users/base/cs_user_periodicity.c \
$(top_srcdir)/users/base/cs_user_postprocess.c \
$(top_srcdir)/users/base/cs_user_saturne_coupling.c \
$(top_srcdir)/users/base/cs_user_solver.c \
$(top_srcdir)/users/base/cs_user_syrthes_coupling.c

libcsbase_la_SOURCES += $(baseusers)

basedir = $(pkgdatadir)/users/base
base_DATA = $(baseusers) \
$(top_srcdir)/users/base/user_modules.f90

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

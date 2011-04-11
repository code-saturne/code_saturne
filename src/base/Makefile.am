#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009-2010 EDF S.A., France
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
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/mei \
-I$(top_srcdir)/src/gui \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS) $(BLAS_CPPFLAGS)

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

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
cs_ale.h \
cs_ast_coupling.h \
cs_base.h \
cs_base_fortran.h \
cs_benchmark.h \
cs_blas.h \
cs_boundary_conditions.h \
cs_calcium.h \
cs_coupling.h \
cs_ext_neighborhood.h \
cs_gradient.h \
cs_grid.h \
cs_halo.h \
cs_io.h \
cs_join.h \
cs_join_intersect.h \
cs_join_merge.h \
cs_join_mesh.h \
cs_join_perio.h \
cs_join_post.h \
cs_join_set.h \
cs_join_split.h \
cs_join_update.h \
cs_join_util.h \
cs_map.h \
cs_matrix.h \
cs_mesh_coherency.h \
cs_mesh_connect.h \
cs_mesh.h \
cs_mesh_halo.h \
cs_mesh_quality.h \
cs_mesh_quantities.h \
cs_mesh_save.h \
cs_mesh_solcom.h \
cs_mesh_warping.h \
cs_multigrid.h \
cs_numbering.h \
cs_opts.h \
cs_parall.h \
cs_perio.h \
cs_post.h \
cs_preprocessor_data.h \
cs_prototypes.h \
cs_proxy_comm.h \
cs_renumber.h \
cs_restart.h \
cs_sat_coupling.h \
cs_search.h \
cs_selector.h \
cs_sles.h \
cs_sort.h \
cs_syr3_comm.h \
cs_syr3_coupling.h \
cs_syr3_messages.h \
cs_syr4_coupling.h \
cs_syr_coupling.h \
cs_time_plot.h \
cs_tpar1d.h \
cs_ventil.h \
tcpumx.h \
tremai.h

# Library source files

noinst_LTLIBRARIES = libcshot.la \
                     libcscore.la \
                     libcsbase.la

###

libcshot_la_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_HOT)
libcshot_la_FCFLAGS = \
$(FCMODINCLUDE)$(top_builddir)/src \
$(FCFLAGS_DBG) $(FCFLAGS_HOT)
libcshot_la_SOURCES = \
cs_benchmark.c \
cs_blas.c \
cs_matrix.c \
cs_sles.c \
gradmc.f90 \
gradrc.f90 \
prods2.f90 \
prods3.f90 \
prodsc.f90 \
promav.f90
libcshot_la_LDFLAGS = -no-undefined

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
armtps.f90 \
autmgr.f90 \
bilsc2.f90 \
calhyd.f90 \
calmom.f90 \
caltri.f90 \
catsma.f90 \
causta.f90 \
clca66.f90 \
clipke.f90 \
clprij.f90 \
clpsca.f90 \
clptrg.f90 \
clptur.f90 \
clpv2f.f90 \
clsyvt.f90 \
codits.f90 \
condli.f90 \
cou1di.f90 \
cou1do.f90 \
coupbi.f90 \
coupbo.f90 \
covofi.f90 \
cregeo.f90 \
crstgr.f90 \
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
cs_gradient.c \
cs_grid.c \
cs_halo.c \
csinit.f90 \
cs_join.c \
cs_join_intersect.c \
cs_join_merge.c \
cs_join_mesh.c \
cs_join_perio.c \
cs_join_post.c \
cs_join_set.c \
cs_join_split.c \
cs_join_update.c \
cs_join_util.c \
cs_mesh.c \
cs_mesh_coherency.c \
cs_mesh_connect.c \
cs_mesh_halo.c \
cs_mesh_quality.c \
cs_mesh_quantities.c \
cs_mesh_save.c \
cs_mesh_solcom.c \
cs_mesh_warping.c \
cs_multigrid.c \
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
divmas.f90 \
divrij.f90 \
dttvar.f90 \
dvvpst.f90 \
ecrava.f90 \
ecrhis.f90 \
ecrlis.f90 \
findpt.f90 \
grdcel.f90 \
haltyp.f90 \
hturbp.f90 \
idrbla.f90 \
impini.f90 \
iniini.f90 \
inimas.f90 \
inipst.f90 \
initi1.f90 \
initi2.f90 \
iniusi.f90 \
iniva0.f90 \
inivar.f90 \
invers.f90 \
iprbla.f90 \
itrgrp.f90 \
itrmas.f90 \
lecamo.f90 \
lecamp.f90 \
lecamx.f90 \
ledgeo.f90 \
letgeo.f90 \
majgeo.f90 \
matrdt.f90 \
matrix.f90 \
memale.f90 \
memcli.f90 \
memcs1.f90 \
memcs2.f90 \
memdis.f90 \
memdtv.f90 \
memdyp.f90 \
memfin.f90 \
memini.f90 \
memkep.f90 \
memkom.f90 \
memnav.f90 \
mempdc.f90 \
memphy.f90 \
memrij.f90 \
memsca.f90 \
memt1d.f90 \
memtri.f90 \
memtsm.f90 \
memv2f.f90 \
memvor.f90 \
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
projts.f90 \
proxav.f90 \
recvmc.f90 \
reseps.f90 \
resolp.f90 \
resrij.f90 \
resssg.f90 \
resv2f.f90 \
rijech.f90 \
rijthe.f90 \
scalai.f90 \
schtmp.f90 \
stdtcl.f90 \
strdep.f90 \
strhis.f90 \
strini.f90 \
strpre.f90 \
tcpumx.c \
tdesi1.f90 \
testel.f90 \
tremai.c \
tridim.f90 \
tsepdc.f90 \
tstvec.f90 \
turbke.f90 \
turbkw.f90 \
turent.f90 \
turrij.f90 \
typecl.f90 \
undscr.f90 \
vandri.f90 \
varpos.f90 \
vectds.f90 \
vericl.f90 \
verini.f90 \
verlon.f90 \
vert1d.f90 \
viscfa.f90 \
visdyn.f90 \
vislmg.f90 \
visort.f90 \
vissec.f90 \
vissma.f90 \
vissst.f90 \
viswal.f90 \
vor2cl.f90 \
vordep.f90 \
vorimp.f90 \
vorin0.f90 \
vorini.f90 \
vorlgv.f90 \
vorpre.f90 \
vortex.f90 \
vorver.f90 \
vorvit.f90 \
zufall.f90


baseusers = \
$(top_srcdir)/users/base/usalcl.f90 \
$(top_srcdir)/users/base/usalin.f90 \
$(top_srcdir)/users/base/usaste.f90 \
$(top_srcdir)/users/base/usclim.f90 \
$(top_srcdir)/users/base/usdpst.f90 \
$(top_srcdir)/users/base/ushist.f90 \
$(top_srcdir)/users/base/usini1.f90 \
$(top_srcdir)/users/base/usiniv.f90 \
$(top_srcdir)/users/base/uskpdc.f90 \
$(top_srcdir)/users/base/usmpst.f90 \
$(top_srcdir)/users/base/usnpst.f90 \
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
$(top_srcdir)/users/base/cs_user_saturne_coupling.c \
$(top_srcdir)/users/base/cs_user_solver.c \
$(top_srcdir)/users/base/cs_user_syrthes_coupling.c

libcsbase_la_SOURCES += $(baseusers)

basedir = $(pkgdatadir)/users/base
base_DATA = $(baseusers) \
$(top_srcdir)/users/base/user_modules.f90

clean-local:
	-rm -f *__genmod.f90 *__genmod.mod

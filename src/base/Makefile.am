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

# Main part

AM_CPPFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/ctwr \
-I$(top_srcdir)/include/lagr \
@BFT_CPPFLAGS@ @FVM_CPPFLAGS@ @MPI_CPPFLAGS@ @LIBXML2_CPPFLAGS@ @BLAS_CPPFLAGS@

AM_CFLAGS = @CFLAGS_DBG@ @CFLAGS_OPT@

AM_FCFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/pprt \
-I$(top_srcdir)/include/atmo \
-I$(top_srcdir)/include/cfbl \
-I$(top_srcdir)/include/cogz \
-I$(top_srcdir)/include/cplv \
-I$(top_srcdir)/include/ctwr \
-I$(top_srcdir)/include/elec \
-I$(top_srcdir)/include/fuel \
-I$(top_srcdir)/include/mati \
-I$(top_srcdir)/include/lagr \
-I$(top_srcdir)/include/rayt \
@FCFLAGS_DBG@ @FCFLAGS_OPT@

AM_LDFLAGS = @LDFLAGS_DBG@ @LDFLAGS_OPT@

# Public header files (to be installed)

saturneincludedir = $(includedir)
saturneinclude_HEADERS = \
$(top_srcdir)/include/base/alaste.h \
$(top_srcdir)/include/base/albase.h \
$(top_srcdir)/include/base/alstru.h \
$(top_srcdir)/include/base/cs_ale.h \
$(top_srcdir)/include/base/cs_ast_coupling.h \
$(top_srcdir)/include/base/cs_base.h \
$(top_srcdir)/include/base/cs_benchmark.h \
$(top_srcdir)/include/base/cs_blas.h \
$(top_srcdir)/include/base/cs_boundary_conditions.h \
$(top_srcdir)/include/base/cs_calcium.h \
$(top_srcdir)/include/base/cs_coupling.h \
$(top_srcdir)/include/base/cs_ext_neighborhood.h \
$(top_srcdir)/include/base/cs_gradient.h \
$(top_srcdir)/include/base/cs_grid.h \
$(top_srcdir)/include/base/cs_gui.h \
$(top_srcdir)/include/base/cs_gui_matisse.h \
$(top_srcdir)/include/base/cs_gui_radiative_transfer.h \
$(top_srcdir)/include/base/cs_gui_util.h \
$(top_srcdir)/include/base/cs_halo.h \
$(top_srcdir)/include/base/cs_io.h \
$(top_srcdir)/include/base/cs_matrix.h \
$(top_srcdir)/include/base/cs_mesh_coherency.h \
$(top_srcdir)/include/base/cs_mesh_connect.h \
$(top_srcdir)/include/base/cs_mesh.h \
$(top_srcdir)/include/base/cs_mesh_halo.h \
$(top_srcdir)/include/base/cs_mesh_quality.h \
$(top_srcdir)/include/base/cs_mesh_quantities.h \
$(top_srcdir)/include/base/cs_mesh_solcom.h \
$(top_srcdir)/include/base/cs_mesh_warping.h \
$(top_srcdir)/include/base/cs_multigrid.h \
$(top_srcdir)/include/base/cs_numbering.h \
$(top_srcdir)/include/base/cs_opts.h \
$(top_srcdir)/include/base/cs_parall.h \
$(top_srcdir)/include/base/cs_perio.h \
$(top_srcdir)/include/base/cs_post.h \
$(top_srcdir)/include/base/cs_preprocessor_data.h \
$(top_srcdir)/include/base/cs_prototypes.h \
$(top_srcdir)/include/base/cs_proxy_comm.h \
$(top_srcdir)/include/base/cs_renumber.h \
$(top_srcdir)/include/base/cs_restart.h \
$(top_srcdir)/include/base/cs_sat_coupling.h \
$(top_srcdir)/include/base/cs_selector.h \
$(top_srcdir)/include/base/cs_sles.h \
$(top_srcdir)/include/base/cs_syr3_comm.h \
$(top_srcdir)/include/base/cs_syr3_coupling.h \
$(top_srcdir)/include/base/cs_syr3_messages.h \
$(top_srcdir)/include/base/cs_syr4_coupling.h \
$(top_srcdir)/include/base/cs_syr_coupling.h \
$(top_srcdir)/include/base/cstnum.h \
$(top_srcdir)/include/base/cs_tpar1d.h \
$(top_srcdir)/include/base/cstphy.h \
$(top_srcdir)/include/base/cs_ventil.h \
$(top_srcdir)/include/base/dimens.h \
$(top_srcdir)/include/base/dimfbr.h \
$(top_srcdir)/include/base/entsor.h \
$(top_srcdir)/include/base/ihmpre.h \
$(top_srcdir)/include/base/mltgrd.h \
$(top_srcdir)/include/base/numvar.h \
$(top_srcdir)/include/base/optcal.h \
$(top_srcdir)/include/base/parall.h \
$(top_srcdir)/include/base/paramx.h \
$(top_srcdir)/include/base/period.h \
$(top_srcdir)/include/base/pointe.h \
$(top_srcdir)/include/base/tcpumx.h \
$(top_srcdir)/include/base/tremai.h \
$(top_srcdir)/include/base/vector.h \
$(top_srcdir)/include/base/vortex.h

# Library source files

lib_LTLIBRARIES = libsaturne.la
libsaturne_la_SOURCES =
libsaturne_la_LIBADD = \
$(top_builddir)/src/atmo/libcsatmo.la \
$(top_builddir)/src/cfbl/libcscfbl.la \
$(top_builddir)/src/cogz/libcscogz.la \
$(top_builddir)/src/cplv/libcscplv.la \
$(top_builddir)/src/ctwr/libcsctwr.la \
$(top_builddir)/src/elec/libcselec.la \
$(top_builddir)/src/fuel/libcsfuel.la \
$(top_builddir)/src/lagr/libcslagr.la \
$(top_builddir)/src/mati/libcsmati.la \
$(top_builddir)/src/pprt/libcspprt.la \
$(top_builddir)/src/rayt/libcsrayt.la
libsaturne_la_LDFLAGS = -no-undefined

####

noinst_LTLIBRARIES = libcshot.la \
                     libcsbase.la

libsaturne_la_LIBADD += libcshot.la
libcshot_la_CFLAGS = @CFLAGS_DBG@ @CFLAGS_HOT@
libcshot_la_FCFLAGS = \
-I$(top_srcdir)/include/base \
-I$(top_srcdir)/include/pprt \
-I$(top_srcdir)/include/atmo \
-I$(top_srcdir)/include/cfbl \
-I$(top_srcdir)/include/cogz \
-I$(top_srcdir)/include/cplv \
-I$(top_srcdir)/include/ctwr \
-I$(top_srcdir)/include/elec \
-I$(top_srcdir)/include/fuel \
-I$(top_srcdir)/include/mati \
-I$(top_srcdir)/include/lagr \
-I$(top_srcdir)/include/rayt \
@FCFLAGS_DBG@ @FCFLAGS_HOT@
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
libcshot_la_LDFLAGS = @LDFLAGS_DBG@ @LDFLAGS_HOT@ -no-undefined

###

#noinst_LTLIBRARIES = libcsbase.la
libsaturne_la_LIBADD += libcsbase.la
libcsbase_la_SOURCES = \
alelap.f90 \
alemaj.f90 \
altycl.f90 \
armtps.f90 \
autmgr.f90 \
bilsc2.f90 \
calgeo.f90 \
calhyd.f90 \
calmom.f90 \
caltri.f90 \
catsma.f90 \
causta.f90 \
clca66.f90 \
cldijp.f90 \
cldipo.f90 \
clipke.f90 \
clprij.f90 \
clpsca.f90 \
clptrg.f90 \
clptur.f90 \
clpv2f.f90 \
clsurn.f90 \
clsyvt.f90 \
clvolc.f90 \
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
cs_base.c \
cs_boundary_conditions.c \
cs_calcium.c \
csclli.f90 \
cs_coupling.c \
cs_ext_neighborhood.c \
csflsh.f90 \
cs_gradient.c \
cs_grid.c \
cs_gui.c \
cs_gui_matisse.c \
cs_gui_radiative_transfer.c \
cs_gui_util.c \
cs_halo.c \
csinit.f90 \
cs_interpr_util.f90 \
cs_io.c \
cs_main.c \
cs_mesh.c \
cs_mesh_coherency.c \
cs_mesh_connect.c \
cs_mesh_halo.c \
cs_mesh_quality.c \
cs_mesh_quantities.c \
cs_mesh_solcom.c \
cs_mesh_warping.c \
cs_multigrid.c \
cs_numbering.c \
cs_opts.c \
cs_parall.c \
cs_perio.c \
cs_post.c \
cs_post_util.f90 \
cs_preprocessor_data.c \
cs_proxy_comm.c \
cs_renumber.c \
cs_restart.c \
cs_sat_coupling.c \
cs_selector.c \
cs_syr3_comm.c \
cs_syr3_coupling.c \
cs_syr3_messages.c \
cs_syr4_coupling.c \
cs_syr_coupling.c \
cs_syr_util.f90 \
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
memcbr.f90 \
memclg.f90 \
memcli.f90 \
memdis.f90 \
memdtv.f90 \
memdyp.f90 \
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

libcsbase_la_SOURCES += \
$(top_srcdir)/users/base/usalcl.f90 \
$(top_srcdir)/users/base/usalin.f90 \
$(top_srcdir)/users/base/usaste.f90 \
$(top_srcdir)/users/base/usclim.f90 \
$(top_srcdir)/users/base/usdpst.f90 \
$(top_srcdir)/users/base/ushist.f90 \
$(top_srcdir)/users/base/usini1.f90 \
$(top_srcdir)/users/base/usiniv.f90 \
$(top_srcdir)/users/base/uskpdc.f90 \
$(top_srcdir)/users/base/usmodg.f90 \
$(top_srcdir)/users/base/usmpst.f90 \
$(top_srcdir)/users/base/usnpst.f90 \
$(top_srcdir)/users/base/usphyv.f90 \
$(top_srcdir)/users/base/usproj.f90 \
$(top_srcdir)/users/base/uspt1d.f90 \
$(top_srcdir)/users/base/usruet.f90 \
$(top_srcdir)/users/base/ussmag.f90 \
$(top_srcdir)/users/base/usstru.f90 \
$(top_srcdir)/users/base/ussyrc.f90 \
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
$(top_srcdir)/users/base/usvpst.f90

basedir = $(pkgdatadir)/users/base
base_DATA = \
$(top_srcdir)/users/base/usalcl.f90 \
$(top_srcdir)/users/base/usalin.f90 \
$(top_srcdir)/users/base/usaste.f90 \
$(top_srcdir)/users/base/usclim.f90 \
$(top_srcdir)/users/base/usdpst.f90 \
$(top_srcdir)/users/base/ushist.f90 \
$(top_srcdir)/users/base/usini1.f90 \
$(top_srcdir)/users/base/usiniv.f90 \
$(top_srcdir)/users/base/uskpdc.f90 \
$(top_srcdir)/users/base/usmodg.f90 \
$(top_srcdir)/users/base/usmpst.f90 \
$(top_srcdir)/users/base/usnpst.f90 \
$(top_srcdir)/users/base/usphyv.f90 \
$(top_srcdir)/users/base/usproj.f90 \
$(top_srcdir)/users/base/uspt1d.f90 \
$(top_srcdir)/users/base/usruet.f90 \
$(top_srcdir)/users/base/ussmag.f90 \
$(top_srcdir)/users/base/usstru.f90 \
$(top_srcdir)/users/base/ussyrc.f90 \
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
$(top_srcdir)/users/base/usvpst.f90

####

bin_PROGRAMS = cs_solver
cs_solver_SOURCES =
cs_solver_LDADD = \
$(top_builddir)/src/base/libsaturne.la \
@FVM_LDFLAGS@ @FVM_LIBS@ \
@BFT_LDFLAGS@ @BFT_LIBS@ \
@MPI_LDFLAGS@ @MPI_LIBS@ \
@LIBXML2_LDFLAGS@ @LIBXML2_LIBS@ \
@BLAS_LDFLAGS@ @BLAS_LIBS@

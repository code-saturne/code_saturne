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

# Xfig files to process

FIGFILES = \
bordortho.fig \
cell_cog.fig \
cell_cog_nonorth.fig \
condlimite.fig \
domaine1d.fig \
entree.fig \
facebord.fig \
face_quant.fig \
face_split_1.fig \
face_split_2.fig \
face_split_delaunay_1.fig \
face_split_delaunay_2.fig \
face_split_delaunay_crit.fig \
face_split_main.fig \
face_surf.fig \
facesym.fig \
facette.fig \
fluxbord.fig \
hf.fig \
join_edge_inter_3d_eps.fig \
join_edge_inter_3d.fig \
join_edge_inter_3d_type.fig \
join_merge_1.fig \
join_merge_2.fig \
join_merge_3.fig \
join_merge_4.fig \
join_merge_5.fig \
join_non_planar.fig \
join_overlap_3d_1.fig \
join_overlap_3d_2.fig \
join_overlap_3d_3.fig \
join_perio.fig \
join_possible.fig \
join_principle.fig \
join_simplify_1.fig \
join_simplify_2.fig \
join_tolerance.fig \
LW2.fig \
LW4.fig \
Pf2.fig \
Pf.fig \
support.fig \
TONOx.fig \
Yf0.fig \
Yf1.fig \
Yf2.fig \
Yf3.fig \
Yf4.fig \
Yf.fig

# Other files

OTHERFILES = \
clthermique.dat \
clthermique.pdf \
clthermique.py \
clthermique.xmgr \
EK_model.png \
loilog.pdf \
loilog.py \
loilog.xmgr

# Distributed files

EXTRA_DIST = $(FIGFILES) $(OTHERFILES)

# Generated documents (pdf files)

PDF_FIGFILES = $(FIGFILES:.fig=.pdf)

# New suffixes and targets

SUFFIXES = .fig .pdf

# On Windows, fig2dev may be broken and generate a file named "-" instead
# of the target. Thus, we try to renamed such files if they exist.

.fig.pdf:
	"$(FIG2DEV)" -L pdf $< > $@
	@if test -f - ; then \mv -f - $@ ; fi

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local: $(PDF_FIGFILES)

clean-local:
	-rm -f $(PDF_FIGFILES)

distclean-local: clean-local

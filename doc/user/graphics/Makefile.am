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

# Xfig files to process

FIGFILES = \
cs_components.fig \
join_periodic.fig \
join_tolerance.fig \
vortex.fig

# Other files

OTHERFILES = \
halo.pdf \
rota_perio_parall.jpg \
ale1.png \
ale2.png \
ale3.png \
ale4.png \
ale5.png \
atmo1.png \
cl-coal1.png \
density_law.png \
global-Lag3.png \
GUI_10.png \
GUI_11.png \
GUI_12.png \
GUI_13.png \
GUI_14.png \
GUI_15.png \
GUI_16.png \
GUI_17.png \
GUI_19.png \
GUI_20.png \
GUI_21.png \
GUI_22.png \
GUI_23.png \
GUI_24.png \
GUI_25.png \
GUI_26.png \
GUI_27.png \
GUI_28.png \
GUI_29.png \
GUI_2.png \
GUI_3.png \
GUI_4.png \
GUI_5.png \
GUI_6.png \
GUI_7.png \
GUI_8.png \
GUI_9.png \
head_loss0.png \
head_loss3.png \
head_loss5.png \
homogeneous-coal1.png \
homogeneous-coal2.png \
homogeneous-coal3.png \
homogeneous-coal4.png \
homogeneous-coal5.png \
join.png \
Lag-cl.png \
Lag-global_settings2.png \
Lag-global_settings3.png \
Lag-output.png \
Lag-Statistics.png \
new_1.png \
periodicity.png \
rad_transf0.png \
rad_transf1.png \
rad_transf2.png \
rad_transf3.png \
rad_transf4.png \
rad_transf5.png \
rad_transf6.png \
rad_transf7.png \
solution_control-coal.png \
syrthes_coupling1.png \
time_average-coal.png \
V-2.png \
V-39.png \
V-3.png \
V-60.png

# Distributed files

EXTRA_DIST = $(FIGFILES) $(OTHERFILES)

# Generated documents (pdf files)

PDF_FIGFILES = $(FIGFILES:.fig=.pdf)

# New suffixes and targets

SUFFIXES = .fig .pdf

.fig.pdf:
	fig2dev -L pdf $< > $@

# One redefines here the standard pdf targets of automake
# so that they handle generation from LaTeX documents.

pdf-local: $(PDF_FIGFILES)

clean-local:
	-rm -f $(PDF_FIGFILES)

distclean-local: clean-local

## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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

FIGFILES =

# Other files

OTHERFILES = \
c2_p1.jpg c2_p2.jpg c2_p3.jpg c2_p4.jpg c2_p5.jpg c2_p6.jpg c2_p7.jpg \
c5_config.jpg cas1_t_10.jpg cas1_t_1.jpg cas1_t_20.jpg cas1_t_29.jpg \
case3_p10.jpg case3_p1.jpg case3_p2.jpg case3_p3.jpg case3_p4.jpg \
case3_p5.jpg case3_p6.jpg case3_p7.jpg case3_p8.jpg case3_p9.jpg \
case4_p1.jpg case4_p2.jpg case4_p3.jpg case4_p4.jpg case4_p5.jpg \
case4_p6.jpg case5_01.jpg case5_02.jpg case5_03.jpg case5_04.jpg \
case5_05.jpg case5_06.jpg case5_07.jpg case5_08.jpg case5_09.jpg \
color_Snimbus.jpg density_law.png evotemp_case3.png \
fig01.jpg fig02.jpg fig03.jpg fig04.jpg fig05.jpg fig06.jpg fig07.jpg \
fig08.jpg fig09.jpg \
head_loss0.png head_loss1.png head_loss2.png head_loss3.png head_loss4.png \
head_loss5.png \
moytemp.png probe9.png \
V-10.png V-11.png V-12.png V-13.png V-14.png V-15.png V-16.png V-17.png \
V-18.png V-1.png V-20.png V-21.png V-22.png V-23.png V-24.png V-25.png \
V-27.png V-28.png V-29.png V-2.png V-30.png V-31.png V-32.png V-33.png \
V-34bis.png V-34.png V-35.png V-36.png V-37.png V-38.png V-39.png V-3.png \
V-40.png V-41.png V-43.png V-44.png V-45.png V-46.png V-47.png V-48.png \
V-49.png V-4.png V-50.png V-51.png V-52.png V-53.png V-54.png V-56.png \
V-57.png V-58.png V-59.png V-5.png V-60bis.png V-60.png V-61.png V-62.png \
V-63.png V-64.png V-65.png V-66.png V-67.png V-68.png V-69bis.png V-69.png \
V-6.png V-70.png V-71.png V-72.png V-73bis.png V-73.png V-74bis.png V-74.png \
V-75.png V-76bis.png V-76.png V-7.png V-8.png V-9.png \
viscosity_law.png \
case3.cmd evotemp_case3.agr probe9.agr

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

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

FIGFILES =

# Other files

OTHERFILES = \
c1_capture01.jpg c1_capture02.jpg c1_capture03.jpg c1_capture04.jpg \
c1_capture05.jpg c1_capture06.jpg c1_capture07.jpg c1_capture08.jpg \
c1_capture09.jpg c1_capture10.jpg c1_capture11.jpg c1_capture12.jpg \
c1_capture13.jpg c1_capture14.jpg c1_capture15.jpg c1_capture16.jpg \
c1_capture17.jpg c1_capture18.jpg c1_capture19.jpg c1_capture20.jpg \
c1_capture21.jpg c1_capture22.jpg c1_capture23.jpg c1_capture24.jpg \
c1_capture25.jpg c1_capture26.jpg c1_capture27.jpg c1_capture28.jpg \
c1_capture29.jpg c1_capture30.jpg c1_capture31.jpg c1_capture32.jpg \
c1_capture33.jpg c1_capture34.jpg c1_capture35.jpg c1_capture36.jpg \
c1_capture37.jpg c1_capture38.jpg c1_capture39.jpg c1_capture40.jpg \
c1_capture41.jpg c1_capture42.jpg c1_capture43.jpg \
c2_capture01.jpg c2_capture02.jpg c2_capture03.jpg c2_capture04.jpg \
c2_capture05.jpg c2_capture06.jpg c2_capture07.jpg c2_capture08.jpg \
c2_capture09.jpg c2_capture10.jpg c2_capture11.jpg c2_capture12.jpg \
c2_capture13.jpg c2_capture14.jpg c2_capture152.jpg c2_capture15.jpg \
c2_capture16.jpg c2_capture17.jpg c2_capture18.jpg c2_capture19.jpg \
c2_capture20.jpg c2_capture212.jpg c2_capture21.jpg c2_capture22.jpg \
c2_capture23.jpg c2_capture24.jpg c2_capture25.jpg c2_capture26.jpg \
c2_capture27.jpg c2_capture28.jpg c2_capture29.jpg c2_capture30.jpg \
c2_p1.jpg c2_p2.jpg c2_p3.jpg c2_p4.jpg c2_p5.jpg c2_p6.jpg c2_p7.jpg \
c3_capture01.jpg c3_capture02.jpg c3_capture03.jpg c3_capture04.jpg \
c3_capture05.jpg c3_capture06.jpg c4_capture01.jpg c4_capture02.jpg \
c5_capture01.jpg c5_capture02.jpg c5_capture03.jpg c5_capture04.jpg \
c5_capture05.jpg c5_capture06.jpg c5_capture07.jpg c5_capture08.jpg \
c5_capture09.jpg c5_capture10.jpg c5_capture11.jpg c5_capture12.jpg \
c5_capture13.jpg c5_config.jpg \
capture10.jpg capture11.jpg capture12.jpg capture13.jpg \
cas1_t_10.jpg cas1_t_1.jpg cas1_t_20.jpg cas1_t_29.jpg \
case3_p10.jpg case3_p1.jpg case3_p2.jpg case3_p3.jpg case3_p4.jpg \
case3_p5.jpg case3_p6.jpg case3_p7.jpg case3_p8.jpg case3_p9.jpg \
case4_p1.jpg case4_p2.jpg case4_p3.jpg case4_p4.jpg case4_p5.jpg \
case4_p6.jpg case5_01.jpg case5_02.jpg case5_03.jpg case5_04.jpg \
case5_05.jpg case5_06.jpg case5_07.jpg case5_08.jpg case5_09.jpg \
color_Snimbus.jpg evotemp_case3.png \
fig01.jpg fig02.jpg fig03.jpg fig04.jpg fig05.jpg fig06.jpg \
fig07.jpg fig08.jpg fig09.jpg moytemp.png probe9.png

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

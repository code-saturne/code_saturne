/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*-----------------------------------------------------------------------------*/

/*!
  \page user_coupling Code coupling

  \section cs_user_coupling_h_intro Introduction

  C user functions for code coupling.

  Several functions are present in the file, each specific to different
  kind of couplings.

  \section cs_user_coupling_h_cs_user_coupling  Global options for coupling

  One can define global options for coupling with the user function
  \ref cs_user_coupling . These options allow defining the time step synchronization policy,
  as well as a time step multiplier.

  \snippet cs_user_coupling.c coupling_1

  \section cs_user_coupling_h_cs_user_syrthes_coupling Code coupling with SYRTHES

  The \ref cs_user_syrthes_coupling subroutine defines a or multiple couplings
  between Code_Saturne and SYRTHES by calling the \ref cs_syr_coupling_define
  function for each coupling to add.

  The following lines of code show different examples of coupling with SYRTHES.

  \subsection cs_user_coupling_h_cs_user_syrthes_coupling_example_1 Example 1

  \snippet cs_user_coupling-syrthes.c coupling_syrthes_1

  \subsection cs_user_coupling_h_cs_user_syrthes_coupling_example_2 Example 2

  \snippet cs_user_coupling-syrthes.c coupling_syrthes_2

  \subsection cs_user_coupling_h_cs_user_syrthes_coupling_example_3 Example 3

  \snippet cs_user_coupling-syrthes.c coupling_syrthes_3

  \section cs_user_coupling_h_cs_user_saturne_coupling Code coupling with other instances of Code_Saturne

  The \ref cs_user_saturne_coupling allows one to couple different instances of
  Code_Saturne by calling the \ref cs_sat_coupling_define function for each
  coupling to add.

  Two examples are provided hereafter.

  \subsection cs_user_coupling_h_cs_user_saturne_coupling_example_1 Example 1

  \snippet cs_user_coupling-saturne.c coupling_saturne_1

  \subsection cs_user_coupling_h_cs_user_saturne_coupling_example_2 Example 2

  \snippet cs_user_coupling-saturne.c coupling_saturne_2

*/

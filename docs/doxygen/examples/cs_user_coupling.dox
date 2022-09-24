/*============================================================================
 * code_saturne documentation page
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
  \ref cs_user_parameters. using \ref cs_coupling_set_sync_flag
  allows defining the time step
  synchronization policy:

  \snippet cs_user_parameters-coupling.c coupling_ts

  A time step multiplier between coupled tools may also be defined.
  The apparent time step for the current instance times (as viewed by
  coupled codes) is equal to the true time step times this multiplier.

  When coupling with SYRTHES, it is recommended to use the same multiplier
  here as for the thermal variable time step (this is not automated,
  so as to allow for more advanced combinations if necessary, so the user
  should ensure this when using a time step multiplier). For example:

  \snippet cs_user_parameters-coupling.c coupling_1

  \section cs_user_coupling_h_cs_user_syrthes_coupling Code coupling with SYRTHES

  The \ref cs_user_syrthes_coupling subroutine defines a or multiple couplings
  between code_saturne and SYRTHES by calling the \ref cs_syr_coupling_define
  function for each coupling to add.

  The following lines of code show different examples of coupling with SYRTHES.

  \subsection cs_user_coupling_h_cs_user_syrthes_coupling_example_1 Example 1

  \snippet cs_user_coupling-syrthes.c coupling_syrthes_1

  \subsection cs_user_coupling_h_cs_user_syrthes_coupling_example_2 Example 2

  \snippet cs_user_coupling-syrthes.c coupling_syrthes_2

  \subsection cs_user_coupling_h_cs_user_syrthes_coupling_example_3 Example 3

  \snippet cs_user_coupling-syrthes.c coupling_syrthes_3

  \section cs_user_coupling_h_cs_user_saturne_coupling Code coupling with other instances of code_saturne

  The \ref cs_user_saturne_coupling allows one to couple different instances of
  code_saturne by calling the \ref cs_sat_coupling_define function for each
  coupling to add.

  Two examples are provided hereafter.

  \subsection cs_user_coupling_h_cs_user_saturne_coupling_example_1 Example 1

  \snippet cs_user_coupling-saturne.c coupling_saturne_1

  \subsection cs_user_coupling_h_cs_user_saturne_coupling_example_2 Example 2

  \snippet cs_user_coupling-saturne.c coupling_saturne_2

*/

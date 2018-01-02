/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
  \page les_inflow Generation of synthetic turbulence at LES inlets
 
  \section cs_user_les_infow_intro Introduction

  This example contains 3 subroutines : 
  - cs_user_les_inflow_init: definition of global caracteristics of synthetic turbulence inlets
  - cs_user_les_inflow_define: definition of the caracteristics of the synthetic turbulence 
  - cs_user_les_inflow_advanced: accurate definition of target statistics at inlet


  \section cs_user_les_inflow_init Global caracteristics of synthetic turbulence inlets

  \subsection cs_user_les_inflow_purpose Purpose

  Generation of synthetic turbulence at LES inlets.
  Definition of global caracteristics of synthetic turbulence inlets:
  nent and isuisy might be defined.

  - nent = Number of inlets
  - isuisy = 1: Reading of the LES inflow module restart file
  - isuisy = 0: not activated (synthetic turbulence reinitialized)


  \subsection cs_user_les_inflow_loc_var_dec1 Local variables declaration

  No local variable.


  \subsection cs_user_les_inflow_init1 Initializations
  \snippet cs_user_les_inflow-base.f90 init_1


  \section cs_user_les_inflow_define Caracteristics of one specific inlet

  \subsection cs_user_les_inflow_purpose_base Purpose

  Generation of synthetic turbulence at LES inlets
  Definition of the caracteristics of the synthetic turbulence inlet 'nument'
  For each LES inlet, the following parameters might be defined:
   -# Data relatve to the method employed
      - typent indicates the synthetic turbulence method:
          - 0 : laminar, no turbulent fluctations
          - 1 : random gaussian noise
          - 2 : Batten method, based on Fourier mode decomposition
          - 3 : Synthetic Eddy Method (SEM)
      - nelent indicates the number of "entities" relative to the method
        (useful only for the Batten method and the SEM):
          - for Batten : number of Fourier modes of the turbulent fluctuations
          - for SEM    : number of synthetic eddies building the fluctuations
      - iverbo indicates the verbosity level (listing)
          - = 0  no specific output
          - > 0 additionnal output (only for SEM)
   -# Data relative to the LES inflow boundary faces
       - nfbent: number of boundary faces of the LES inflow
       - lfbent: list of boundary faces of the LES inflow
   -# Data relative to the flow
       - vitent(3): reference mean velocity vector
       - enrent: reference turbulent kinetic energy
       - dspent: reference dissipation rate
       - Note:
         - dspent useful only for typent = 2 (Batten) or typent = 3 (SEM).
         - Strictly positive values are required for enrent and dspent.
         - Accurate specification of the statistics of the flow at LES inlet
           can be made via the user subroutine cs_user_les_inflow_advanced.

  \subsection cs_user_les_inflow_loc_var_dec2 Local variables declaration

  No local variable.

  \subsection cs_user_les_inflow_init2 Initializations

  First synthetic turbulence inlet: the Batten Method is used
  for boundary faces of color '1'.

  \snippet cs_user_les_inflow-base.f90 init_21

  Second synthetic turbulence inlet: the Synthetic Eddy Method is used
  for the boundary faces verifying a geometric criterion.

  \snippet cs_user_les_inflow-base.f90 init_22


  \section cs_user_les_inflow_advanced Accurate specification of target statistics at inlet

  \subsection cs_user_les_inflow_purpose_advanced Purpose

  Generation of synthetic turbulence at LES inlets.
  Accurate definition of mean velocity, Reynolds stresses and dissipation
  rate for each boundary faces of the synthetic turbulence inlet 'nument'. 

  Usage:
   - uvwent(ndim,nfbent) : mean velocity vector
   - rijent(   6,nfbent) : Reynolds stresses!
   - epsent(     nfbent) : dissipation rate


  \subsection cs_user_les_inflow_loc_var_dec3 Local variables declaration

  \snippet cs_user_les_inflow-base.f90 loc_var_dec3

  \subsection cs_user_les_inflow_example_1 Example 1

  Mean velocity, Reynolds stresses an dissipation are deduced from a wall law for
  the first synthetic turbulence inlet, 
  - no refining of the statistics of the flow 
  is provided for the second synthetic turbulence inlet.
  
  \snippet cs_user_les_inflow-base.f90 example_1
  
  \subsection cs_user_les_inflow_example_2 Example 2

   Reynolds stresses and dissipation at the inlet are computed
   using the turbulence intensity and standard laws for
   a circular pipe for the first synthetic turbulence inlet,
   - no refining of the statistics of the flow is provided for the
   other synthetic turbulence inlet.
  
  \snippet cs_user_les_inflow-base.f90 example_2
  
*/

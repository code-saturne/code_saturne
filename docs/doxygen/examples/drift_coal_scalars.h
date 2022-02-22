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
  \page drift_coal_scalars Data setting for drift scalars for coal combustion
 

  \section drift_coal_scalars_intro Introduction

  This page provides an example of code blocks that may be used
  to perform a calculation with drift scalars for coal combustion.

  \section cs_user_physical_properties-coal_drift Physical properties

  \subsection drift_coal_scalars_loc_var Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_physical_properties-coal_drift.f90 loc_var_dec

  \subsection drift_coal_scalars_init Initialization and finalization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_physical_properties-coal_drift.f90 init

  In theory Fortran 95 deallocates locally-allocated arrays automatically,
  but deallocating arrays in a symmetric manner to their allocation is good
  practice, and it avoids using a different logic for C and Fortran.

  \subsection drift_coal_scalars_body Body

  Here is the corresponding code:
  
  \snippet cs_user_physical_properties-coal_drift.f90 example_1
  
*/

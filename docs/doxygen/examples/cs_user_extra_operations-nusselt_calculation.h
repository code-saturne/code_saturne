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

  \page cs_user_extra_operations-nusselt_calculation Calculation of a local Nusselt number

  \brief This function is called at the end of each time step, and has a very general purpose (\c i.e. anything that does not have another dedicated user subroutine)


  \section cs_f_user_extra_operations cs_f_user_extra_operations subroutine

  \subsection loc_var_f_user Local variables

  \snippet cs_user_extra_operations-nusselt_calculation.f90 loc_var_f_user

  \subsection nusselt_number Calculation of the Nusselt number

  \snippet cs_user_extra_operations-nusselt_calculation.f90  nusselt_number

  \subsubsection compute_nusselt Compute value reconstructed at I' for boundary faces
  
  \snippet cs_user_extra_operations-nusselt_calculation.f90 compute_nusselt

   \subsubsection general_nusselt General cas (for non-orthogonal meshes)
  \snippet cs_user_extra_operations-nusselt_calculation.f90 gen_nusselt

   Allocate a work array for the gradient calculation

  \snippet  cs_user_extra_operations-nusselt_calculation.f90  allocate_nusselt

  Compute gradient

  \snippet  cs_user_extra_operations-nusselt_calculation.f90 gradient_nusselt

  Compute reconstructed value in boundary cells

 \snippet  cs_user_extra_operations-nusselt_calculation.f90  value_nusselt

  Free memory

\snippet  cs_user_extra_operations-nusselt_calculation.f90 free_nusselt

\subsubsection orthogonal_nusselt Case of orthogonal meshes 

\snippet  cs_user_extra_operations-nusselt_calculation.f90 else_nusselt 

Compute reconstructed value 

\note Here, we assign the non-reconstructed value.

\snippet  cs_user_extra_operations-nusselt_calculation.f90 value_ortho_nusselt

Calculation of the bulk temperature

\snippet  cs_user_extra_operations-nusselt_calculation.f90  bulk_nusselt

\section sortc2  sortc2 subroutine 
\subsection loc_var_sortc2 Local variables
\snippet  cs_user_extra_operations-nusselt_calculation.f90 loc_var_sortc2
 
\subsection body_sortc2 Body
\snippet cs_user_extra_operations-nusselt_calculation.f90 body_sortc2 

*/  

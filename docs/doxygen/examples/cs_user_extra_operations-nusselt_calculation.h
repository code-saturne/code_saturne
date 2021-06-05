/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

  \page cs_user_extra_operations_examples_nusselt_calculation_p Calculation of a local Nusselt number

  \brief This function is called at the end of each time step, and has a very general purpose (\c i.e. anything that does not have another dedicated user subroutine)


  \section cs_f_user_extra_operations cs_f_user_extra_operations subroutine

  \subsection loc_var_f_user Local variables

  \snippet cs_user_extra_operations-nusselt_calculation.c loc_var_f_user

  \subsection nusselt_number Calculation of the Nusselt number

  \snippet cs_user_extra_operations-nusselt_calculation.c  nusselt_number

  \subsubsection compute_nusselt Compute value reconstructed at I' for boundary faces

  \snippet cs_user_extra_operations-nusselt_calculation.c compute_nusselt

  \subsubsection general_nusselt General case (for non-orthogonal meshes)
  Allocate a work array for the gradient calculation, then
  compute gradient, then compute reconstructed value in boundary cells,
  and then free memory
  \snippet cs_user_extra_operations-nusselt_calculation.c gen_nusselt

  \subsubsection orthogonal_nusselt Case of orthogonal meshes

  Compute boundary value without reconstruction
  \snippet cs_user_extra_operations-nusselt_calculation.c else_nusselt

  \note Here, we assign the non-reconstructed value.

  Open file to print values and perform parallel operations to
  broadcast values to all ranks.
  \snippet  cs_user_extra_operations-nusselt_calculation.c value_ortho_nusselt

  Calculation of the bulk temperature, finalize the Nusselt number, print it and
  free memory

  \snippet cs_user_extra_operations-nusselt_calculation.c  bulk_nusselt

\section sortc2 Utility function to sort global data

\subsection body_sortc2 Body
\snippet cs_user_extra_operations-nusselt_calculation.c body_sortc2

*/

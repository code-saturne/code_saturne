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

  \page cs_user_extra_operations_examples_nusselt_calculation_p Calculation of a local Nusselt number

  This is an example of \ref cs_user_extra_operations function which computes a Nusselt number.
  balances on specified zones.

  \subsection loc_var_f_user Local variables

  Faces can be selected either using a selector (as is done here) or directly
  using boundary zone element ids if an appropriate zone has been defined.

  \snippet cs_user_extra_operations-nusselt_calculation.c loc_var_f_user

  \subsection nusselt_number Calculation of the Nusselt number

  \snippet cs_user_extra_operations-nusselt_calculation.c  nusselt_number

  \subsubsection compute_nusselt Reconstruct value at selected boundary faces

  Allocate a local array for the selected boundary faces.
  \snippet cs_user_extra_operations-nusselt_calculation.c compute_nusselt

  \subsubsection general_nusselt General case (for non-orthogonal meshes)
  \snippet cs_user_extra_operations-nusselt_calculation.c gen_nusselt

  \subsubsection orthogonal_nusselt Case of orthogonal meshes

  Compute boundary value without reconstruction
  \snippet cs_user_extra_operations-nusselt_calculation.c else_nusselt

  \note Here, we assign the non-reconstructed value.

  Open file to print values and broadcast values to all parallel ranks.
  Values are ordered by the \c xabs array
  values provided to the \ref cs_parall_allgather_ordered_r, so this function
  is needed even when not running in parallel.

  \snippet  cs_user_extra_operations-nusselt_calculation.c value_ortho_nusselt

  Compute the bulk temperature, finalize the Nusselt number, print it and
  free memory not already freed before:

  \snippet cs_user_extra_operations-nusselt_calculation.c  bulk_nusselt

*/

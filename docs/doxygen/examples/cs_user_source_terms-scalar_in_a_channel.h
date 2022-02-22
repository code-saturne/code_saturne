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

  \page cs_user_source_terms-scalar_in_a_channel Examples of data settings for source terms with scalar in a channel

\brief Additional right-hand side source terms for scalar equations (user
 scalars and specific physics scalars) with  the \ref cs_user_source_terms
user-defined function.

\section loc_var_scal Local variables and initialization

\snippet cs_user_source_terms-scalar_in_a_channel.c st_meta

\section thermal_scalar_only Only apply to thermal scalar

\snippet cs_user_source_terms-scalar_in_a_channel.c thermal_scalar_only

\section func_body Function body

 Map required fields
\snippet cs_user_source_terms-scalar_in_a_channel.c map_fields

Compute bulk mean velocity
\snippet cs_user_source_terms-scalar_in_a_channel.c bulk_mean_velocity

Compute source terms; we want to impose a total flux of 1 Watt.
 \snippet cs_user_source_terms-scalar_in_a_channel.c scalar_st

*/

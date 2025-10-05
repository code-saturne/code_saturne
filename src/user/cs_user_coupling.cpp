/*============================================================================
 * Code couplings definition with SYRTHES, code_saturne., and CATHARE.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Define couplings with other instances of code_saturne.
 *
 * This is done by calling the \ref cs_sat_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_saturne_coupling
void
cs_user_saturne_coupling(void)
{
}

/*----------------------------------------------------------------------------*/
/*
 * Define couplings with SYRTHES code.
 *
 * This is done by calling the \ref cs_syr_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_syrthes_coupling
void
cs_user_syrthes_coupling(void)
{
}

/*----------------------------------------------------------------------------*/
/*
 * Compute a volume exchange coefficient for SYRTHES couplings.
 *
 * \param[in]   coupling_id   Syrthes coupling id
 * \param[in]   syrthes_name  name of associated Syrthes instance
 * \param[in]   n_elts        number of associated cells
 * \param[in]   elt_ids       associated cell ids
 * \param[out]  h_vol         associated exchange coefficient (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_syrthes_coupling_volume_h
void
cs_user_syrthes_coupling_volume_h
(
  [[maybe_unused]] int               coupling_id,
  [[maybe_unused]] const char       *syrthes_name,
  [[maybe_unused]] cs_lnum_t         n_elts,
  [[maybe_unused]] const cs_lnum_t   elt_ids[],
  [[maybe_unused]] cs_real_t         h_vol[]
)
{
}

/*----------------------------------------------------------------------------*/
/*
 * Define couplings with CATHARE code.
 *
 * This is done by calling the \ref cs_sys_coupling_add function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_cathare_coupling
void
cs_user_cathare_coupling(void)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

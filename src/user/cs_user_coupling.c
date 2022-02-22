/*============================================================================
 * Code couplings definition with SYRTHES and code_saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of code_saturne
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_coupling.c
 *
 * \brief Code couplings definition with SYRTHES and code_saturne.
 *
 * See \ref user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with other instances of code_saturne.
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
/*!
 * \brief Define couplings with SYRTHES code.
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
/*!
 * \brief Compute a volume exchange coefficient for SYRTHES couplings.
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
cs_user_syrthes_coupling_volume_h(int               coupling_id,
                                  const char       *syrthes_name,
                                  cs_lnum_t         n_elts,
                                  const cs_lnum_t   elt_ids[],
                                  cs_real_t         h_vol[])
{
  CS_UNUSED(coupling_id);
  CS_UNUSED(syrthes_name);
  CS_UNUSED(n_elts);
  CS_UNUSED(elt_ids);
  CS_UNUSED(h_vol);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

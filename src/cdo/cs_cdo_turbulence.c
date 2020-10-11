/*============================================================================
 * Routines to handle the resolution of the turbulence modelling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_turbulence.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_cdo_turbulence.c
 *
 *  \brief  Routines to handle the resoltion of the turbulence modelling
 *          within the CDO framework
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_CDO_TURBULENCE_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the structure managing the turbulence modelling
 *
 * \return a pointer to a new allocated cs_cdo_turbulence_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_turbulence_t *
cs_cdo_turbulence_create(void)
{
  cs_cdo_turbulence_t  *turb = NULL;

  BFT_MALLOC(turb, 1, cs_cdo_turbulence_t);

  /* The following structures are shared with the Legacy part */
  turb->model_param = NULL;
  turb->rans_modelling = NULL;
  turb->les_modelling = NULL;
  turb->reference_values = NULL;

  return turb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the structure managing the turbulence modelling
 *
 * \param[in, out]  p_turb_struct   pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_turbulence_free(cs_cdo_turbulence_t   **p_turb_struct)
{
  cs_cdo_turbulence_t  *turb = *p_turb_struct;

  /* The following structures are shared with the Legacy part. So, there is
   * no need to free the structure here */

  turb->model_param = NULL;
  turb->rans_modelling = NULL;
  turb->les_modelling = NULL;
  turb->reference_values = NULL;

  BFT_FREE(turb);
  *p_turb_struct = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the structure managing the turbulence modelling
 *
 * \param[in, out]  turb   pointer to the structure to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_turbulence_init(cs_cdo_turbulence_t   *turb)
{
  /* The following structures are shared with the Legacy part. So, there is
   * no need to free the structure here */

  turb->model_param = cs_get_glob_turb_model();
  turb->rans_modelling = cs_get_glob_turb_rans_model();
  turb->les_modelling = cs_get_glob_turb_les_model();
  turb->reference_values = cs_get_glob_turb_ref_values();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

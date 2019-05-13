#ifndef __CS_GAS_MIX_H__
#define __CS_GAS_MIX_H__

/*============================================================================
 * Base gas mix data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Gas mix descriptor */
/*-------------------------------*/

typedef struct {

  int           n_species;       /* number of species in the gas mix */
  int          *sp_id_to_f_id;   /* array of species f_ids */

} cs_gas_mix_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main physical constants structure */

extern const cs_gas_mix_t  *cs_glob_gas_mix;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a species field to the gas mix (set of fields).
 *
 * \param[in]   f_id   field id of an already created scalar model field
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_species(int f_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free array mapping gas mix species ids to field ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GAS_MIX_H__ */

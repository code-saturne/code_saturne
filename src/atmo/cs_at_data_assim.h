#ifndef __CS_AT_DATA_ASSIM_H__
#define __CS_AT_DATA_ASSIM_H__

/*============================================================================
 * Data assimilation
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_at_opt_interp.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize data assimilation structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build operators.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_build_ops(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all structures linked to data assimilation.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log data assimilation.
 *
 * \param[in]  ms         pointer to measures set structure
 * \param[in]  oi         pointer to optimal interpolation structure
 * \param[in]  f          pointer to a field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_log(cs_measures_set_t   *ms,
                     cs_at_opt_interp_t  *oi,
                     cs_field_t          *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute source terms for a given variable and add them up to source
 *        term arrays.
 *
 * \param[in]  f_id     field id of variable
 * \param[in]  exp_st   array containing the explicit part of the source term
 * \param[in]  imp_st   array containing the implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_at_data_assim_source_term(int        f_id,
                             cs_real_t *exp_st,
                             cs_real_t *imp_st);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_AT_DATA_ASSIM_H__ */

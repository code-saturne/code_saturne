#ifndef __CS_FUNCTION_DEFAULT_H__
#define __CS_FUNCTION_DEFAULT_H__

/*============================================================================
 * Base predefined function objects.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_function.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Predefined function objects */

typedef enum {

  CS_FUNCTION_CELL_RANK_ID,      /*!< cell MPI rank id (integer) */
  CS_FUNCTION_B_FACE_RANK_ID     /*!< boundary face MPI rank id (integer) */

} cs_function_predefined_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define functions based on code_saturne case setup.
 */
/*----------------------------------------------------------------------------*/

void
cs_function_default_define(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access a function whose data values will be computed
 *        using the a predefined evaluation function.
 *
 * \param[in]   location_id  base associated mesh location id

 * \return  pointer to the associated function object in case of success,
 *          or NULL in case of error
 */
/*----------------------------------------------------------------------------*/

cs_function_t *
cs_function_define_mpi_rank_id(cs_mesh_location_type_t  location_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FUNCTION_DEFAULT_H__ */

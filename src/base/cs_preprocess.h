#ifndef __CS_PREPROCESS_H__
#define __CS_PREPROCESS_H__

/*============================================================================
 * Handle successive preprocessing operations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define all mesh preprocessing operations.
 *----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_define(void);

/*----------------------------------------------------------------------------
 * Apply all mesh preprocessing operations.
 *
 * parameters:
 *   halo_type  <->  type of halo (standard or extended)
 *----------------------------------------------------------------------------*/

void
cs_preprocess_mesh(cs_halo_type_t   halo_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply numbering changes to ignore selected boundary faces.
 *
 * \param[in, out]  m         pointer to mesh structure
 * \param[in, out]  mq        pointer to mesh quantities structure
 * \param[in]       n_faces   number of selected faces
 * \param[in]       face_ids  ids of selected faces
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_selected_b_faces_ignore(cs_mesh_t             *m,
                                           cs_mesh_quantities_t  *mq,
                                           cs_lnum_t              n_faces,
                                           const cs_lnum_t        face_ids[]);

/*----------------------------------------------------------------------------
 * Update fortran arrays relative to the global mesh.
 *----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_update_fortran(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map some mesh arrays for use on device.
 *
 * More elements may be mapped dependin on which arrays are used in
 * accelerated algorithms.
 *
 * \param[in]  alloc_mode  chosen allocation mode
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_update_device(cs_alloc_mode_t  alloc_mode);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PREPROCESS_H__ */

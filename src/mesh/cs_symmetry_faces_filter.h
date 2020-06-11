#ifndef __CS_SYMMETRY_FACES_FILTER_H__
#define __CS_SYMMETRY_FACES_FILTER_H__

/*============================================================================
 * Filter symmetry faces whose effects cancel out
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_symmetry_faces_filter.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_symmetry_faces_filter.c
 * \brief Filter symmetry faces whose effects cancel out
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Filter selected faces whose effects should cancel out.
 *
 * This function simply checks if the sum of associated cell face normals
 * cancels out, and deselects faces for which this is not verified..
 *
 * \param[in]       m          pointer to mesh
 * \param[in]       mq         pointer to mesh quantities
 * \param[in, out]  n_faces    number of selected boundary faces
 * \param[in, out]  face_ids   ids of selected boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_symmetry_faces_filter_cancel(const cs_mesh_t             *m,
                                const cs_mesh_quantities_t  *mq,
                                cs_lnum_t                   *n_faces,
                                cs_lnum_t                    face_ids[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYMMETRY_FACES_FILTER__ */

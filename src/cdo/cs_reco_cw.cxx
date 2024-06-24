/*============================================================================
 * Functions to handle the cell-wise reconstruction of fields relying on the
 * cell mesh structure
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_reco_cw.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct a scalar value at the cell center.
 *        array is assumed to be interlaced and of size stride*n_vertices
 *
 * \param[in]      stride    number of values for each vertex
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      array     array of values
 * \param[in, out] reco      reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_stride_v2c(int                        stride,
                      const cs_cell_mesh_t      *cm,
                      const cs_real_t           *array,
                      cs_real_t                 *reco)
{
  if (array == nullptr)
    return;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));
  assert(reco != nullptr);

  for (int k = 0; k < stride; k++)
    reco[k] = 0;

  for (short int v = 0; v < cm->n_vc; v++) {

    const double  wvc = cm->wvc[v];
    const cs_real_t  *a = array + stride*cm->v_ids[v];

    for (int k = 0; k < stride; k++)
      reco[k] += wvc * a[k];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct a scalar value at the cell center.
 *        array is scanned thanks to the c2v connectivity. Pointer is already
 *        located at the beginning of the cell sequence, i.e. a shift equal to
 *        stride*c2v->idx[cm->c_id] has been done.
 *        array is assumed to be interlaced
 *
 * \param[in]      stride    number of values for each vertex
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      array     array of values for each couple (v,c)
 * \param[in, out] reco      reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_stride_vbyc2c(int                        stride,
                         const cs_cell_mesh_t      *cm,
                         const cs_real_t           *array,
                         cs_real_t                 *reco)
{
  if (array == nullptr)
    return;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));
  assert(reco != nullptr);

  for (int k = 0; k < stride; k++)
    reco[k] = 0;

  for (short int v = 0; v < cm->n_vc; v++) {
    const double  wvc = cm->wvc[v];
    for (int k = 0; k < stride; k++)
      reco[k] += wvc * array[stride*v + k];
  }
}

/*----------------------------------------------------------------------------*/
END_C_DECLS

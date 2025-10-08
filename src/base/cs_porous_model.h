#ifndef __CS_POROUS_MODEL_H__
#define __CS_POROUS_MODEL_H__

/*============================================================================
 * Porous model management
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

#include "fvm/fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

typedef struct {
  fvm_nodal_t *ib_mesh;
  int mesh_id;
  bool activate_post;
} cs_porous_model_extra_faces;

/*============================================================================
 * Global variables
 *============================================================================*/

extern cs_porous_model_extra_faces *cs_glob_porous_model_extra_faces;

/* Choice of the porous model */
extern int cs_glob_porous_model;

/* Specific mesh quantities associated with porous model */
extern cs_mesh_quantities_t  *cs_glob_mesh_quantities_f;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map fluid mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_map_mesh_quantites_f_and_compute(void);

/*----------------------------------------------------------------------------
 * Compute fluid volumes and fluid surfaces in addition to cell volumes
 * and surfaces.
 *
 * parameters:
 *   porous_model <-- porous model option (> 0 for porosity)
 *----------------------------------------------------------------------------*/

void
cs_porous_model_set_model(int  porous_model);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize disable_flag
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_init_disable_flag(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set (unset) has_disable_flag
 *
 * \param[in]  flag   1: on, 0: off
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_set_has_disable_flag(int  flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Init fluid quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_init_fluid_quantities(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute solid quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_mesh_quantities_update(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic computation of the face porosity and factors.
 *
 * This is useful for the integral porous model.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_auto_face_porosity(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Preprocess the fluid surfaces.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_fluid_surfaces_preprocessing(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Penalize porosity and fluid surfaces.
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_clip(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Convert cell array to boundary array
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_convert_cell_to_boundary(const cs_lnum_t   n_ib_cells,
                                         const cs_lnum_t   ibcell_cells[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize porous model arrays
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_postprocess_meshes(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-processes the immersed boundary (ib) planes for display
 *         on paraview.
 *
 * \param[in]       n_ib_cells         ib cell number
 * \param[in]       n_ib_cells_filt    ib cell number after filtering
 * \param[in]       ib_cells_filt      indices of the filtered cells
 * \param[in]       n_glob_vtx         total vertex number
 * \param[in]       ibcell_cells       connectivity ib_cell->cells
 * \param[in]       ibcell_cells_filt  connectivity ib_cell->cells after filtering
 * \param[in]       vtx_ids            vertex ids on both sides of a IB vertex
 *                                     (v0<v1)
 * \param[in]       w_vtx_idx          ib vertex indexes
 * \param[in]       face_vertex_idx    vertex indexes of the ib faces
 * \param[in]       w_vtx              ib vertex coordinates
 */
/*----------------------------------------------------------------------------*/

void
cs_porous_model_post_immmersed_plane(const cs_lnum_t      n_ib_cells,
                                     const cs_lnum_t      n_ib_cells_filt,
                                     const cs_lnum_t      ib_cells_filt[],
                                     const cs_lnum_t      n_glob_vtx,
                                     const cs_lnum_t      ibcell_cells[],
                                     const cs_lnum_t      vtx_ids[][2],
                                     const cs_lnum_t      w_vtx_idx[],
                                     const cs_lnum_t      face_vertex_idx[],
                                     const cs_real_t      w_vtx[][3]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POROUS_MODEL_H__ */

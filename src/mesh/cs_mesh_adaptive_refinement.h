#ifndef __CS_MESH_ADAPTIVE_REFINEMENT_H__
#define __CS_MESH_ADAPTIVE_REFINEMENT_H__

/*============================================================================
 * Adaptive Mesh refinement.
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer for the refinement criteria computation of the adaptive 
 * mesh refinement algorithm.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * parameters:
 *   input <-- pointer to optional (untyped) value or structure.
 *   vals  --> pointer to refinement indicator
 *----------------------------------------------------------------------------*/

typedef void
(cs_amr_indicator_t) (const void  *input,
                      cs_lnum_t   *vals);

/*----------------------------------------------------------------------------
 * Function pointer for field reallocation and interpolation when refining.
 *
 * parameters:
 *   m        <-- pointer to the mesh
 *   location <-- fields location
 *   n_i_elts <-- old number of elements
 *   o2n_idx  <-- old to new element index
 *   new_idx  <-- Newly created elements index (interior faces) or nullptr
 *   o_cog    <-- old mesh cells center of gravity
 *   n_cog    <-- new mesh cells center of gravity
 *----------------------------------------------------------------------------*/

typedef void
(cs_amr_refinement_interpolation_t) (cs_mesh_t       *mesh,
                                     int             location_id,
                                     const cs_lnum_t n_i_elts,
                                     const cs_lnum_t o2n_idx[],
                                     const cs_lnum_t new_idx[],
                                     const cs_real_3_t o_cog[],
                                     const cs_real_3_t n_cog[],
                                     const cs_real_t measure[]);

/*----------------------------------------------------------------------------
 * Function pointer for field reallocation and interpolation when coarsening.
 *
 * parameters:
 *   m        <-- pointer to the mesh
 *   location <-- fields location
 *   n_i_elts <-- old number of elements
 *   n2o_idx  <-- new to old index
 *   n2o      <-- new to old values
 *   measure  <-- elements measure (volume for cells, surface for faces, NULL
                  for vertices
 *----------------------------------------------------------------------------*/

typedef void
(cs_amr_coarsening_interpolation_t) (cs_mesh_t       *mesh,
                                     int             location_id,
                                     const cs_lnum_t n_i_elts,
                                     const cs_lnum_t n2o_idx[],
                                     const cs_lnum_t n2o[],
                                     cs_real_t       measure[]);

/*----------------------------------------------------------------------------
 * Structure containing general information on the adaptive meshing process
 *----------------------------------------------------------------------------*/

typedef struct {

  bool                                is_set;       /* Idicates if AMR 
                                                        is active or not*/
  int                                 n_layers;     /* Number of refine layers
                                                        around user zone*/
  int                                 n_freq;       /* Frequency at which an AMR
                                                       step is performed*/
  cs_amr_indicator_t                  *indic_func;  /* Associated indicator
                                                         computation function*/
  const void                          *indic_input; /* pointer to optional 
                                                         (untyped) value or
                                                         structure */
  cs_amr_refinement_interpolation_t *fields_interp_refinement_func;

  cs_amr_coarsening_interpolation_t *fields_interp_coarsening_func;

  int                                interpolation; /*0, P0 interpolation (default)
                                                      1, simple gradient-based 
                                                      interpolation for fields 
                                                      on cells (ignored for 
                                                      boundaries) */

} cs_amr_info_t ;

/*=============================================================================
 * Static global variables
 *============================================================================*/

extern cs_amr_info_t  *cs_glob_amr_info;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
cs_adaptive_refinement_update_gradients(void);

void
cs_adaptive_refinement_free_gradients(void);

void
cs_adaptive_refinement_define(int                n_layers,
                              int                n_freq,
                              cs_amr_indicator_t *indic_func,
                              const void         *indic_input,
                              int                interpolation);

void
cs_adaptive_refinement_step(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* ___CS_MESH_ADAPTIVE_REFINEMENT_H__ */

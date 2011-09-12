/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_SYR3_COUPLING_H__
#define __CS_SYR3_COUPLING_H__

/*============================================================================
 * Syrthes 3 coupling
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_syr3_comm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

/* Structure associated to Syrthes coupling */

typedef struct _cs_syr3_coupling_t  cs_syr3_coupling_t;

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of SYRTHES couplings.
 *
 * returns:
 *   number of SYRTHES couplings
 *----------------------------------------------------------------------------*/

int
cs_syr3_coupling_n_couplings(void);

/*----------------------------------------------------------------------------
 * Get pointer to SYRTHES coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of SYRTHES coupling
 *
 * returns:
 *   pointer to SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

cs_syr3_coupling_t *
cs_syr3_coupling_by_id(int coupling_id);

/*----------------------------------------------------------------------------
 * Get communicator associated with SYRTHES coupling
 *
 * parameters:
 *   syr_coupling <-- coupling structure associated with SYRTHES
 *
 * returns:
 *   pointer to communicator
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_coupling_get_comm(const cs_syr3_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Get number of vertices in coupled mesh
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_syr3_coupling_get_n_vertices(const cs_syr3_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Get number of associated coupled faces in main mesh
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_syr3_coupling_get_n_faces(const cs_syr3_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Get local list of coupled faces
 *
 * parameters:
 *   syr_coupling    <-- SYRTHES coupling structure
 *   coupl_face_list --> List of coupled faces (1 to n)
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_get_face_list(const cs_syr3_coupling_t  *syr_coupling,
                               cs_lnum_t                  coupl_face_list[]);

/*----------------------------------------------------------------------------
 * Create a syr3_coupling_t structure.
 *
 * parameters:
 *   dim                <-- spatial mesh dimension
 *   ref_axis           <-- reference axis
 *   face_sel_criterion <-- criterion for selection of boundary faces
 *   syr_name           <-- SYRTHES application name
 *   syr_proc_rank      <-- syrthes process rank for MPI
 *   comm_type          <-- communicator type
 *   verbosity          <-- verbosity level
 *   visualization      <-- visualization output level
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_add(int                 dim,
                     int                 ref_axis,
                     const char         *face_sel_criterion,
                     const char         *syr_name,
                     int                 syr_proc_rank,
                     int                 verbosity,
                     int                 visualization);

/*----------------------------------------------------------------------------
 * Initialize communicator for Syrthes coupling
 *
 * parameters:
 *   syr_coupling     <-- SYRTHES coupling structure
 *   syr_id           <-- SYRTHRS coupling id
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_init_comm(cs_syr3_coupling_t  *syr_coupling,
                           int                  syr_id);

/*----------------------------------------------------------------------------
 * Destroy cs_syr3_coupling_t structures
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_all_destroy(void);

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_init_mesh(cs_syr3_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Interpolate a vertex field to an element-centered field
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   vtx_values   <-- values defined on vertices
 *   elt_values   <-> values defined on elements
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_vtx_to_elt(const cs_syr3_coupling_t  *syr_coupling,
                            const cs_real_t           *vtx_values,
                            cs_real_t                 *elt_values);

/*----------------------------------------------------------------------------
 * Interpolate an element-centered field to a vertex field.
 *
 * The size of vtx_values array must be twice the number of vertices.
 * The first half gets values and the second half is used as a working array.
 * The two parts must be contiguous in parallel mode for MPI transfers.
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   elt_values   <-> array of values defined on elements
 *   n_vtx_values <-- number of values defined on vertices
 *   vtx_values   <-> array of values defined on vertices
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_elt_to_vtx(const cs_syr3_coupling_t  *syr_coupling,
                            const cs_real_t           *elt_values,
                            cs_lnum_t                  n_vertices,
                            cs_real_t                 *vtx_values);

/*----------------------------------------------------------------------------
 * Update post-processing variables of a Syrthes coupling
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   step         <-- 0: var = wall temperature
 *                    1: var = fluid temperature
 *                    2: var = exchange coefficient
 *   var          <-- Pointer to variable values
 *----------------------------------------------------------------------------*/

void
cs_syr3_coupling_post_var_update(cs_syr3_coupling_t *syr_coupling,
                                 int                 step,
                                 const cs_real_t    *var);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR3_COUPLING_H__ */

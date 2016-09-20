#ifndef __CS_INTERNAL_COUPLING_H__
#define __CS_INTERNAL_COUPLING_H__

/*============================================================================
 * Internal coupling: coupling for one instance of Code_Saturne
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/


/* Internal coupling structure definition */

typedef struct {

  /* Space dimension */
  int dim;

  /* Selection criterias for coupled domains and coupling interface */
  char *criteria_cells_1;
  char *criteria_cells_2;
  char *criteria_juncture;

  /* local = faces_1, distant = faces_2 */
  ple_locator_t* locator_1;
  /* local = faces_2, distant = faces_1 */
  ple_locator_t* locator_2;

  cs_lnum_t n_1; /* Number of faces in group 1 */
  cs_lnum_t n_2; /* Number of faces in group 2 */

  cs_lnum_t *faces_1; /* Coupling boundary faces of group 1, numbered 1..n   */
  cs_lnum_t *faces_2; /* Coupling boundary faces of group 2  numbered 1..n   */

  cs_lnum_t n_dist_1; /* Number of faces in dist_loc_1 */
  cs_lnum_t n_dist_2; /* Number of faces in dist_loc_2 */

  cs_lnum_t* dist_loc_1; /* Distant boundary faces associated with locator_1 */
  cs_lnum_t* dist_loc_2; /* Distant boundary faces associated with locator_2 */

  /* face i is coupled in this entity if coupled_faces[i] = true */
  bool *coupled_faces;

  /* hint seen from 1 and 2 respectively */
  cs_real_t *hint_1, *hint_2;
  /* hext seen from 1 and 2 respectively */
  cs_real_t *hext_1, *hext_2;

  /* Geometrical weights around coupling interface */
  cs_real_t *gweight_1, *gweight_2;

  /* IJ vectors */
  cs_real_3_t* ij_1;
  cs_real_3_t* ij_2;

  /* OF vectors  */
  cs_real_3_t* ofij_1;
  cs_real_3_t* ofij_2;

  /* Calculation parameters */
  cs_real_t thetav;
  int       idiff;

  /* Gradient reconstruction */
  cs_real_33_t* cocgb_s_lsq;
  cs_real_33_t* cocgb_s_it;
  cs_real_33_t* cocg_s_it;

  /* User information */
  char* namesca;

} cs_internal_coupling_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize coupling criteria from strings.
 *
 * parameters:
 *   criteria_cells_1  <-- string criteria for the first group of cells
 *   criteria_cells_2  <-- string criteria for the second group of cells
 *   criteria_juncture <-- string criteria for the juncture, which is a
 *                         group of faces
 *   cpl               --> pointer to coupling structure to initialize
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_criteria_initialize(const char   criteria_cells_1[],
                                         const char   criteria_cells_2[],
                                         const char   criteria_juncture[],
                                         cs_internal_coupling_t  *cpl);

ple_locator_t *
cs_internal_coupling_create_locator(cs_internal_coupling_t  *cpl);

/*----------------------------------------------------------------------------
 * Initialize locators using selection criteria.
 *
 * parameters:
 *   cpl <-> pointer to coupling structure to modify
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_locators_initialize(cs_internal_coupling_t  *cpl);

/*----------------------------------------------------------------------------
 * Destruction of all internal coupling related structures.
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_finalize(void);

/*----------------------------------------------------------------------------
 * Return the coupling associated with a given coupling_id.
 *
 * parameters:
 *   coupling_id <-> id associated with a coupling entity
 *----------------------------------------------------------------------------*/

cs_internal_coupling_t *
cs_internal_coupling_by_id(int coupling_id);

/*----------------------------------------------------------------------------
 * Exchange quantities from distant to local
 *
 * parameters:
 *   cpl       <-- pointer to coupling entity
 *   stride    <-- Stride (e.g. 1 for double, 3 for interleaved coordinates)
 *   distant_1 <-- Distant values 1, size coupling->n_dist_1
 *   distant_2 <-- Distant values 2, size coupling->n_dist_2
 *   local_1   --> Local values 1, size coupling->n_1
 *   local_2   --> Local values 2, size coupling->n_2
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_var(const cs_internal_coupling_t  *cpl,
                                  int                            stride,
                                  cs_real_t                      distant_1[],
                                  cs_real_t                      distant_2[],
                                  cs_real_t                      local_1[],
                                  cs_real_t                      local_2[]);

/*----------------------------------------------------------------------------
 * Exchange variable between groups using cell id
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
 *   stride   <-- number of values (non interlaced) by entity
 *   tab      <-- variable exchanged
 *   local_1  --> local data for group 1
 *   local_2  --> local data for group 2
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_cell_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local_1[],
                                         cs_real_t                      local_2[]);

/*----------------------------------------------------------------------------
 * Exchange variable between groups using face id
 *
 * parameters:
 *   cpl     <-- pointer to coupling entity
 *   stride  <-- number of values (non interlaced) by entity
 *   tab     <-- variable exchanged
 *   local_1 --> local data for group 1
 *   local_2 --> local data for group 2
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_face_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local_1[],
                                         cs_real_t                      local_2[]);

/*----------------------------------------------------------------------------
 * Modify LSQ COCG matrix to include internal coupling
 *
 * parameters:
 *   coupling <-- pointer to coupling entity
 *   cocg            <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                           cs_real_33_t                   cocg[]);

/*----------------------------------------------------------------------------
 * Modify iterative COCG matrix to include internal coupling
 *
 * parameters:
 *   coupling <-- pointer to coupling entity
 *   cocg            <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_it_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                          cs_real_33_t                   cocg[]);

/*----------------------------------------------------------------------------
 * Initialize internal coupling related structures.
 *
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize(void);

/*----------------------------------------------------------------------------
 * Compute and exchange ij vectors
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_ij(const cs_internal_coupling_t  *cpl);

/*----------------------------------------------------------------------------
 * Add internal coupling rhs contribution for LSQ gradient calculation
 *
 * parameters:
 *   cpl       <-- pointer to coupling entity
 *   c_weight <-- weighted gradient coefficient variable, or NULL
 *   rhsv     <-> rhs contribution modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_rhs(const cs_internal_coupling_t  *cpl,
                             const cs_real_t                c_weight[],
                             cs_real_4_t                    rhsv[]);

/*----------------------------------------------------------------------------
 * Add internal coupling rhs contribution for iterative gradient calculation
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
 *   c_weight <-- weighted gradient coefficient variable, or NULL
 *   grad     <-- pointer to gradient
 *   pvar     <-- pointer to variable
 *   rhs      <-> pointer to rhs contribution
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_iter_rhs(const cs_internal_coupling_t  *cpl,
                              const cs_real_t                c_weight[],
                              cs_real_3_t          *restrict grad,
                              const cs_real_t                pvar[],
                              cs_real_3_t                    rhs[]);

/*----------------------------------------------------------------------------
 * Modify matrix-vector product in case of internal coupling
 *
 * parameters:
 *   exclude_diag <-- extra diagonal flag
 *   matrix       <-- matrix m in m * x = y
 *   x            <-- vector x in m * x = y
 *   y            <-> vector y in m * x = y
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_spmv_contribution(bool              exclude_diag,
                                       void             *input,
                                       const cs_real_t  *restrict x,
                                       cs_real_t        *restrict y);

/*----------------------------------------------------------------------------
 * Add contribution from coupled faces (internal coupling) to polynomial
 * preconditionning.
 *
 * This function is common to most solvers
 *
 * parameters:
 *   input  <-- input
 *   ad     <-> diagonal part of linear equation matrix
 *----------------------------------------------------------------------------*/

void
cs_matrix_preconditionning_add_coupling_contribution(void       *input,
                                                     cs_real_t  *ad);

/*----------------------------------------------------------------------------
 * Return pointers to coupling components
 *
 * parameters:
 *   cpl             <-- pointer to coupling entity
 *   n1              --> NULL or pointer to component n1
 *   n2              --> NULL or pointer to component n2
 *   fac_1[]         --> NULL or pointer to component fac_1[]
 *   fac_2[]         --> NULL or pointer to component fac_2[]
 *   n_dist_1        --> NULL or pointer to component n_dist_1
 *   n_dist_2        --> NULL or pointer to component n_dist_2
 *   dist_loc_1[]    --> NULL or pointer to component dist_loc_1[]
 *   dist_loc_2[]    --> NULL or pointer to component dist_loc_2[]
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_coupled_faces(const cs_internal_coupling_t  *cpl,
                                   cs_lnum_t                     *n1,
                                   cs_lnum_t                     *n2,
                                   cs_lnum_t                     *fac_1[],
                                   cs_lnum_t                     *fac_2[],
                                   cs_lnum_t                     *n_dist_1,
                                   cs_lnum_t                     *n_dist_2,
                                   cs_lnum_t                     *dist_loc_1[],
                                   cs_lnum_t                     *dist_loc_2[]);

/*----------------------------------------------------------------------------
 * Print informations about the given coupling entity
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_print(const cs_internal_coupling_t  *c);

/*----------------------------------------------------------------------------
 * Print informations about all coupling entities
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_dump(void);

/*----------------------------------------------------------------------------
 * Update components hint_* and hext_* using hbord
 *   in the coupling entity associated with given field_id
 *
 * parameters:
 *   field_id <-- id of the field
 *   hbord    <-- array used to update hint_* and hext_*
 *----------------------------------------------------------------------------*/

void
cs_ic_set_exchcoeff(const int         field_id,
                    const cs_real_t  *hbord);

/*----------------------------------------------------------------------------
 * Define coupling entity using given criterias
 *
 * parameters:
 *   field_id   <-- id of the field
 *   volume_1[] <-- string criteria for the first group of cells
 *   volume_2[] <-- string criteria for the second group of cells
 *   juncture[] <-- string criteria for the juncture, which is a
 *                  group of faces
 *----------------------------------------------------------------------------*/

int
cs_internal_coupling_add_entity(int        field_id,
                                const char volume_1[],
                                const char volume_2[],
                                const char juncture[]);

/*----------------------------------------------------------------------------
 * Add contribution from coupled faces (internal coupling) to initialisation
 * for iterative scalar gradient calculation
 *
 * parameters:
 *   cpl      <-- pointer to coupling entity
 *   c_weight <-- weighted gradient coefficient variable, or NULL
 *   pvar     <-- variable
 *   grad     <-> gradient
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_initial_contribution(const cs_internal_coupling_t  *cpl,
                                          const cs_real_t                c_weight[],
                                          const cs_real_t                pvar[],
                                          cs_real_3_t          *restrict grad);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INTERNAL_COUPLING_H__ */

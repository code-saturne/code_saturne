#ifndef __CS_INTERNAL_COUPLING_H__
#define __CS_INTERNAL_COUPLING_H__

/*============================================================================
 * Internal coupling: coupling for one instance of Code_Saturne
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_defs.h"

#include "cs_base.h"
#include "cs_matrix_assembler.h"
#include "cs_mesh.h"
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

  /* Locator + tag for exchanging variables */
  ple_locator_t   *locator;
  int             *c_tag;

  /* Selection criteria for coupled domains */
  char  *cells_criteria;
  char  *faces_criteria;

  cs_lnum_t  n_local; /* Number of faces */
  cs_lnum_t *faces_local; /* Coupling boundary faces, numbered 0..n-1 */

  cs_lnum_t  n_distant; /* Number of faces in faces_distant */
  cs_lnum_t *faces_distant; /* Distant boundary faces associated with locator */

  /* face i is coupled in this entity if coupled_faces[i] = true */
  bool *coupled_faces;

  /* Geometrical weights around coupling interface */
  cs_real_t *g_weight;

  /* IJ vectors */
  cs_real_3_t *ci_cj_vect;

  /* OF vectors  */
  cs_real_3_t *offset_vect;

  /* Gradient reconstruction */
  cs_real_33_t *cocgb_s_lsq;
  cs_real_33_t *cocg_it;

  /* User information */
  char *namesca;

} cs_internal_coupling_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of defined internal couplings.
 *
 * \return  number of internal couplings
 */
/*----------------------------------------------------------------------------*/

int
cs_internal_coupling_n_couplings(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupling volume using given selection criteria.
 *
 * Then, this volume must be seperated from the rest of the domain with a wall.
 *
 * \param[in, out] mesh            pointer to mesh structure to modify
 * \param[in]      criteria_cells  criteria for the first group of cells
 * \param[in]      criteria_faces  criteria for faces to be joined
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add(cs_mesh_t   *mesh,
                         const char   criteria_cells[],
                         const char   criteria_faces[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupling volume using given criteria. Then, this volume will
 * be separated from the rest of the domain with thin walls.
 *
 * \param[in, out] mesh            pointer to mesh structure to modify
 * \param[in]      criteria_cells  criteria for the first group of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_volume(cs_mesh_t  *mesh,
                                const char criteria_cells[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Impose wall BCs to internal coupled faces if not yet defined.
 *
 *   \param[in, out]     bc_type       face boundary condition type
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_bcs(int         bc_type[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destruction of all internal coupling related structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the coupling associated with a given coupling_id.
 *
 * \param[in]  coupling_id  associated with a coupling entity
 *
 * \return pointer to associated coupling structure
 */
/*----------------------------------------------------------------------------*/

cs_internal_coupling_t *
cs_internal_coupling_by_id(int coupling_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange quantities from distant to local
 * (update local using distant).
 *
 * \param[in]  cpl     pointer to coupling entity
 * \param[in]  stride  stride (e.g. 1 for double, 3 for interleaved coordinates)
 * \param[in]  distant distant values, size coupling->n_distant
 * \param[out] local   local values, size coupling->n_local
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_var(const cs_internal_coupling_t  *cpl,
                                  int                            stride,
                                  cs_real_t                      distant[],
                                  cs_real_t                      local[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange variable between groups using cell id.
 *
 * \param[in]  cpl     pointer to coupling entity
 * \param[in]  stride  number of values (non interlaced) by entity
 * \param[in]  tab     variable exchanged
 * \param[out] local   local data
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_cell_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange variable between groups using face id.
 *
 * \param[in]  cpl     pointer to coupling entity
 * \param[in]  stride  number of values (non interlaced) by entity
 * \param[in]  tab     variable exchanged
 * \param[out] local   local data
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_exchange_by_face_id(const cs_internal_coupling_t  *cpl,
                                         int                            stride,
                                         const cs_real_t                tab[],
                                         cs_real_t                      local[]);

/*----------------------------------------------------------------------------
 * Modify LSQ COCG matrix to include internal coupling
 *
 * parameters:
 *   cpl  <-- pointer to coupling entity
 *   cocg <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                           cs_real_33_t                   cocg[]);

/*----------------------------------------------------------------------------
 * Modify LSQ COCG matrix to include internal coupling
 * when diffusivity is a tensor
 *
 * parameters:
 *   cpl  <-- pointer to coupling entity
 *   c_weight  <-- weigthing coefficients
 *   cocg <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_cocg_weighted(const cs_internal_coupling_t  *cpl,
                                       const cs_real_t               *c_weight,
                                       cs_real_33_t                   cocg[]);

/*----------------------------------------------------------------------------
 * Modify iterative COCG matrix to include internal coupling
 *
 * parameters:
 *   cpl  <-- pointer to coupling entity
 *   cocg <-> cocg matrix modified
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_it_cocg_contribution(const cs_internal_coupling_t  *cpl,
                                          cs_real_33_t                   cocg[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup internal coupling related parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize internal coupling related structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add internal coupling rhs contribution for LSQ gradient calculation
 *
 * \param[in]       cpl       pointer to coupling entity
 * \param[in]       c_weight  weighted gradient coefficient variable, or NULL
 * \param[in]       w_stride  stride of weighting coefficient
 * \param[in, out]  rhsv      pointer to rhs contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_scalar_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    const int                      w_stride,
    cs_real_4_t                    rhsv[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add internal coupling rhs contribution for LSQ gradient calculation
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       c_weight weighted gradient coefficient variable, or NULL
 * \param[in]       w_stride stride of weighting coefficient
 * \param[in]       pvar     pointer to variable
 * \param[in, out]  rhs      pointer to rhs contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_vector_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    const int                      w_stride,
    const cs_real_3_t              pvar[],
    cs_real_33_t                   rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add internal coupling rhs contribution for LSQ gradient calculation
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       c_weight weighted gradient coefficient variable, or NULL
 * \param[in]       w_stride stride of weighting coefficient
 * \param[in]       pvar     pointer to variable
 * \param[in, out]  rhs      pointer to rhs contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_lsq_tensor_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    const int                      w_stride,
    const cs_real_6_t              pvar[],
    cs_real_63_t                   rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add internal coupling rhs contribution for iterative gradient
 * calculation
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       c_weight weighted gradient coefficient variable, or NULL
 * \param[in]       grad     pointer to gradient
 * \param[in]       pvar     pointer to variable
 * \param[in, out]  rhs      pointer to rhs contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_iterative_scalar_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    cs_real_3_t          *restrict grad,
    const cs_real_t                pvar[],
    cs_real_3_t                    rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add internal coupling rhs contribution for iterative vector gradient
 * calculation
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       c_weight weighted gradient coefficient variable, or NULL
 * \param[in]       grad     pointer to gradient
 * \param[in]       pvar     pointer to variable
 * \param[in, out]  rhs      pointer to rhs contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_iterative_vector_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    cs_real_33_t         *restrict grad,
    const cs_real_3_t              pvar[],
    cs_real_33_t                   rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add internal coupling rhs contribution for iterative tensor gradient
 * calculation
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       c_weight weighted gradient coefficient variable, or NULL
 * \param[in]       grad     pointer to gradient
 * \param[in]       pvar     pointer to variable
 * \param[in, out]  rhs      pointer to rhs contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_iterative_tensor_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    cs_real_63_t         *restrict grad,
    const cs_real_6_t              pvar[],
    cs_real_63_t                   rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add internal coupling contribution for reconstruction of the
 * gradient of a scalar.
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       r_grad   pointer to reconstruction gradient
 * \param[in, out]  grad     pointer to gradient to be reconstructed var
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_reconstruct_scalar_gradient(
    const cs_internal_coupling_t  *cpl,
    cs_real_3_t          *restrict r_grad,
    cs_real_3_t                    grad[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add internal coupling contribution for reconstruction of the
 * gradient of a vector.
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       r_grad   pointer to reconstruction gradient
 * \param[in, out]  grad     pointer to gradient to be reconstructed var
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_reconstruct_vector_gradient(
    const cs_internal_coupling_t  *cpl,
    cs_real_33_t         *restrict r_grad,
    cs_real_33_t                   grad[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add internal coupling contribution for reconstruction of the
 *         gradient of a symmetric tensor.
 *
 * \param[in]       cpl      pointer to coupling entity
 * \param[in]       r_grad   pointer to reconstruction gradient
 * \param[in, out]  grad     pointer to gradient to be reconstructed var
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_reconstruct_tensor_gradient(
    const cs_internal_coupling_t  *cpl,
    cs_real_63_t         *restrict r_grad,
    cs_real_63_t                   grad[]);

/*----------------------------------------------------------------------------
 * Addition to matrix-vector product in case of internal coupling.
 *
 * parameters:
 *   exclude_diag <-- extra diagonal flag
 *   f            <-- associated field pointer
 *   x            <-- vector x in m * x = y
 *   y            <-> vector y in m * x = y
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_spmv_contribution(bool               exclude_diag,
                                       const cs_field_t  *f,
                                       const cs_real_t   *restrict x,
                                       cs_real_t         *restrict y);

/*----------------------------------------------------------------------------
 * Add coupling term coordinates to matrix assembler.
 *
 * parameters:
 *   coupling_id
 *   r_g_id   <-- global row ids (per cell)
 *   ma       <-> matrix assembler
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_matrix_add_ids(int                     coupling_id,
                                    const cs_gnum_t        *r_g_id,
                                    cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------
 * Add coupling terms to matrix values assembly.
 *
 * parameters:
 *   f        <-- associated field
 *   db_size  <-- diagonal block size
 *   eb_size  <-- extra-diagonal block size
 *   r_g_id   <-- global row ids (per cell)
 *   mav      <-> matrix values assembler
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_matrix_add_values(const cs_field_t              *f,
                                       cs_lnum_t                      db_size,
                                       cs_lnum_t                      eb_size,
                                       const cs_gnum_t                r_g_id[],
                                       cs_matrix_assembler_values_t  *mav);

/*----------------------------------------------------------------------------
 * Return pointers to coupling components
 *
 * parameters:
 *   cpl             <-- pointer to coupling entity
 *   n_local         --> NULL or pointer to component n_local
 *   faces_local     --> NULL or pointer to component faces_local
 *   n_distant       --> NULL or pointer to component n_distant
 *   faces_distant   --> NULL or pointer to component faces_distant
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_coupled_faces(const cs_internal_coupling_t  *cpl,
                                   cs_lnum_t                     *n_local,
                                   cs_lnum_t                     *faces_local[],
                                   cs_lnum_t                     *n_distant,
                                   cs_lnum_t                     *faces_distant[]);

/*----------------------------------------------------------------------------
 * Log information about a given internal coupling entity
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_log(const cs_internal_coupling_t  *cpl);

/*----------------------------------------------------------------------------
 * Print informations about all coupling entities
 *
 * parameters:
 *   cpl <-- pointer to coupling entity
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_dump(void);

/*----------------------------------------------------------------------------
 * Add preprocessing operations required by coupling volume using given
 * criteria.
 *
 * The volume is seperated from the rest of the domain with inserted
 * boundaries.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure to modify
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_preprocess(cs_mesh_t   *mesh);

/*----------------------------------------------------------------------------
 * Define face to face mappings for internal couplings.
 *
 * parameters:
 *   mesh           <-> pointer to mesh structure to modify
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_map(cs_mesh_t   *mesh);

/*----------------------------------------------------------------------------
 * Define coupling entity using given criteria.
 *
 * parameters:
 *   f_id       <-- id of the field
 *----------------------------------------------------------------------------*/

void
cs_internal_coupling_add_entity(int        f_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add contribution from coupled faces (internal coupling) to
 * initialisation for iterative scalar gradient calculation.
 *
 * \param[in]      cpl       pointer to coupling entity
 * \param[in]      c_weight  weighted gradient coefficient variable, or NULL
 * \param[in]      pvar      variable
 * \param[in, out] grad      gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize_scalar_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    const cs_real_t                pvar[],
    cs_real_3_t          *restrict grad);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add contribution from coupled faces (internal coupling) to
 * initialisation for iterative vector gradient calculation
 *
 * \param[in]       cpl       pointer to coupling entity
 * \param[in]       c_weight  weighted gradient coefficient variable, or NULL
 * \param[in]       pvar      variable
 * \param[in, out]  grad      gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize_vector_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    const cs_real_3_t              pvar[],
    cs_real_33_t         *restrict grad);


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add contribution from coupled faces (internal coupling) to
 * initialisation for iterative symmetric tensor gradient calculation
 *
 * \param[in]       cpl       pointer to coupling entity
 * \param[in]       c_weight  weighted gradient coefficient variable, or NULL
 * \param[in, out]  pvar      variable
 * \param[in, out]  grad      gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_internal_coupling_initialize_tensor_gradient(
    const cs_internal_coupling_t  *cpl,
    const cs_real_t                c_weight[],
    const cs_real_6_t              pvar[],
    cs_real_63_t         *restrict grad);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update internal coupling coefficients of the field of the
 * given id using given boundary exchange coefficients passed by face id.
 *
 * \param[in] field_id  field id
 * \param[in] hbnd      boundary exchange coefficients passed by face id
 */
/*----------------------------------------------------------------------------*/

void
cs_ic_field_set_exchcoeff(const int         field_id,
                          const cs_real_t  *hbnd);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get distant data using face id at all coupling faces for a given
 * field id.
 *
 * \param[in]  field_id    field id
 * \param[in]  stride      number of values (interlaced) by entity
 * \param[in]  tab_distant exchanged data by face id
 * \param[out] tab_local   local data by face id
 */
/*----------------------------------------------------------------------------*/

void
cs_ic_field_dist_data_by_face_id(const int         field_id,
                                 int               stride,
                                 const cs_real_t   tab_distant[],
                                 cs_real_t         tab_local[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INTERNAL_COUPLING_H__ */

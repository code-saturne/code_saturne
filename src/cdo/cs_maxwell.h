#ifndef __CS_MAXWELL_H__
#define __CS_MAXWELL_H__

/*============================================================================
 * Header to handle the maxwell module with CDO schemes
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_time_step.h"
#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_maxwell.h

  \brief Structure and routines handling the Maxwell module dedicated to
         the resolution of electro-magnetic equations

*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Generic name given to fields and equations related to this module */

#define CS_MAXWELL_ESTATIC_EQNAME   "electrostatic"
#define CS_MAXWELL_EFIELD_NAME      "electric_field"
#define CS_MAXWELL_DFIELD_NAME      "electric_induction"

#define CS_MAXWELL_MSTATIC_EQNAME   "magnetostatic"
#define CS_MAXWELL_MFIELD_NAME      "magnetic_field"
#define CS_MAXWELL_BFIELD_NAME      "magnetic_induction"

#define CS_MAXWELL_JEFFECT_NAME     "joule_effect"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * @name Flags specifying the modelling for the Maxwell module
 * @{
 *
 * These fields are considered:
 * H (A/m) the magnetic field
 * E (V/m) the electric field
 * B (Vs/m^2) the magnetic induction or magnetic flux density
 * D (As/m^2) the electric induction or electric flux density
 *
 * Moreover, one considers
 * rho_e the electric charge density
 * j the electric source current (or current for short)
 *
 * One assumes a quasi-static approximation of the maxwell equations
 * i.e. ones solves the Ampere equation: curl(H) = j
 * and one solves the Faraday equation: curl(E) = -partial_t B
 *
 * \def CS_MAXWELL_MODEL_ELECTROSTATIC
 * \brief Solve the equation -div(epsilon grad(phi)) = rho_e + BCs which
 *        results from div(D) = rho_e, D = epsilon E and E = -grad(phi)
 *
 * \def CS_MAXWELL_MODEL_MAGNETOSTATIC
 * \brief Solve the system curl(H) = j, div(mu H) = 0 + BCs which yields
 *        curl( 1/mu curl(A)) = j and div(A) = 0 (Coulomb gauge) + BCs
 *        if one sets: mu.H = curl(A)
 */

#define CS_MAXWELL_MODEL_ELECTROSTATIC      (1 << 0) /* 1 */
#define CS_MAXWELL_MODEL_MAGNETOSTATIC      (1 << 1) /* 2 */

/*!
 * @}
 * @name Flags specifying options for the Maxwell module
 * @{
 *
 * \def CS_MAXWELL_A_PHI_FORMULATION
 * \brief Vector-valued potential (A) and scalar-valued potential (phi) are
 *        introduced to solve the Maxwell equations. One introduces the
 *        following identities: curl(A) = B and E = -partial_t A - grad(phi)
 *
 * \def CS_MAXWELL_JOULE_EFFECT
 * \brief Take into account the Joule effect
 */

#define CS_MAXWELL_A_PHI_FORMULATION       (1 << 0) /* 1 */
#define CS_MAXWELL_JOULE_EFFECT            (1 << 1) /* 2 */

/*!
 * @}
 */

typedef struct _maxwell_t  cs_maxwell_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if the computation of Maxwell equations is activated
 */
/*----------------------------------------------------------------------------*/

bool
cs_maxwell_is_activated(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the future computation of the Maxwell equations
 *
 * \param[in]      model         type of modelling
 * \param[in]      options       flag to handle optional parameters
 *
 * \return a pointer to a new allocated Maxwell structure
 */
/*----------------------------------------------------------------------------*/

cs_maxwell_t *
cs_maxwell_activate(cs_flag_t    model,
                    cs_flag_t    options);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the Maxwell module
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_maxwell_t *
cs_maxwell_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup equations/properties related to the Maxwell module
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_init_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize the setup stage for equations related to the Maxwell module
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_finalize_setup(const cs_cdo_connect_t       *connect,
                          const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a summary of the Maxwell module
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve if needed the steady-state equations related to the Maxwell
 *         module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_compute_steady_state(const cs_mesh_t              *mesh,
                                const cs_time_step_t         *time_step,
                                const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve equations related to the Maxwell module
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_compute(const cs_mesh_t              *mesh,
                   const cs_time_step_t         *time_step,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update/initialize the Maxwell module according to the settings
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts         pointer to a cs_time_step_t structure
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_update(const cs_mesh_t             *mesh,
                  const cs_cdo_connect_t      *connect,
                  const cs_cdo_quantities_t   *quant,
                  const cs_time_step_t        *ts,
                  bool                         cur2prev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the Maxwell module
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_extra_op(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the Maxwell module.
 *         Prototype of this function is fixed since it is a function pointer
 *         defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_maxwell_extra_post(void                      *input,
                      int                        mesh_id,
                      int                        cat_id,
                      int                        ent_flag[5],
                      cs_lnum_t                  n_cells,
                      cs_lnum_t                  n_i_faces,
                      cs_lnum_t                  n_b_faces,
                      const cs_lnum_t            cell_ids[],
                      const cs_lnum_t            i_face_ids[],
                      const cs_lnum_t            b_face_ids[],
                      const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MAXWELL_H__ */

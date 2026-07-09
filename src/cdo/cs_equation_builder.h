#ifndef CS_EQUATION_BUILDER_H
#define CS_EQUATION_BUILDER_H

/*============================================================================
 * Functions to handle the equation builder structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "alge/cs_matrix.h"
#include "cdo/cs_cdo_bc.h"
#include "cdo/cs_cdo_connect.h"
#include "cdo/cs_cdo_local.h"
#include "cdo/cs_cdo_quantities.h"
#include "cdo/cs_cdo_system.h"
#include "cdo/cs_enforcement.h"
#include "cdo/cs_equation_param.h"
#include "cdo/cs_flag.h"
#include "cdo/cs_source_term.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _equation_builder_t  cs_equation_builder_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic function prototype for a hook during the cellwise building
 *        of the linear system
 *        Enable an advanced user to get a fine control of the discretization
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context to cast for this discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] context     pointer to a context structure
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure (mass matrix)
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure (diffusion)
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_hook_t)
(
 [[maybe_unused]] const cs_equation_param_t   *eqp,
 [[maybe_unused]] const cs_equation_builder_t *eqb,
 [[maybe_unused]] const void                  *eqc,
 [[maybe_unused]] const cs_cell_mesh_t        *cm,
 [[maybe_unused]] void                        *context,
 [[maybe_unused]] cs_hodge_t                  *mass_hodge,
 [[maybe_unused]] cs_hodge_t                  *diff_hodge,
 [[maybe_unused]] cs_cell_sys_t               *csys,
 [[maybe_unused]] cs_cell_builder_t           *cb
);

/*! \struct cs_equation_builder_t
 *  \brief Store common elements used when building an algebraic system
 *  related to an equation
 */

struct _equation_builder_t {

  bool        init_step;  /*!< true if this is the initialization step */

  /*!
   * @name Flags to know what to build and how to build such terms
   * @{
   */

  cs_eflag_t  msh_flag;   /*!< Flag storing which quantities to build in a
                           *   \ref cs_cell_mesh_t structure for all cells */
  cs_eflag_t  bdy_flag;   /*!< Flag storing which quantities to build in a
                           *   \ref cs_cell_mesh_t structure for a boundary
                           *   cell */
  cs_eflag_t  src_flag;   /*!< Flag storing which quantities to build in a
                           *   \ref cs_cell_mesh_t structure for the specific
                           *   computation of source terms */
  cs_flag_t   sys_flag;   /*!< Metadata related to the sytem */

  /*!
   * @}
   * @name Metadata related to associated physical properties
   * @{
   */

  bool   diff_pty_uniform;      /*!< Is diffusion property uniform ? */
  bool   curlcurl_pty_uniform;  /*!< Is curl-curl property uniform ? */
  bool   graddiv_pty_uniform;   /*!< Is grad-div property uniform ? */
  bool   time_pty_uniform;      /*!< Is time property uniform ? */
  bool   reac_pty_uniform[CS_CDO_N_MAX_REACTIONS]; /*!< Is each reaction
                                                    * property uniform ? */

  /*!
   * @}
   * @name Source terms
   * @{
   */

  cs_mask_t   *source_mask;  /*!< null if no source term or one source term
                              * is defined. Allocated to n_cells in order to
                              * know in each cell which source term has to be
                              * computed */

  /*! \var compute_source
   * Pointer to functions which compute the value of the source term
   */

  cs_source_term_cellwise_t  **compute_source;

  /*!
   * @}
   * @name Helper structure to build the matrix and manage arrays of DoFs
   * @{
   */

  /*! \var  system_helper
   *  Pointer to the associated system helper structure
   */

  cs_cdo_system_helper_t     *system_helper;

  /*!
   * @}
   * @name Enforcement of degrees of freedom (DoFs)
   * @{
   */

  /*! \var  enforced_values
   *  Array of values to enforced
   */

  cs_real_t                  *enforced_values;

  /*!
   * @}
   * @name Incremental solving
   * @{
   *
   * \var increment
   * array of values for the last computed increment. Only allocated if an
   * incremental solving has been requested.
   *
   * \var incremental_algo
   * Structure which handles the incremental algorithm
   */

  cs_real_t              *increment;
  cs_iter_algo_t         *incremental_algo;

  /*!
   * @}
   * @name Boundary conditions
   * @{
   *
   * \var face_bc
   * face_bc should not change during the simulation.
   * The case of a definition of the BCs which changes of type during the
   * simulation is possible but not implemented.
   * You just have to call the initialization step each time the type of BCs
   * is modified to define an updated \ref cs_cdo_bc_face_t structure.
   */

  cs_cdo_bc_face_t   *face_bc;  /*!< Information about boundary conditions
                                     applied to faces */

  cs_real_t          *dir_values; /*!< Array storing the Dirichlet values at
                                     DoFs */

  /*!
   * @}
   * @name User hook
   * @{
   *
   * \var hook_context
   * Pointer to a shared structure (the lifecycle of this structure is not
   * managed by the current cs_equation_builder_t structure)
   *
   * \var hook_function
   * Function pointer associated to a predefined prototype
   *
   * This function allows one to modify the cellwise system (matrix and rhs)
   * before applying the time scheme and the static condensation (if needed)
   * and the strong/penalized enforcement of boundary conditions.
   *
   * This is useful to add a term in the equation like an advanced source term
   * without the need to allocate an array and with an access to the local
   * structure such as the local cell mesh, the cell builder and high-level
   * structures related to an equation
   */

  void                        *hook_context;
  cs_equation_build_hook_t    *hook_function;

  /*!
   * @}
   * @name Performance monitoring
   * @{
   *
   * Monitoring the efficiency of the algorithm used to manipulate/build
   * an equation.
   */

  cs_timer_counter_t     tcb; /*!< Cumulated elapsed time for building the
                               *   current system */
  cs_timer_counter_t     tcs; /*!< Cumulated elapsed time for solving the
                               *   current system */
  cs_timer_counter_t     tce; /*!< Cumulated elapsed time for computing
                               *   all extra operations (post, balance,
                               *   fluxes...) */

  /*! @} */
};

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the flag to give for building a cs_cell_mesh_t structure
 *
 * \param[in] cell_flag  flag related to the current cell
 * \param[in] eqb        pointer to a cs_equation_builder_t structure
 *
 * \return the flag to set for the current cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_eflag_t
cs_equation_builder_cell_mesh_flag(cs_flag_t                     cell_flag,
                                   const cs_equation_builder_t  *eqb)
{
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)
    return eqb->msh_flag | eqb->src_flag | eqb->bdy_flag;
  else
    return eqb->msh_flag | eqb->src_flag;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

void
cs_equation_builder_update_default_flags(cs_eflag_t msh_flag,
                                         cs_eflag_t bdy_flag,
                                         cs_eflag_t src_flag);

cs_equation_builder_t *
cs_equation_builder_create(const cs_equation_param_t *eqp,
                           const cs_mesh_t           *mesh);

const cs_matrix_t *
cs_equation_builder_get_matrix(const cs_equation_builder_t *builder,
                               int                          block_id);

const cs_range_set_t *
cs_equation_builder_get_range_set(const cs_equation_builder_t *builder,
                                  int                          block_id);

void
cs_equation_builder_free(cs_equation_builder_t  **p_builder);

void
cs_equation_builder_reset(cs_equation_builder_t  *eqb);

void
cs_equation_builder_apply_default_flags(cs_equation_builder_t  *eqb);

void
cs_equation_builder_log_performance(size_t                       compute_amount,
                                    const cs_equation_param_t   *eqp,
                                    const cs_equation_builder_t *eqb);

bool
cs_equation_builder_set_reaction_pty_cw(const cs_equation_param_t     *eqp,
                                        const cs_equation_builder_t   *eqb,
                                        const cs_cell_mesh_t          *cm,
                                        cs_cell_builder_t             *cb);

void
cs_equation_builder_init_properties(const cs_equation_param_t     *eqp,
                                    const cs_equation_builder_t   *eqb,
                                    cs_hodge_t                    *diff_hodge,
                                    cs_cell_builder_t             *cb);

void
cs_equation_builder_enforce_dofs(const cs_equation_builder_t *eqb,
                                 cs_cell_builder_t           *cb,
                                 cs_cell_sys_t               *csys);

void
cs_equation_builder_enforce_block_dofs(const cs_equation_builder_t *eqb,
                                       cs_cell_builder_t           *cb,
                                       cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/

#endif /* CS_EQUATION_BUILDER_H */

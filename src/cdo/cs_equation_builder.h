#ifndef __CS_EQUATION_BUILDER_H__
#define __CS_EQUATION_BUILDER_H__

/*============================================================================
 * Functions to handle the equation builder structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_system.h"
#include "cs_enforcement.h"
#include "cs_equation_param.h"
#include "cs_flag.h"
#include "cs_matrix.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _equation_builder_t  cs_equation_builder_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Generic function prototype for a hook during the cellwise building
 *          of the linear system
 *          Enable an advanced user to get a fine control of the discretization
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
(cs_equation_build_hook_t)(const cs_equation_param_t     *eqp,
                           const cs_equation_builder_t   *eqb,
                           const void                    *eqc,
                           const cs_cell_mesh_t          *cm,
                           void                          *context,
                           cs_hodge_t                    *mass_hodge,
                           cs_hodge_t                    *diff_hodge,
                           cs_cell_sys_t                 *csys,
                           cs_cell_builder_t             *cb);

/*! \struct cs_equation_builder_t
 *  \brief Store common elements used when building an algebraic system
 *  related to an equation
 */

struct _equation_builder_t {

  bool         init_step;    /*!< true if this is the initialization step */

  /*!
   * @name Flags to know what to build and how to build such terms
   * @{
   */

  cs_eflag_t   msh_flag;     /*!< Information related to what to build in a
                              *   \ref cs_cell_mesh_t structure for a generic
                              *   cell */
  cs_eflag_t   bd_msh_flag;  /*!< Information related to what to build in a
                              *   \ref cs_cell_mesh_t structure for a cell close
                              *   to the boundary */
  cs_eflag_t   st_msh_flag;  /*!< Information related to what to build in a
                              *   \ref cs_cell_mesh_t structure when only the
                              *   source term has to be built */
  cs_flag_t    sys_flag;     /*!< Information related to the sytem */

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

  cs_mask_t   *source_mask;  /*!< NULL if no source term or one source term
                              * is defined. Allocated to n_cells in order to
                              * know in each cell which source term has to be
                              * computed */

  /*! \var compute_source
   * Pointer to functions which compute the value of the source term
   */

  cs_source_term_cellwise_t  *compute_source[CS_N_MAX_SOURCE_TERMS];

  /*!
   * @}
   * @name Helper structure to build the matrix and manage arrays of DoFs
   * @{
   */

  cs_cdo_system_helper_t     *system_helper;

  /*! \var  */
  /*!
   * @}
   * @name Enforcement of degrees of freedom (DoFs)
   * @{
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
 * Inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the flag to give for building a cs_cell_mesh_t structure
 *
 * \param[in]  cell_flag   flag related to the current cell
 * \param[in]  eqb         pointer to a cs_equation_builder_t structure
 *
 * \return the flag to set for the current cell
 */
/*----------------------------------------------------------------------------*/

static inline cs_eflag_t
cs_equation_builder_cell_mesh_flag(cs_flag_t                      cell_flag,
                                   const cs_equation_builder_t   *eqb)
{
  cs_eflag_t  _flag = eqb->msh_flag | eqb->st_msh_flag;

  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)
    _flag |= eqb->bd_msh_flag;

  return _flag;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a new structure to handle the building of algebraic system
 *         related to a cs_equation_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_builder_t *
cs_equation_builder_create(const cs_equation_param_t   *eqp,
                           const cs_mesh_t             *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the range set structure associated to a builder structure
 *        for the block defined in block_id in the system helper structure
 *
 * \param[in, out]  builder      pointer to a cs_equation_builder_t
 * \param[in]       block_id     id of the block to consider
 *
 * \return a pointer to a cs_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_t *
cs_equation_builder_get_matrix(const cs_equation_builder_t  *builder,
                               int                           block_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the range set structure associated to a builder structure
 *        for the block defined in block_id in the system helper structure
 *
 * \param[in, out]  builder      pointer to a cs_equation_builder_t
 * \param[in]       block_id     id of the block to consider
 *
 * \return a pointer to a cs_range_set structure
 */
/*----------------------------------------------------------------------------*/

const cs_range_set_t *
cs_equation_builder_get_range_set(const cs_equation_builder_t  *builder,
                                  int                           block_id);
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_builder_t structure
 *
 * \param[in, out]  p_builder  pointer of pointer to the cs_equation_builder_t
 *                             structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_free(cs_equation_builder_t  **p_builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free some members of a cs_equation_builder_t structure
 *
 * \param[in, out]  eqb   pointer to the cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_reset(cs_equation_builder_t  *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a message in the performance output file related to the
 *          monitoring of equation
 *
 * \param[in]  eqp    pointer to a set of equation parameters
 * \param[in]  eqb    pointer to an equation builder  structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_log_performance(const cs_equation_param_t     *eqp,
                                    const cs_equation_builder_t   *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize all reaction properties.
 *        This function is shared across all CDO schemes. The cs_cell_builder_t
 *        structure stores the computed property values. If the property is
 *        uniform, a first call to the function
 *        cs_equation_builder_init_properties has to be done before the loop on
 *        cells
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      eqb      pointer to a cs_equation_builder_t structure
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out] cb       pointer to a \ref cs_cell_builder_t structure
 *
 * \return true if the reaction property is not equal to zero
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_builder_set_reaction_pty_cw(const cs_equation_param_t     *eqp,
                                        const cs_equation_builder_t   *eqb,
                                        const cs_cell_mesh_t          *cm,
                                        cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize all properties potentially useful to build the algebraic
 *         system. This function is shared across all CDO schemes.
 *         The \ref cs_cell_builder_t structure stores property values related
 *         to the reaction term, unsteady term and grad-div term.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in, out] diff_hodge   pointer to the diffusion hodge structure
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_init_properties(const cs_equation_param_t     *eqp,
                                    const cs_equation_builder_t   *eqb,
                                    cs_hodge_t                    *diff_hodge,
                                    cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account the enforcement of internal DoFs. Apply an
 *          algebraic manipulation. Update members of the cs_cell_sys_t
 *          structure related to the internal enforcement.
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_enforce_dofs(const cs_equation_builder_t     *eqb,
                                 cs_cell_builder_t               *cb,
                                 cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take into account the enforcement of internal DoFs. Case of matrices
 *        defined by blocks. Apply an algebraic manipulation. Update members
 *        of the cs_cell_sys_t structure related to the internal enforcement.
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_enforce_block_dofs(const cs_equation_builder_t   *eqb,
                                       cs_cell_builder_t             *cb,
                                       cs_cell_sys_t                 *csys);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_BUILDER_H__ */

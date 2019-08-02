#ifndef __CS_EQUATION_COMMON_H__
#define __CS_EQUATION_COMMON_H__

/*============================================================================
 * Routines to handle common equation features for building algebraic system
 * in CDO schemes
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

#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_param.h"
#include "cs_flag.h"
#include "cs_matrix.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Strategy of synchronization of values shared across several cells
 * This applies to vertices, edges and faces
 *
 * CS_EQUATION_SYNC_ZERO_VALUE
 * If zero is a possible value then set this value, otherwise one takes
 * the mean-value
 *
 * CS_EQUATION_SYNC_MEAN_VALUE
 * Compute the mean-value across values to set
 *
 */

#define  CS_EQUATION_SYNC_ZERO_VALUE    1
#define  CS_EQUATION_SYNC_MEAN_VALUE    2

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_equation_builder_t
 *  \brief Store common elements used when building an algebraic system
 *  related to an equation
 */

typedef struct {

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

  /*!
   * @}
   * @name Performance monitoring
   * @{
   *
   * Monitoring the efficiency of the algorithm used to manipulate/build
   * an equation.
   */

  cs_timer_counter_t     tcb; /*!< Cumulated elapsed time for building the
                               *   current system: tcb >= tcd+tca+tcr+tcs+tcs */
  cs_timer_counter_t     tcd; /*!< Cumulated elapsed time for building
                               *   diffusion terms */
  cs_timer_counter_t     tca; /*!< Cumulated elapsed time for building
                               *   advection terms */
  cs_timer_counter_t     tcr; /*!< Cumulated elapsed time for building
                               *   reaction terms */
  cs_timer_counter_t     tcs; /*!< Cumulated elapsed time for building
                               *   source terms */
  cs_timer_counter_t     tce; /*!< Cumulated elapsed time for computing
                               *   all extra operations (post, balance,
                               *   fluxes...) */

  /*! @} */

} cs_equation_builder_t;

/*
 * Structure used to store information generated during the analysis
 * of the balance of each term of an equation
 */
typedef struct {

  /* where balance is computed: primal vertices or primal cells */
  cs_flag_t       location;
  cs_lnum_t       size;
  cs_real_t      *balance;

  /* Balance for each main term */
  cs_real_t      *unsteady_term;
  cs_real_t      *reaction_term;
  cs_real_t      *diffusion_term;
  cs_real_t      *advection_term;
  cs_real_t      *source_term;
  cs_real_t      *boundary_term;

} cs_equation_balance_t;

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
cs_equation_cell_mesh_flag(cs_flag_t                      cell_flag,
                           const cs_equation_builder_t   *eqb)
{
  cs_eflag_t  _flag = eqb->msh_flag | eqb->st_msh_flag;

  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)
    _flag |= eqb->bd_msh_flag;

  return _flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the diffusion property inside a cell and its related quantities
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      c_id    id of the cell to deal with
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_equation_set_diffusion_property(const cs_equation_param_t  *eqp,
                                   const cs_lnum_t             c_id,
                                   const cs_real_t             t_eval,
                                   const cs_flag_t             c_flag,
                                   cs_cell_builder_t          *cb)
{
  cs_property_get_cell_tensor(c_id,
                              t_eval,
                              eqp->diffusion_property,
                              eqp->diffusion_hodge.inv_pty,
                              cb->dpty_mat);

  if (cs_property_is_isotropic(eqp->diffusion_property))
    cb->dpty_val = cb->dpty_mat[0][0];

  /* Set additional quantities in case of more advanced way of enforcing the
     essential BCs */
  if (c_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {
    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
      cs_math_33_eigen((const cs_real_t (*)[3])cb->dpty_mat,
                       &(cb->eig_ratio),
                       &(cb->eig_max));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the curl-curl property inside a cell and its related quantities
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      c_id    id of the cell to deal with
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_equation_set_curlcurl_property(const cs_equation_param_t  *eqp,
                                  const cs_lnum_t             c_id,
                                  const cs_real_t             t_eval,
                                  const cs_flag_t             c_flag,
                                  cs_cell_builder_t          *cb)
{
  cs_property_get_cell_tensor(c_id,
                              t_eval,
                              eqp->curlcurl_property,
                              eqp->curlcurl_hodge.inv_pty,
                              cb->cpty_mat);

  if (cs_property_is_isotropic(eqp->curlcurl_property))
    cb->cpty_val = cb->cpty_mat[0][0];

  /* Set additional quantities in case of more advanced way of enforcing
     essential BCs */
  if (c_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {
    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
      cs_math_33_eigen((const cs_real_t (*)[3])cb->cpty_mat,
                       &(cb->eig_ratio),
                       &(cb->eig_max));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the diffusion property inside a cell and its related quantities.
 *         Cellwise version using a cs_cell_mesh_t structure
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_equation_set_diffusion_property_cw(const cs_equation_param_t   *eqp,
                                      const cs_cell_mesh_t        *cm,
                                      const cs_real_t              t_eval,
                                      const cs_flag_t              c_flag,
                                      cs_cell_builder_t           *cb)
{
  cs_property_tensor_in_cell(cm,
                             eqp->diffusion_property,
                             t_eval,
                             eqp->diffusion_hodge.inv_pty,
                             cb->dpty_mat);

  if (cs_property_is_isotropic(eqp->diffusion_property))
    cb->dpty_val = cb->dpty_mat[0][0];

  /* Set additional quantities in case of more advanced way of enforcing the
     essential BCs */
  if (c_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {
    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
      cs_math_33_eigen((const cs_real_t (*)[3])cb->dpty_mat,
                       &(cb->eig_ratio),
                       &(cb->eig_max));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the curl-curl property inside a cell and its related quantities.
 *         Cellwise version using a cs_cell_mesh_t structure
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      c_flag    flag related to this cell
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_equation_set_curlcurl_property_cw(const cs_equation_param_t    *eqp,
                                     const cs_cell_mesh_t         *cm,
                                     const cs_real_t               t_eval,
                                     const cs_flag_t               c_flag,
                                     cs_cell_builder_t            *cb)
{
  cs_property_tensor_in_cell(cm,
                             eqp->curlcurl_property,
                             t_eval,
                             eqp->curlcurl_hodge.inv_pty,
                             cb->cpty_mat);

  if (cs_property_is_isotropic(eqp->curlcurl_property))
    cb->cpty_val = cb->cpty_mat[0][0];

  /* Set additional quantities in case of more advanced way of enforcing the
     essential BCs */
  if (c_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {
    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
      cs_math_33_eigen((const cs_real_t (*)[3])cb->cpty_mat,
                       &(cb->eig_ratio),
                       &(cb->eig_max));
  }
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  quant        pointer to additional mesh quantities struct.
 * \param[in]  time_step    pointer to a time step structure
 * \param[in]  eb_flag      metadata for Edge-based schemes
 * \param[in]  fb_flag      metadata for Face-based schemes
 * \param[in]  vb_flag      metadata for Vertex-based schemes
 * \param[in]  vcb_flag     metadata for Vertex+Cell-basde schemes
 * \param[in]  hho_flag     metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_init(const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *quant,
                        const cs_time_step_t         *time_step,
                        cs_flag_t                     eb_flag,
                        cs_flag_t                     fb_flag,
                        cs_flag_t                     vb_flag,
                        cs_flag_t                     vcb_flag,
                        cs_flag_t                     hho_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according to the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a new structure to handle the building of algebraic system
 *         related to an cs_equation_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_builder_t *
cs_equation_init_builder(const cs_equation_param_t   *eqp,
                         const cs_mesh_t             *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_builder_t structure
 *
 * \param[in, out]  p_builder  pointer of pointer to the cs_equation_builder_t
 *                             structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_free_builder(cs_equation_builder_t  **p_builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the residual normalization at the cellwise level according
 *         to the type of renormalization
 *         Case of scalar-valued system
 *
 * \param[in]      type            type of renormalization
 * \param[in]      vol_c           cell volume
 * \param[in]      csys            pointer to a cs_cell_sys_t structure
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the  residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_cw_scal_res_normalization(cs_param_resnorm_type_t    type,
                                      cs_real_t                  vol_c,
                                      const cs_cell_sys_t       *csys,
                                      const cs_real_t            weight[],
                                      cs_real_t                 *normalization);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the residual normalization at the cellwise level according
 *         to the type of renormalization
 *         Case of vector-valued system
 *
 * \param[in]      type            type of renormalization
 * \param[in]      vol_c           cell volume
 * \param[in]      csys            pointer to a cs_cell_sys_t structure
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the  residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_cw_vect_res_normalization(cs_param_resnorm_type_t    type,
                                      cs_real_t                  vol_c,
                                      const cs_cell_sys_t       *csys,
                                      const cs_real_t            weight[],
                                      cs_real_t                 *normalization);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage to compute of the renormalization coefficient for the
 *         the residual norm of the linear system
 *
 * \param[in]      type            type of renormalization
 * \param[in]      rhs_size        size of the rhs array
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the  residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_res_normalization(cs_param_resnorm_type_t    type,
                                   cs_lnum_t                  rhs_size,
                                   const cs_real_t            rhs[],
                                   cs_real_t                 *normalization);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare a linear system and synchronize buffers to handle parallelism.
 *        Transfer a mesh-based description of arrays x0 and rhs into an
 *        algebraic description for the linear system in x and b.
 *
 * \param[in]      stride   stride to apply to the range set operations
 * \param[in]      x_size   size of the vector unknows (scatter view)
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      rset     pointer to a range set structure
 * \param[in, out] x        array of unknows (in: initial guess)
 * \param[in, out] b        right-hand side
 *
 * \returns the number of non-zeros in the matrix
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_equation_prepare_system(int                     stride,
                           cs_lnum_t               x_size,
                           const cs_matrix_t      *matrix,
                           const cs_range_set_t   *rset,
                           cs_real_t              *x,
                           cs_real_t              *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO scheme
 *
 * \param[in]  n_scatter_dofs local number of DoFs (may be != n_gather_elts)
 * \param[in]  eqp            pointer to a cs_equation_param_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  rs             pointer to a cs_range_set_t structure
 * \param[in]  normalization  value used for the residual normalization
 * \param[in, out] x          solution of the linear system (in: initial guess)
 * \param[in, out] b          right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_solve_scalar_system(cs_lnum_t                     n_scatter_dofs,
                                const cs_equation_param_t    *eqp,
                                const cs_matrix_t            *matrix,
                                const cs_range_set_t         *rset,
                                cs_real_t                     normalization,
                                cs_real_t                    *x,
                                cs_real_t                    *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a message in the performance output file related to the
 *          monitoring of equation
 *
 * \param[in]  eqname    pointer to the name of the current equation
 * \param[in]  eqb       pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_write_monitoring(const char                    *eqname,
                             const cs_equation_builder_t   *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize all properties for an algebraic system
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure (diffusion
 *                           property is stored inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_properties(const cs_equation_param_t     *eqp,
                            const cs_equation_builder_t   *eqb,
                            const cs_real_t                t_eval,
                            cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize all properties for a given cell when building the
 *         algebraic system. If the property is uniform, a first call has to
 *         be done before the loop on cells
 *         Call \ref cs_eqution_init_properties for instance
 *
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in]      eqb        pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in]      cell_flag  flag related to the current cell
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 *                            (diffusion property is stored inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_properties_cw(const cs_equation_param_t     *eqp,
                               const cs_equation_builder_t   *eqb,
                               const cs_real_t                t_eval,
                               const cs_flag_t                cell_flag,
                               const cs_cell_mesh_t          *cm,
                               cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the diffusion property inside a cell and its related quantities
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      c_id    id of the cell to deal with
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_diffusion_property(const cs_equation_param_t   *eqp,
                                   cs_lnum_t                    c_id,
                                   cs_real_t                    t_eval,
                                   cs_flag_t                    c_flag,
                                   cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the diffusion property inside a cell and its related quantities.
 *         Cellwise version using a cs_cell_mesh_t structure
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_diffusion_property_cw(const cs_equation_param_t     *eqp,
                                      const cs_cell_mesh_t          *cm,
                                      cs_real_t                      t_eval,
                                      cs_flag_t                      c_flag,
                                      cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account the enforcement of internal DoFs. Apply an
 *          algebraic manipulation
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforced_internal_dofs(const cs_equation_param_t       *eqp,
                                   cs_cell_builder_t               *cb,
                                   cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account the enforcement of internal DoFs. Case of matrices
 *          defined by blocks. Apply an algebraic manipulation.
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforced_internal_block_dofs(const cs_equation_param_t       *eqp,
                                         cs_cell_builder_t               *cb,
                                         cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the 2*n_cells
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_tmpbuf(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_equation_get_tmpbuf_size(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_equation_balance_t structure
 *
 * \param[in]  location   where the balance is performed
 * \param[in]  size       size of arrays in the structure
 *
 * \return  a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_balance_t *
cs_equation_balance_create(cs_flag_t    location,
                           cs_lnum_t    size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a cs_equation_balance_t structure
 *
 * \param[in, out] b     pointer to a cs_equation_balance_t to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_balance_reset(cs_equation_balance_t   *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize balance terms if this is a parallel computation
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in, out] b         pointer to a cs_equation_balance_t to rsync
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_balance_sync(const cs_cdo_connect_t    *connect,
                         cs_equation_balance_t     *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_balance_t structure
 *
 * \param[in, out]  p_balance  pointer to the pointer to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_balance_destroy(cs_equation_balance_t   **p_balance);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the definition to consider at each vertex
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2v_idx   index array  to define
 * \param[in, out]  def2v_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_definitions_at_vertices(const cs_cdo_connect_t  *connect,
                                         int                      n_defs,
                                         cs_xdef_t              **defs,
                                         cs_lnum_t                def2v_idx[],
                                         cs_lnum_t                def2v_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the definition to consider at each face
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2f_idx   index array  to define
 * \param[in, out]  def2f_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_definitions_at_faces(const cs_cdo_connect_t  *connect,
                                      int                      n_defs,
                                      cs_xdef_t              **defs,
                                      cs_lnum_t                def2f_idx[],
                                      cs_lnum_t                def2f_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mean-value across ranks at each vertex
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       dim         number of entries for each vertex
 * \param[in]       counter     number of occurences on this rank
 * \param[in, out]  values      array to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_vertex_mean_values(const cs_cdo_connect_t     *connect,
                                    int                         dim,
                                    int                        *counter,
                                    cs_real_t                  *values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_COMMON_H__ */

#ifndef __CS_EQUATION_COMMON_H__
#define __CS_EQUATION_COMMON_H__

/*============================================================================
 * Routines to handle common equation features for building algebraic system
 * in CDO schemes
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_time.h"
#include "cs_domain.h"
#include "cs_equation_param.h"
#include "cs_flag.h"
#include "cs_matrix.h"
#include "cs_range_set.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_equation_builder_t
 *  \brief Store common elements used when building an algebraic system
 *  related to an equation
 */

typedef struct {

  /*!
   * @name Flags to know what to build and how to build such terms
   * @{
   */

  cs_flag_t    msh_flag;     /*!< Information related to what to build in a
                              *   \ref cs_cell_mesh_t structure for a generic
                              *   cell */
  cs_flag_t    bd_msh_flag;  /*!< Information related to what to build in a
                              *   \ref cs_cell_mesh_t structure for a cell close
                              *   to the boundary */
  cs_flag_t    st_msh_flag;  /*!< Information related to what to build in a
                              *   \ref cs_cell_mesh_t structure when only the
                              *   source term has to be built */
  cs_flag_t    sys_flag;     /*!< Information related to the sytem */

  /*!
   * @}
   * @name Metadata related to associated physical properties
   * @{
   */

  bool   diff_pty_uniform;      /*!< Is diffusion property uniform ? */
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
   * is modified to define an updated \ref cs_cdo_bc_t structure.
   */

  cs_cdo_bc_t           *face_bc; /*!< list of faces sorted by type of BCs */

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

static inline cs_flag_t
cs_equation_cell_mesh_flag(cs_flag_t                      cell_flag,
                           const cs_equation_builder_t   *eqb)
{
  cs_flag_t  _flag = eqb->msh_flag | eqb->st_msh_flag;

  if (cell_flag & CS_FLAG_BOUNDARY)
    _flag |= eqb->bd_msh_flag;

  return _flag;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 * \param[in]  quant         pointer to additional mesh quantities struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  cc            pointer to a cs_domain_cdo_context_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_allocate(const cs_cdo_connect_t          *connect,
                            const cs_cdo_quantities_t       *quant,
                            const cs_time_step_t            *time_step,
                            const cs_domain_cdo_context_t   *cc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \param[in]  cc    pointer to a structure storing CDO/HHO metadata
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_free(const cs_domain_cdo_context_t   *cc);

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
 * \brief Prepare a linear system and synchronize buffers to handle parallelism.
 *        Transfer a mesh-based description of arrays x0 and rhs into an
 *        algebraic description for the linear system in x and b.
 *
 * \param[in]      stride   stride to apply to the range set operations
 * \param[in]      x_size   size of the vector unknows (scatter view)
 * \param[in]      x0       pointer to an array (unknows to compute)
 * \param[in]      rhs      pointer to an array (right-hand side)
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      rset     pointer to a range set structure
 * \param[in, out] p_x      pointer of pointer to the linear solver unknows
 * \param[in, out] p_rhs    pointer of pointer to the right-hand side
 *
 * \returns the number of non-zeros in the matrix
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_equation_prepare_system(int                   stride,
                           cs_lnum_t             x_size,
                           const cs_real_t      *x0,
                           const cs_real_t      *rhs,
                           const cs_matrix_t    *matrix,
                           cs_range_set_t       *rset,
                           cs_real_t            *p_x[],
                           cs_real_t            *p_rhs[]);

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
 * \brief  Set members of the cs_cell_sys_t structure related to the boundary
 *         conditions. Only the generic part is done here. The remaining part
 *         is performed in _init_cell_system() for each scheme
 *
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in, out] csys      pointer to a cs_cell_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_cell_sys_bc(const cs_equation_builder_t   *eqb,
                             const cs_cell_mesh_t          *cm,
                             cs_cell_sys_t                 *csys);

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
                            cs_real_t                      t_eval,
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
 * \brief  Assemble a cellwise system into the global algebraic system
 *
 * \param[in]      csys         cellwise view of the algebraic system
 * \param[in]      rset         pointer to a cs_range_set_t structure
 * \param[in, out] mav          pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix(const cs_cell_sys_t            *csys,
                            const cs_range_set_t           *rset,
                            cs_matrix_assembler_values_t   *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system defined by blocks into the global
 *         algebraic system
 *
 * \param[in]      csys         cellwise view of the algebraic system
 * \param[in]      rset         pointer to a cs_range_set_t structure
 * \param[in]      n_x_dofs     number of DoFs per entity (= size of the block)
 * \param[in, out] mav          pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_block_matrix(const cs_cell_sys_t            *csys,
                                  const cs_range_set_t           *rset,
                                  int                             n_x_dofs,
                                  cs_matrix_assembler_values_t   *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the connectivity vertex->vertices for the local rank
 *
 * \return  a pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_adjacency_t *
cs_equation_get_v2v_index(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the connectivity face->faces for the local rank
 *
 * \return  a pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_adjacency_t *
cs_equation_get_f2f_index(void);

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

END_C_DECLS

#endif /* __CS_EQUATION_COMMON_H__ */

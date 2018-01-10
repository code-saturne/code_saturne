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

/* Store common elements used when building an algebraic system related to
   an equation */
typedef struct {

  /* Shortcut to know what to build */
  cs_flag_t    msh_flag;     // Information related to cell mesh
  cs_flag_t    bd_msh_flag;  // Information related to cell mesh (boundary)
  cs_flag_t    st_msh_flag;  // Information related to cell mesh (source term)
  cs_flag_t    sys_flag;     // Information related to the sytem

  /* Metadata related to associated properties */
  bool         diff_pty_uniform;
  bool         time_pty_uniform;
  bool         reac_pty_uniform[CS_CDO_N_MAX_REACTIONS];

  /* Source terms */
  cs_mask_t   *source_mask;  /* NULL if at least one source term is not
                                defined for all cells (size = n_cells) */

  /* Pointer to functions which compute the value of the source term */
  cs_source_term_cellwise_t  *compute_source[CS_N_MAX_SOURCE_TERMS];

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.
  */

  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs

  /* Monitoring the efficiency of the algorithm used to manipulate/build
     an equation builder. */
  cs_timer_counter_t               tcb; /* Cumulated elapsed time for building
                                           the current system */
  /* tcb >= tcd + tca + tcr + tcs */
  cs_timer_counter_t               tcd; /* Cumulated elapsed time for building
                                           diffusion terms */
  cs_timer_counter_t               tca; /* Cumulated elapsed time for building
                                           advection terms */
  cs_timer_counter_t               tcr; /* Cumulated elapsed time for building
                                           reaction terms */
  cs_timer_counter_t               tcs; /* Cumulated elapsed time for building
                                           source terms */

  cs_timer_counter_t               tce; /* Cumulated elapsed time for computing
                                           all extra operations (post, balance,
                                           fluxes...) */

} cs_equation_builder_t;

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
cs_equation_get_cell_mesh_flag(cs_flag_t                      cell_flag,
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
 * \param[in]  connect          pointer to a cs_cdo_connect_t structure
 * \param[in]  quant            pointer to additional mesh quantities struct.
 * \param[in]  time_step        pointer to a time step structure
 * \param[in]  cdo_context      pointer to a cs_domain_cdo_context_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_allocate(const cs_cdo_connect_t         *connect,
                            const cs_cdo_quantities_t      *quant,
                            const cs_time_step_t           *time_step,
                            const cs_domain_cdo_context_t  *cdo_context);

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
 * \param[in, out] tpty_val  pointer to the value for the time property
 * \param[in, out] rpty_vals pointer to the values for reaction properties
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure (diffusion
 *                           property is stored inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_properties(const cs_equation_param_t     *eqp,
                            const cs_equation_builder_t   *eqb,
                            double                        *tpty_val,
                            double                        *rpty_vals,
                            cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the diffusion property inside a cell and its related quantities
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      c_id    id of the cell to deal with
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_diffusion_property(const cs_equation_param_t     *eqp,
                                   cs_lnum_t                      c_id,
                                   cs_flag_t                      c_flag,
                                   cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the diffusion property inside a cell and its related quantities.
 *         Cellwise version using a cs_cell_mesh_t structure
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_diffusion_property_cw(const cs_equation_param_t     *eqp,
                                      const cs_cell_mesh_t          *cm,
                                      cs_flag_t                      c_flag,
                                      cs_cell_builder_t             *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system related to cell vertices into the global
 *         algebraic system
 *
 * \param[in]       csys      cellwise view of the algebraic system
 * \param[in]       rset      pointer to a cs_range_set_t structure on vertices
 * \param[in]       eqp       pointer to a cs_equation_param_t structure
 * \param[in, out]  rhs       array storing the right-hand side
 * \param[in, out]  sources   array storing the contribution of source terms
 * \param[in, out]  mav       pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_v(const cs_cell_sys_t            *csys,
                       const cs_range_set_t           *rset,
                       const cs_equation_param_t      *eqp,
                       cs_real_t                      *rhs,
                       cs_real_t                      *sources,
                       cs_matrix_assembler_values_t   *mav);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system related to cell faces into the global
 *         algebraic system
 *
 * \param[in]      csys         cellwise view of the algebraic system
 * \param[in]      rset         pointer to a cs_range_set_t structure
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      n_face_dofs  number of DoFs for each face
 * \param[in, out] rhs          array storing the right-hand side
 * \param[in, out] mav          pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_f(const cs_cell_sys_t            *csys,
                       const cs_range_set_t           *rset,
                       const cs_equation_param_t      *eqp,
                       int                             n_face_dofs,
                       cs_real_t                      *rhs,
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

END_C_DECLS

#endif /* __CS_EQUATION_COMMON_H__ */

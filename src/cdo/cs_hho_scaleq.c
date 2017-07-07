/*============================================================================
 * Build an algebraic system for scalar conv./diff. eq. with Hybrid High Order
 * space discretization
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_scheme_geometry.h"
#include "cs_equation_common.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_search.h"
#include "cs_sla.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hho_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_HHO_SCALEQ_DBG  0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for HHO discretization */

struct _cs_hho_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* System size */
  cs_lnum_t  n_faces;
  cs_lnum_t  n_cells;
  cs_lnum_t  n_dofs;
  int        n_max_fcbyc;  // n_max_fbyc + 1

  /* Shortcut to know what to build */
  cs_flag_t    msh_flag;     // Information related to cell mesh
  cs_flag_t    bd_msh_flag;  // Information related to cell mesh (boundary)
  cs_flag_t    st_msh_flag;  // Information related to cell mesh (source term)
  cs_flag_t    sys_flag;     // Information related to the sytem

  /* Store the matrix to invert after assembling and static condensation for
     upper left part
     Store in the right and left bottom part quantities attached to the system
     in order to recover the value at each cell centers */
  cs_sla_hmatrix_t      *hybrid_storage;

  /* Store the values of the field at cell centers */
  cs_real_t             *face_values;

  /* Common members for all terms */
  double                *loc_vals; // local temporary values

  /* Source terms */
  cs_real_t     *source_terms;
  cs_mask_t     *source_mask;  /* NULL if at least one source term is not
                                  defined for all cells (size = n_cells) */

  /* Pointer to functions which compute the value of the source term */
  cs_source_term_cellwise_t   *compute_source[CS_N_MAX_SOURCE_TERMS];

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to border vertices
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs, the treatment is more complex since the
     contributions of these BCs to a dual cell related to a border vertex is
     computed from the contribution of different portions of primal faces
     (these primal border faces form a closure of the dual cell).
     These contributions are computed on the fly.

   */

  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs
  double                *dir_val; // size = vtx_dir->n_nhmg_elts

  /* In case of a weak enforcement of the Dirichlet BCs */
  cs_lnum_t             *c2bcbf_idx;  // size: n_cells + 1
  cs_lnum_t             *c2bcbf_ids;  // cell --> border faces ids

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t  **cs_hho_cell_systems = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;


/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                  const cs_cdo_connect_t       *connect,
                                  const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_initialize(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers and generic structures related to HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_finalize(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_hho_scaleq_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_hho_scaleq_init(const cs_equation_param_t   *eqp,
                   const cs_mesh_t             *mesh)
{
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_HHO && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued HHO equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_b_faces = connect->n_faces[1];
  const cs_lnum_t  n_i_faces = connect->n_faces[2];
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_hho_scaleq_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_hho_scaleq_t);

  /* Shared pointers */
  b->eqp = eqp;

  /* System dimension */
  b->n_faces = connect->n_faces[0];
  b->n_cells = n_cells;
  b->n_dofs = b->n_faces; // * CS_HHO_N_FACE_DOFS[scheme_order]
  b->n_max_fcbyc = connect->n_max_fbyc + 1;

  /* Store a direct access to which term one has to compute
     High-level information on how to build the current system */
  b->sys_flag = 0;
  if (eqp->flag & CS_EQUATION_DIFFUSION)
    b->sys_flag |= CS_FLAG_SYS_DIFFUSION;
  if (eqp->flag & CS_EQUATION_CONVECTION)
    b->sys_flag |= CS_FLAG_SYS_ADVECTION;
  if (eqp->flag & CS_EQUATION_REACTION)
    b->sys_flag |= CS_FLAG_SYS_REACTION;
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    b->sys_flag |= CS_FLAG_SYS_TIME;
  if (eqp->n_source_terms > 0)
    b->sys_flag |= CS_FLAG_SYS_SOURCETERM;

  /* Flag to indicate what to build in a cell mesh */
  b->msh_flag = CS_CDO_LOCAL_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  b->bd_msh_flag = CS_CDO_LOCAL_FE;

  /* Default intialization */
  b->st_msh_flag = 0;
  b->source_terms = NULL;

  if (b->sys_flag & CS_FLAG_SYS_SOURCETERM) {

    /* Default intialization */
    b->st_msh_flag = cs_source_term_init(CS_SPACE_SCHEME_HHO,
                                         eqp->n_source_terms,
                                         (const cs_xdef_t **)eqp->source_terms,
                                         b->compute_source,
                                         &(b->sys_flag),
                                         &(b->source_mask));

    BFT_MALLOC(b->source_terms, b->n_dofs, cs_real_t);
#pragma omp parallel for if (b->n_dofs > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_dofs; i++)
      b->source_terms[i] = 0;

  } /* There is at least one source term */

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_hho_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_hho_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_hho_scaleq_free(void   *builder)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  if (b == NULL)
    return b;

  /* eqp is only shared. This structure is freed later. */

  /* Last free */
  BFT_FREE(b);
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in, out] builder     pointer to a cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_compute_source(void            *builder)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  if (eqp->n_source_terms == 0)
    return;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a HHO scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_hho_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_build_system(const cs_mesh_t         *mesh,
                           const cs_real_t         *field_val,
                           double                   dt_cur,
                           void                    *builder,
                           cs_real_t               *rhs,
                           cs_matrix_t             *matrix)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;


}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values required for hybrid discretization
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_update_field(const cs_real_t     *solu,
                           const cs_real_t     *rhs,
                           void                *builder,
                           cs_real_t           *field_val)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t  *)builder;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at faces (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in]  builder    pointer to a builder structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

double *
cs_hho_scaleq_get_face_values(const void          *builder)
{
  const cs_hho_scaleq_t  *b = (const cs_hho_scaleq_t  *)builder;

  if (b == NULL)
    return NULL;
  else
    return b->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_extra_op(const char            *eqname,
                       const cs_field_t      *field,
                       void                  *builder)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t  *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  // TODO
  CS_UNUSED(eqp);
  CS_UNUSED(field);
  CS_UNUSED(eqname);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS

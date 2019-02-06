#ifndef __CS_CDOFB_NAVSTO_H__
#define __CS_CDOFB_NAVSTO_H__

/*============================================================================
 * Routines shared among all face-based schemes for the discretization of the
 * Navier--Stokes system
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_field.h"
#include "cs_math.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_navsto_param.h"
#include "cs_time_step.h"
#include "cs_sdm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence vector associated to the current cell.
 *         WARNING: mind that, differently form the original definition, the
 *         result here is not divided by the cell volume
 *
 * \param[in]      cm         pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out] div        array related to the divergence operator
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_cdofb_navsto_divergence_vect(const cs_cell_mesh_t  *cm,
                                cs_real_t              div[])
{
  /* D(\hat{u}) = \frac{1}{|c|} \sum_{f_c} \iota_{fc} u_f.f
   * But, when integrating [[ p, q ]]_{P_c} = |c| p_c q_c
   * Thus, the volume in the divergence drops
   */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_real_t  i_f = cm->f_sgn[f] * pfq.meas;

    cs_real_t  *_div_f = div + 3*f;
    _div_f[0] = i_f * pfq.unitv[0];
    _div_f[1] = i_f * pfq.unitv[1];
    _div_f[2] = i_f * pfq.unitv[2];

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add contribution related to the pressure if Nitsche's method for the
 *         boundary conditions is requested
 *
 * \param[in]       eqp         pointer to \ref cs_equation_param_t structure
 * \param[in]       cm          pointer to \ref cs_cell_mesh_t structure
 * \param[in]       prs_c       value of the pressure at the current cell
 * \param[in, out]  csys        pointer to \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_cdofb_navsto_pressure_nitsche(const cs_equation_param_t *eqp,
                                 const cs_cell_mesh_t      *cm,
                                 const cs_real_t            prs_c,
                                 cs_cell_sys_t             *csys)
{
  /* Boundary condition contribution to the algebraic system
   * Operations that have to be performed BEFORE the static condensation */
  if ((csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)           &&
      cs_equation_param_has_diffusion(eqp)                        &&
      (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
       eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
     ) {
    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */
      const short int  f = csys->_f_ids[i];

      if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {
        const cs_quant_t pfq = cm->face[f];
        const cs_real_t f_prs = pfq.meas * prs_c;
        cs_real_t *f_rhs = csys->rhs + 3*f;
        f_rhs[0] -= f_prs * pfq.unitv[0];
        f_rhs[1] -= f_prs * pfq.unitv[1];
        f_rhs[2] -= f_prs * pfq.unitv[2];
      } /* If Dirichlet boundary face */

    } /* Loop on i */

  } /* Boundary cell */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence of a cell using the \ref cs_cdo_quantities_t
 *         structure
 *
 * \param[in]     c_id         cell id
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_dof        values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdofb_navsto_cell_divergence(const cs_lnum_t               c_id,
                                const cs_cdo_quantities_t    *quant,
                                const cs_adjacency_t         *c2f,
                                const cs_real_t              *f_dof);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add contribution related to the pressure if Nitsche's method for the
 *         boundary conditions (Dirichlet or Sliding) is requested
 *
 * \param[in]       eqp         pointer to \ref cs_equation_param_t structure
 * \param[in]       cm          pointer to \ref cs_cell_mesh_t structure
 * \param[in]       prs_c       value of the pressure at the current cell
 * \param[in, out]  csys        pointer to \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_pressure_nitsche(const cs_equation_param_t *eqp,
                                 const cs_cell_mesh_t      *cm,
                                 const cs_real_t            prs_c,
                                 cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the grad-div part to the local matrix (i.e. for the current
 *         cell)
 *
 * \param[in]      n_fc       local number of faces for the current cell
 * \param[in]      zeta       scalar coefficient for the grad-div operator
 * \param[in]      div        divergence
 * \param[in, out] mat        local system matrix to update
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_add_grad_div(short int          n_fc,
                             const cs_real_t    zeta,
                             const cs_real_t    div[],
                             cs_sdm_t          *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the pressure values
 *
 * \param[in]       nsp     pointer to a \ref cs_navsto_param_t structure
 * \param[in]       quant   pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]       ts      pointer to a \ref cs_time_step_t structure
 * \param[in, out]  pr      pointer to the pressure \ref cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_pressure(const cs_navsto_param_t     *nsp,
                              const cs_cdo_quantities_t   *quant,
                              const cs_time_step_t        *ts,
                              cs_field_t                  *pr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a zero-mean
 *         average
 *
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_set_zero_mean_pressure(const cs_cdo_quantities_t  *quant,
                                       cs_real_t                   values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform extra-operation related to Fb schemes when solving
 *         Navier-Stokes.
 *         - Compute the mass flux accross the boundaries.
 *
 * \param[in]  nsp        pointer to a \ref cs_navsto_param_t struct.
 * \param[in]  quant      pointer to a \ref cs_cdo_quantities_t struct.
 * \param[in]  connect    pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  adv_field  pointer to a \ref cs_adv_field_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_extra_op(const cs_navsto_param_t     *nsp,
                         const cs_cdo_quantities_t   *quant,
                         const cs_cdo_connect_t      *connect,
                         const cs_adv_field_t        *adv_field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         an algebraic technique.
 *         One assumes that static condensation has been performed and that
 *         the velocity-block has size 3*n_fc
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_alge(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a penalization technique (with a large coefficient).
 *         One assumes that static condensation has been performed and that
 *         the velocity-block has size 3*n_fc
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_pena(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size 3*(n_fc + 1)
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_weak(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (symmetrized Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size 3*(n_fc + 1)
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_wsym(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a symmetric boundary (treated as a sliding BCs on
 *         the three velocity components.
 *         A weak penalization technique (symmetrized Nitsche) is used.
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has (n_fc + 1) blocks of size 3x3.
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_symmetry(short int                       f,
                  const cs_equation_param_t      *eqp,
                  const cs_cell_mesh_t           *cm,
                  cs_cell_builder_t              *cb,
                  cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account a wall BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          Case of vector-valued CDO Face-based schemes
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_fixed_wall(short int                       f,
                    const cs_equation_param_t      *eqp,
                    const cs_cell_mesh_t           *cm,
                    cs_cell_builder_t              *cb,
                    cs_cell_sys_t                  *csys);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOFB_NAVSTO_H__ */

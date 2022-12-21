#ifndef __CS_CDOCB_SCALEQ_H__
#define __CS_CDOCB_SCALEQ_H__

/*============================================================================
 * Build an algebraic CDO cell-based system for the diffusion equations
 * and solved it as one block (monolithic approach of the flux-potential
 * coupling)
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
#include "cs_cdocb_priv.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_mesh.h"
#include "cs_time_step.h"
#include "cs_domain.h"
/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_cdocb_t  cs_cdocb_scaleq_t;

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the values of the velocity on the faces
 *
 * \param[in] previous        retrieve the previous state (true/false)
 *
 * \return a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/


/*============================================================================
 * Public function prototypes
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the related variables-fields
 *         to previous values
 *         Case of the monolithic coupling algorithm.
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme context structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_init_sharing(const cs_cdo_quantities_t     *quant,
                             const cs_cdo_connect_t        *connect,
                             const cs_time_step_t          *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared pointers with lifecycle dedicated to this file
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_finalize_sharing();

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdocb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdocb_monolithic_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_free_context(void *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdocb_monolithic_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_monolithic_free_scheme_context(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed.
 *         Case of scalar-valued CDO-Cb schemes
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 * \param[in, out] vtx_bc_flag     pointer to an array of BC flag for each vtx
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_setup(cs_real_t                      t_eval,
                      const cs_mesh_t               *mesh,
                      const cs_equation_param_t     *eqp,
                      cs_equation_builder_t         *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Cb schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion term in the
 *          scalar-valued CDO-Cb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_diffusion(const cs_equation_param_t     *eqp,
                          const cs_equation_builder_t   *eqb,
                          const cs_cdocb_scaleq_t       *eqc,
                          const cs_cell_mesh_t          *cm,
                          cs_hodge_t                    *diff_hodge,
                          cs_cell_sys_t                 *csys,
                          cs_cell_builder_t             *cb);


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Stokes or Oseen system with a CDO Cell-based scheme
 *         using a monolithic approach.
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_face_values(void    *context,
                                bool     previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_cell_values(void        *context,
                                bool         previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain
 *         Case of scalar-valued CDO cell-based scheme
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme builder structure
 *
 * \return a pointer to a \ref cs_cdo_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdocb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Neumann BCs on the potential, which is equivelant
 *         to applying Dirichlet BCs on the flux
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_neumann(const cs_equation_param_t      *eqp,
                 const cs_cell_mesh_t           *cm,
                 const cs_face_mesh_t           *fm,
                 const cs_hodge_t               *hodge,
                 cs_cell_builder_t              *cb,
                 cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(fm);
  CS_UNUSED(hodge);

  double  *x_neu = cb->values;
  double  *ax_neu = cb->values + 1;
  cs_sdm_t  *m = csys->mat;

  for (short int f = 0; f < csys->n_dofs; f++) {

    if (   csys->bf_flag[f] & CS_CDO_BC_NEUMANN
        || csys->bf_flag[f] & CS_CDO_BC_HMG_NEUMANN) {

      /* Build x_neu */
      bool  is_non_homogeneous = false; /* Assume homogeneous by default */

      memset(cb->values, 0, 2*sizeof(double));

      if (csys->bf_flag[f] & CS_CDO_BC_NEUMANN) {
        *x_neu = csys->neu_values[f];
        is_non_homogeneous = true;

        csys->rhs[cm->n_fc] -= -*x_neu*cm->f_sgn[f];
      }

      if (is_non_homogeneous) {

        for (int i = 0; i < m->n_rows; i++) {

          if (i == f)
            continue;

          *ax_neu = *x_neu * m->val[i*m->n_rows + f];
          csys->rhs[i] -= *ax_neu;

        }

      } /* Non-homogeneous Dirichlet BC */

      /* Set RHS to the Dirichlet value for the related face */

      csys->rhs[f] = *x_neu;

      /* Second pass: Replace the Dirichlet entry by a one and fill with
       * zero the remaining row and column */

      for (int i = 0; i < m->n_rows; i++) {

        if (i != f) {
          m->val[i*m->n_rows + f] = 0.0;
        }
        else { /* i == f */

          for (int j = 0; j < m->n_cols; j++)
            m->val[f*m->n_rows + j] = 0.0;

          m->val[f*m->n_rows + f] = 1.0;
        }

      } /* row i */
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the potential
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system

 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_dirichlet(const cs_equation_param_t      *eqp,
                   const cs_cell_mesh_t           *cm,
                   const cs_face_mesh_t           *fm,
                   const cs_hodge_t               *hodge,
                   cs_cell_builder_t              *cb,
                   cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(fm);
  CS_UNUSED(hodge);

  cs_sdm_t  *m = csys->mat;

  for (short int f = 0; f < csys->n_dofs; f++) {

    if (   csys->bf_flag[f] & CS_CDO_BC_DIRICHLET
        || csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Set RHS to the Dirichlet value for the related face */

      csys->rhs[f] = -cm->f_sgn[f]*csys->dir_values[f];
    }
  }
}

END_C_DECLS

#endif /* __CS_CDOCB_SCALEQ_H__ */

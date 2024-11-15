#ifndef __CS_CDOCB_SCALEQ_H__
#define __CS_CDOCB_SCALEQ_H__

/*============================================================================
 * Build an algebraic CDO cell-based system for the diffusion equations
 * and solved it as one block (monolithic approach of the flux-potential
 * coupling)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the generic structures for building a CDO-Cb scheme are
 *        allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdocb_scaleq_is_initialized(void);

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
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in] mesh       pointer to the mesh structure
 * \param[in] cdoq       additional CDO mesh quantities
 * \param[in] connect    pointer to a \ref cs_cdo_connect_t struct.
 * \param[in] time_step  pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_init_sharing(const cs_mesh_t            *mesh,
                             const cs_cdo_quantities_t  *cdoq,
                             const cs_cdo_connect_t     *connect,
                             const cs_time_step_t       *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared pointers with lifecycle dedicated to this file
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_finalize_sharing(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a cs_cdocb_scaleq_t structure storing data useful for
 *        building and managing such a scheme
 *
 * \param[in, out] eqp       set of parameters related an equation
 * \param[in]      var_id    id of the variable field
 * \param[in]      bflux_id  id of the boundary flux field
 * \param[in, out] eqb       pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_init_context(cs_equation_param_t    *eqp,
                             int                     var_id,
                             int                     bflux_id,
                             cs_equation_builder_t  *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_cdocb_scaleq_t structure
 *
 * \param[in, out] scheme_context  pointer to a scheme context to free
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_free_context(void  *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the boundary conditions known from the settings
 *        Define an indirection array for the enforcement of internal DoFs
 *        only if needed.
 *        Case of scalar-valued CDO-Cb schemes
 *
 * \param[in]      t_eval  time at which one evaluates BCs
 * \param[in]      mesh    pointer to a cs_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in, out] eqb     pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_setup(cs_real_t                   t_eval,
                      const cs_mesh_t            *mesh,
                      const cs_equation_param_t  *eqp,
                      cs_equation_builder_t      *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the initial values of the variable field taking into account
 *        the boundary conditions.
 *        Case of scalar-valued CDO-Cb schemes.
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
cs_cdocb_scaleq_init_values(cs_real_t                   t_eval,
                            const int                   field_id,
                            const cs_mesh_t            *mesh,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context);

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
 * \brief  Solve the steady-state equation with a CDO cell-based scheme
 *         Scalar-valued diffusion equation up-to-now
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdocb_scaleq_t structure
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
 * \brief  Retrieve an array of values at mesh cellss for the variable field
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
 * \brief  Retrieve an array of values at mesh faces for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_face_values(void        *context,
                                bool         previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux accross (primal) faces
 *         A scalar-valued flux for each face.
 *         Case of scalar-valued CDO-Cb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  diff_flux   values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_diff_flux_faces(const cs_real_t             *values,
                                const cs_equation_param_t   *eqp,
                                cs_real_t                    t_eval,
                                cs_equation_builder_t       *eqb,
                                void                        *context,
                                cs_real_t                   *diff_flux);

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

END_C_DECLS

#endif /* __CS_CDOCB_SCALEQ_H__ */

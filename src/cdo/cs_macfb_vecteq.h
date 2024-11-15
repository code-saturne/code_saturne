#ifndef __CS_MACFB_VECTEQ_H__
#define __CS_MACFB_VECTEQ_H__

/*============================================================================
 * Build an algebraic MAC face-based system for unsteady convection/diffusion
 * reaction of vector-valued equations with source terms
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
#include "cs_cdo_assembly.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_builder.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_macfb_priv.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_restart.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for MAC face-based discretization */

typedef struct _cs_macfb_t cs_macfb_vecteq_t;

/*============================================================================
 * Static inline public function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a MAC scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool cs_macfb_vecteq_is_initialized(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         vector-valued face-based schemes.
 *         Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_init_sharing(const cs_cdo_quantities_t *quant,
                                  const cs_cdo_connect_t    *connect,
                                  const cs_time_step_t      *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   double pointer to a \ref cs_cell_sys_t structure
 * \param[out]  cb     double pointer to a \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_get(cs_cell_sys_t **csys, cs_cell_builder_t **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO face-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_finalize_sharing(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_macfb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void *cs_macfb_vecteq_init_context(cs_equation_param_t   *eqp,
                                   int                    var_id,
                                   int                    bflux_id,
                                   cs_equation_builder_t *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_macfb_vecteq_t structure
 *
 * \param[in, out]  data   pointer to a cs_macfb_vecteq_t structure
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

void *cs_macfb_vecteq_free_context(void *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of vector-valued MAC schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_init_values(cs_real_t                  t_eval,
                                 const int                  field_id,
                                 const cs_mesh_t           *mesh,
                                 const cs_equation_param_t *eqp,
                                 cs_equation_builder_t     *eqb,
                                 void                      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *          The algebraic system for time t^{n+1} is going to be built knowing
 *          previous field at time t^{n} and potentially the field at time
 *          t^{n-1}. Make sure to be consistent between the call to
 *          current_to_previous and the parameters vel_{f}_n/nm1 given
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      val_f_n     face DoFs at time step n
 * \param[in]      val_f_nm1   face DoFs at time step n-1 or null
 * \param[in, out] macb        pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_init_cell_system(const cs_cell_mesh_t        *cm,
                                      const cs_equation_param_t   *eqp,
                                      const cs_equation_builder_t *eqb,
                                      const cs_real_t              val_f_n[],
                                      const cs_real_t              val_f_nm1[],
                                      cs_macfb_builder_t          *macb,
                                      cs_cell_sys_t               *csys,
                                      cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed. This is stored inside eqb
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_setup(cs_real_t                  t_eval,
                           const cs_mesh_t           *mesh,
                           const cs_equation_param_t *eqp,
                           cs_equation_builder_t     *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize stuctures for a gven cell
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      c_id         cell id
 * \param[in]      vel_f_n      velocity face DoFs of the previous time step
 * \param[in, out] cm           pointer to a cellwise view of the mesh
 * \param[in, out] macb         pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_init_build(const cs_cdo_connect_t      *connect,
                                const cs_cdo_quantities_t   *quant,
                                const cs_equation_param_t   *eqp,
                                const cs_equation_builder_t *eqb,
                                const cs_lnum_t              c_id,
                                const cs_real_t              vel_f_n[],
                                cs_cell_mesh_t              *cm,
                                cs_macfb_builder_t          *macb,
                                cs_cell_sys_t               *csys,
                                cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the source term for a vector-valued MAC scheme
 *         and add it to the local rhs
 *
 * \param[in]      cm          pointer to a \ref cs_cell_mesh_t structure
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      t_eval      time at which the source term is evaluated
 * \param[in]      coef        scaling of the time source (for theta schemes)
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] cb          pointer to a \ref cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_sourceterm(const cs_cell_mesh_t      *cm,
                                const cs_equation_param_t *eqp,
                                cs_macfb_builder_t        *macb,
                                const cs_real_t            t_eval,
                                const cs_real_t            coef,
                                cs_equation_builder_t     *eqb,
                                cs_cell_builder_t         *cb,
                                cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the matrix and rhs for a vector-valued MAC scheme
 *          and Euler implicit. Values are added in place
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      cb          pointer to a \ref cs_cell_builder_t structure
 * \param[in]      dt          value of the time step
 * \param[in, out] csys        pointer to a \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_euler_implicit_term(const cs_equation_param_t *eqp,
                                         const cs_cell_mesh_t      *cm,
                                         const cs_macfb_builder_t  *macb,
                                         const cs_cell_builder_t   *cb,
                                         const cs_real_t            dt,
                                         cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion term in the
 *          vector-valued MAC schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      diff_pty    pointer to a cs_property_t structure
 *                              for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_diffusion(const cs_equation_param_t *eqp,
                               const cs_cell_mesh_t      *cm,
                               const cs_macfb_builder_t  *macb,
                               const cs_property_t       *diff_pty,
                               cs_cell_sys_t             *csys,
                               cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the convection term in the
 *          vector-valued MAC-Fb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_advection(const cs_equation_param_t *eqp,
                               const cs_macfb_vecteq_t   *eqc,
                               const cs_cell_mesh_t      *cm,
                               const cs_macfb_builder_t  *macb,
                               cs_cell_sys_t             *csys,
                               cs_cell_builder_t         *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the convection, diffusion,
 *          reaction terms in vector-valued MAC-Fb schemes.
 *          mass_hodge could be set to null if a Voronoi algo. is used.
 *          Otherwise, the mass matrix should be pre-computed.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      macb        pointer to a cs_macfb_builder_t structure
 * \param[in]      diff_pty    pointer to a cs_property_data_t structure
 *                              for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_conv_diff_reac(const cs_equation_param_t   *eqp,
                                    const cs_equation_builder_t *eqb,
                                    const cs_macfb_vecteq_t     *eqc,
                                    const cs_cell_mesh_t        *cm,
                                    const cs_macfb_builder_t    *macb,
                                    const cs_property_data_t    *diff_pty,
                                    cs_cell_sys_t               *csys,
                                    cs_cell_builder_t           *cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with MAC scheme
 *
 * \param[in]      csys         pointer to a cs_cell_sys_t structure
 * \param[in, out] block        pointer to a block structure
 * \param[in, out] rhs          array of values for the rhs
 * \param[in, out] asb          pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_assembly(const cs_cell_sys_t   *csys,
                              cs_cdo_system_block_t *block,
                              cs_real_t             *rhs,
                              cs_cdo_assembly_t     *asb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the variables associated to cells in case of a MAC-fb
 *         scheme. This has to be done after a resolution.
 *
 * \param[in, out] tce       pointer to a timer counter
 * \param[in, out] fld       pointer to a cs_field_t structure
 * \param[in]      cur2prev  true if one performs "current to previous" op.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_update_fields(cs_timer_counter_t *tce,
                                   cs_field_t         *fld,
                                   bool                cur2prev);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector diffusion
 *         equation with a MAC-Fb scheme:
 *           - steady scheme
 *           - implicit Euler scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_solve_steady_implicit(bool                       cur2prev,
                                           const cs_mesh_t           *mesh,
                                           const int                  field_id,
                                           const cs_equation_param_t *eqp,
                                           cs_equation_builder_t     *eqb,
                                           void                      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector diffusion
 *         equation with a MAC scheme and an implicit/explicit theta scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_solve_theta(bool                       cur2prev,
                                 const cs_mesh_t           *mesh,
                                 const int                  field_id,
                                 const cs_equation_param_t *eqp,
                                 cs_equation_builder_t     *eqb,
                                 void                      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_current_to_previous(const cs_equation_param_t *eqp,
                                         cs_equation_builder_t     *eqb,
                                         void                      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_macfb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_extra_post(const cs_equation_param_t *eqp,
                                cs_equation_builder_t     *eqb,
                                void                      *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at mesh cells from the inverse operation
 *         w.r.t. the static condensation (DoF used in the linear system are
 *         located at primal faces)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size 3*n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *cs_macfb_vecteq_get_cell_values(void *context, bool previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size 3*n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *cs_macfb_vecteq_get_face_values(void *context, bool previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_read_restart(cs_restart_t *restart,
                                  const char   *eqname,
                                  void         *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_vecteq_write_restart(cs_restart_t *restart,
                                   const char   *eqname,
                                   void         *scheme_context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MACFB_VECTEQ_H__ */

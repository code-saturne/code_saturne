#ifndef __CS_MACFB_NAVSTO_H__
#define __CS_MACFB_NAVSTO_H__

/*============================================================================
 * Functions shared for MAC face-based schemes for the discretization of the
 * Navier--Stokes system
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
#include "cs_cdo_quantities.h"
#include "cs_cdofb_navsto.h"
#include "cs_field.h"
#include "cs_iter_algo.h"
#include "cs_math.h"
#include "cs_matrix.h"
#include "cs_mesh.h"
#include "cs_navsto_param.h"
#include "cs_sdm.h"
#include "cs_time_plot.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \enum cs_macfb_navsto_boussinesq_type_t
 * \brief Type of algorithm to compute the Boussinesq approximation
 */

typedef cs_cdofb_navsto_boussinesq_type_t cs_macfb_navsto_boussinesq_type_t;

/*! \struct cs_macfb_navsto_builder_t
 *
 * \brief Structure storing additional arrays related to the building of the
 *        Navier-Stokes system.
 *
 * This structure is associated to a cell-wise building in case of MAC
 * face-based scheme.
 */

typedef cs_cdofb_navsto_builder_t cs_macfb_navsto_builder_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute and add a source term to the local RHS.
 *         This is a special treatment to enable source involving face DoFs and
 *         potentially the local discrete divergence/gradient operators.
 *         Examples are gravity term or Bousinesq term(s)
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void(cs_macfb_navsto_source_t)(const cs_navsto_param_t         *nsp,
                                       const cs_cell_mesh_t            *cm,
                                       const cs_macfb_navsto_builder_t *nsb,
                                       cs_cell_sys_t                   *csys);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the way to compute the Boussinesq approximation
 *
 * \param[in] type     type of algorithm to use
 */
/*----------------------------------------------------------------------------*/

void
cs_macfb_navsto_set_boussinesq_algo(cs_macfb_navsto_boussinesq_type_t type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and allocate a local NavSto builder when Fb schemes are used
 *
 * \param[in] nsp         set of parameters to define the NavSto system
 * \param[in] connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a cs_macfb_navsto_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_macfb_navsto_builder_t
cs_macfb_navsto_create_builder(const cs_navsto_param_t *nsp,
                               const cs_cdo_connect_t  *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy the given cs_macfb_navsto_builder_t structure
 *
 * \param[in, out] nsb   pointer to the cs_macfb_navsto_builder_t to free
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_free_builder(cs_macfb_navsto_builder_t *nsb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the members of the cs_macfb_navsto_builder_t structure
 *
 * \param[in]      t_eval     time at which one evaluates the pressure BC
 * \param[in]      nsp        set of parameters to define the NavSto system
 * \param[in]      cm         cellwise view of the mesh
 * \param[in]      macb       macfb builder
 * \param[in]      csys       cellwise view of the algebraic system
 * \param[in]      pr_bc      set of definitions for the presuure BCs
 * \param[in]      bf_type    type of boundaries for all boundary faces
 * \param[in, out] nsb        builder to update
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_define_builder(cs_real_t                  t_eval,
                                    const cs_navsto_param_t   *nsp,
                                    const cs_cell_mesh_t      *cm,
                                    const cs_macfb_builder_t  *macb,
                                    const cs_cell_sys_t       *csys,
                                    const cs_cdo_bc_face_t    *pr_bc,
                                    const cs_boundary_type_t  *bf_type,
                                    cs_macfb_navsto_builder_t *nsb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mass flux playing the role of the advection field in
 *         the Navier-Stokes equations
 *         One considers the mass flux across primal faces which relies on the
 *         velocity vector defined on each face.
 *
 * \param[in]      nsp         set of parameters to define the NavSto system
 * \param[in]      quant       set of additional geometrical quantities
 * \param[in]      face_vel    velocity vectors for each face
 * \param[in, out] mass_flux   array of mass flux values to update (allocated)
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_mass_flux(const cs_navsto_param_t   *nsp,
                               const cs_cdo_quantities_t *quant,
                               const cs_real_t           *face_vel,
                               cs_real_t                 *mass_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence in a cell of a vector-valued array defined at
 *         faces (values are defined both at interior and border faces).
 *         Variant based on the usage of \ref cs_cdo_quantities_t structure.
 *
 * \param[in]     c_id         cell id
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_vals       values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

double cs_macfb_navsto_cell_divergence(const cs_lnum_t            c_id,
                                       const cs_cdo_quantities_t *quant,
                                       const cs_adjacency_t      *c2f,
                                       const cs_real_t           *f_vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the pressure at faces
 *
 * \param[in]       mesh      pointer to a cs_mesh_t structure
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]       ts        pointer to a cs_time_step_t structure
 * \param[in]       nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]       p_cell    value of the pressure inside each cell
 * \param[in, out]  p_face    value of the pressure at each face
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_compute_face_pressure(const cs_mesh_t           *mesh,
                                           const cs_cdo_connect_t    *connect,
                                           const cs_cdo_quantities_t *quant,
                                           const cs_time_step_t      *ts,
                                           const cs_navsto_param_t   *nsp,
                                           const cs_real_t           *p_cell,
                                           cs_real_t                 *p_face);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the velocity at cells
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]       mass_flux scalar-valued mass flux for each face
 * \param[in, out]  u_cell    vector-value velocity velo at each cell
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_compute_cell_velocity(const cs_cdo_connect_t    *connect,
                                           const cs_cdo_quantities_t *quant,
                                           const cs_real_t           *mass_flux,
                                           cs_real_t                 *u_cell);

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

void cs_macfb_navsto_add_grad_div(short int       n_fc,
                                  const cs_real_t zeta,
                                  const cs_real_t div[],
                                  cs_sdm_t       *mat);

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

void cs_macfb_navsto_init_pressure(const cs_navsto_param_t   *nsp,
                                   const cs_cdo_quantities_t *quant,
                                   const cs_time_step_t      *ts,
                                   cs_field_t                *pr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the pressure values when the pressure is defined at
 *         faces
 *
 * \param[in]       nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]       connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]       ts        pointer to a \ref cs_time_step_t structure
 * \param[in, out]  pr_f      pointer to the pressure values at faces
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_init_face_pressure(const cs_navsto_param_t *nsp,
                                        const cs_cdo_connect_t  *connect,
                                        const cs_time_step_t    *ts,
                                        cs_real_t               *pr_f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a mean-value
 *         equal to the reference value
 *
 * \param[in]       nsp       pointer to a cs_navsto_param_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_rescale_pressure_to_ref(const cs_navsto_param_t   *nsp,
                                             const cs_cdo_quantities_t *quant,
                                             cs_real_t values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a zero-mean
 *         average
 *
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_set_zero_mean_pressure(const cs_cdo_quantities_t *quant,
                                            cs_real_t values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform extra-operation related to MAC-Fb schemes when solving
 *         Navier-Stokes. Computation of the following quantities according to
 *         post-processing flags beeing activated.
 *         - The mass flux accross the boundaries.
 *         - The global mass in the computational domain
 *         - The norm of the velocity divergence
 *         - the cellwise mass flux balance
 *         - the kinetic energy
 *         - the velocity gradient
 *         - the pressure gradient
 *         - the vorticity
 *         - the helicity
 *         - the enstrophy
 *         - the stream function
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t struct.
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      quant         pointer to a \ref cs_cdo_quantities_t struct.
 * \param[in]      connect       pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]      ts            pointer to a \ref cs_time_step_t struct.
 * \param[in,out]  time_plotter  pointer to a \ref cs_time_plot_t struct.
 * \param[in]      adv_field     pointer to a \ref cs_adv_field_t struct.
 * \param[in]      mass_flux     scalar-valued mass flux for each face
 * \param[in]      p_cell        scalar-valued pressure in each cell
 * \param[in]      u_face        vector-valued velocity on each face
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_extra_op(const cs_navsto_param_t   *nsp,
                              const cs_mesh_t           *mesh,
                              const cs_cdo_quantities_t *quant,
                              const cs_cdo_connect_t    *connect,
                              const cs_time_step_t      *ts,
                              cs_time_plot_t            *time_plotter,
                              const cs_adv_field_t      *adv_field,
                              const cs_real_t           *mass_flux,
                              const cs_real_t           *p_cell,
                              const cs_real_t           *u_face);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         an algebraic technique.
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_block_dirichlet_alge(short int                  f,
                                   const cs_equation_param_t *eqp,
                                   const cs_cell_mesh_t      *cm,
                                   const cs_property_data_t  *pty,
                                   cs_cell_builder_t         *cb,
                                   cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a penalization technique (with a large coefficient).
 *         One assumes that static condensation has been performed and that
 *         the velocity-block has size 3*n_fc
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_block_dirichlet_pena(short int                  f,
                                   const cs_equation_param_t *eqp,
                                   const cs_cell_mesh_t      *cm,
                                   const cs_property_data_t  *pty,
                                   cs_cell_builder_t         *cb,
                                   cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size 3*(n_fc + 1)
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_block_dirichlet_weak(short int                  fb,
                                   const cs_equation_param_t *eqp,
                                   const cs_cell_mesh_t      *cm,
                                   const cs_property_data_t  *pty,
                                   cs_cell_builder_t         *cb,
                                   cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (symmetrized Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size 3*(n_fc + 1)
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_block_dirichlet_wsym(short int                  fb,
                                   const cs_equation_param_t *eqp,
                                   const cs_cell_mesh_t      *cm,
                                   const cs_property_data_t  *pty,
                                   cs_cell_builder_t         *cb,
                                   cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a boundary defined as 'symmetry' (treated as a
 *         sliding BCs on the three velocity components.)
 *         A weak penalization technique (symmetrized Nitsche) is used.
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_symmetry(short int                  fb,
                       const cs_equation_param_t *eqp,
                       const cs_cell_mesh_t      *cm,
                       const cs_property_data_t  *pty,
                       cs_cell_builder_t         *cb,
                       cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account a wall BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_fixed_wall(short int                  fb,
                         const cs_equation_param_t *eqp,
                         const cs_cell_mesh_t      *cm,
                         const cs_property_data_t  *pty,
                         cs_cell_builder_t         *cb,
                         cs_cell_sys_t             *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one has to do one more non-linear iteration.
 *         Test if performed on the relative norm on the increment between
 *         two iterations
 *
 * \param[in]      nl_algo_type   type of non-linear algorithm
 * \param[in]      pre_iterate    previous state of the mass flux iterate
 * \param[in]      cur_iterate    current state of the mass flux iterate
 * \param[in, out] algo           pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_macfb_navsto_nl_algo_cvg(cs_param_nl_algo_t nl_algo_type,
                            const cs_real_t   *pre_iterate,
                            cs_real_t         *cur_iterate,
                            cs_iter_algo_t    *algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the function pointer computing the source term in the momentum
 *         equation related to the gravity effect (hydrostatic pressure or the
 *         Boussinesq approximation)
 *
 * \param[in]  nsp          set of parameters for the Navier-Stokes system
 * \param[out] p_func       way to compute the gravity effect
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_set_gravity_func(const cs_navsto_param_t   *nsp,
                                      cs_macfb_navsto_source_t **p_func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the gravity effects.
 *         Compute and add the source term to the local RHS.
 *         This is a special treatment since of face DoFs are involved
 *         contrary to the standard case where only the cell DoFs is involved.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_gravity_term(const cs_navsto_param_t         *nsp,
                                  const cs_cell_mesh_t            *cm,
                                  const cs_macfb_navsto_builder_t *nsb,
                                  cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This is the standard case where the face DoFs are used for the
 *         constant part rho0 . g[] and only the cell DoFs are involved for the
 *         remaining part (the Boussinesq approximation).
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_boussinesq_at_cell(const cs_navsto_param_t         *nsp,
                                        const cs_cell_mesh_t            *cm,
                                        const cs_macfb_navsto_builder_t *nsb,
                                        cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This way to compute the Boussinesq approximation relies only on DoFs
 *         at faces. This should enable to keep a stable (no velocity) in case
 *         of a stratified configuration.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_boussinesq_at_face(const cs_navsto_param_t         *nsp,
                                        const cs_cell_mesh_t            *cm,
                                        const cs_macfb_navsto_builder_t *nsb,
                                        cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This way to compute the Boussinesq approximation relies only on DoFs
 *         at faces. This should enable to keep a stable (no velocity) in case
 *         of a stratified configuration.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_boussinesq_by_part(const cs_navsto_param_t         *nsp,
                                        const cs_cell_mesh_t            *cm,
                                        const cs_macfb_navsto_builder_t *nsb,
                                        cs_cell_sys_t                   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the source term for computing the stream function.
 *         This relies on the prototype associated to the generic function
 *         pointer \ref cs_dof_func_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         null or pointer to a structure cast on-the-fly
 * \param[in, out] retval        result of the function. Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_navsto_stream_source_term(cs_lnum_t        n_elts,
                                        const cs_lnum_t *elt_ids,
                                        bool             dense_output,
                                        void            *input,
                                        cs_real_t       *retval);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MACFB_NAVSTO_H__ */

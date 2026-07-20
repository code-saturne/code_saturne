#ifndef CS_CDOFB_NAVSTO_H
#define CS_CDOFB_NAVSTO_H

/*============================================================================
 * Functions shared among all face-based schemes for the discretization of the
 * Navier--Stokes system
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_matrix.h"
#include "base/cs_base.h"
#include "base/cs_field.h"
#include "base/cs_math.h"
#include "base/cs_time_plot.h"
#include "base/cs_time_step.h"
#include "cdo/cs_cdo_connect.h"
#include "cdo/cs_cdo_quantities.h"
#include "cdo/cs_cdo_turbulence.h"
#include "cdo/cs_iter_algo.h"
#include "cdo/cs_navsto_context.h"
#include "cdo/cs_navsto_param.h"
#include "cdo/cs_sdm.h"
#include "mesh/cs_mesh.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * \enum cs_cdofb_navsto_boussinesq_type_t
 * \brief Type of algorithm to compute the Boussinesq approximation
 */

enum cs_cdofb_navsto_boussinesq_type_t {

  /*!
   * \brief Boussinesq approximation relyong on a cell contribution
   *
   * This algorithm uses cell DoFs for the Boussinesq part corresponding to
   * rho0 * beta * (var[c_id] - var0) * g[] while the constant part equal to
   * rho0 * g[] is built in order to be in balance with the pressure gradient
   * (face DoF contributions).
   */

  CS_CDOFB_NAVSTO_BOUSSINESQ_CELL_DOF,

  /*!
   * \brief Boussinesq approximation relyong on face contributions
   *
   * This algorithm uses only face DoFs for the Boussinesq approximation.
   * For the constant part (rho0 * g[]) as well as the variable part
   * rho0 * beta * (var[c_id] - var0) * g[]
   * The aim is to be in balance with the pressure gradient
   */

  CS_CDOFB_NAVSTO_BOUSSINESQ_FACE_DOF

};


/*! \struct cs_cdofb_navsto_builder_t
 *
 * \brief Structure storing additional arrays related to the building of the
 *        Navier-Stokes system.
 *
 * This structure is associated to a cell-wise building in case of CDO
 * face-based scheme.
 */

struct cs_cdofb_navsto_builder_t {

  /*!
   * @name Properties
   * @{
   * \var rho_c
   * Value of the mass density for the current cell
   *
   */

  cs_real_t            rho_c;

  /*!
   * @}
   * @name Operator(s)
   * @{
   *
   * \var div_op
   * Cell-wise divergence operator (in fact div_op = -|c| div).
   * This operator is stored in an array of size 3*n_fc
   *
   */

  cs_real_t           *div_op;

  /*!
   * @}
   * @name Source term for the mass equation
   * @{
   *
   * \var mass_rhs
   * Value of the rhs dedicated to the mass equation
   */

  cs_real_t            mass_rhs;

  /*!
   * @}
   * @name Boundary conditions
   * @{
   *
   * \var bf_type
   * Type of boundary to apply to each face.
   *
   * Array of size n_fc.  Zero is set for interior faces.
   *
   *
   * \var pressure_bc_val
   * Store the value of the pressure on boundary faces.
   *
   * Array of size n_fc. Only useful if a Dirichlet boundary condition is set
   * on a boundary face.
   */

  cs_boundary_type_t  *bf_type;          /* Size: n_fc */
  cs_real_t           *pressure_bc_val;  /* Size: n_fc */

  /*!
   * @}
   */

};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and add a source term to the local RHS.
 *        This is a special treatment to enable the computation of source terms
 *        involving face DoFs and potentially the local discrete divergence or
 *        gradient operators. In the standard case, only the cell DoFs are
 *        involved.
 *        Examples of this function pointer are given for the treatment of the
 *        gravity term or the Bousinesq term(s)
 *
 * \param[in]      nsp   set of parameters to handle the Navier-Stokes system
 * \param[in]      cm    pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb   pointer to a builder structure for the NavSto system
 * \param[in, out] csys  pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdofb_navsto_source_t)
(
 const cs_navsto_param_t                          *nsp,
 const cs_cell_mesh_t                             *cm,
 [[maybe_unused]] const cs_cdofb_navsto_builder_t *nsb,
 cs_cell_sys_t                                    *csys
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

void
cs_cdofb_navsto_set_boussinesq_algo(cs_cdofb_navsto_boussinesq_type_t type);

cs_cdofb_navsto_builder_t
cs_cdofb_navsto_create_builder(const cs_navsto_param_t *nsp,
                               const cs_cdo_connect_t  *connect);

void
cs_cdofb_navsto_free_builder(cs_cdofb_navsto_builder_t *nsb);

void
cs_cdofb_navsto_define_builder(cs_real_t                  t_eval,
                               const cs_navsto_param_t   *nsp,
                               const cs_cell_mesh_t      *cm,
                               const cs_cell_sys_t       *csys,
                               const cs_cdo_bc_face_t    *pr_bc,
                               const cs_boundary_type_t  *bf_type,
                               cs_cdofb_navsto_builder_t *nsb);

void
cs_cdofb_navsto_mass_flux(const cs_navsto_param_t   *nsp,
                          const cs_cdo_quantities_t *cdoq,
                          const cs_real_t           *face_vel,
                          cs_real_t                 *mass_flux);

double
cs_cdofb_navsto_cell_divergence(const cs_lnum_t            c_id,
                                const cs_cdo_quantities_t *cdoq,
                                const cs_adjacency_t      *c2f,
                                const cs_real_t           *f_vals);

void
cs_cdofb_navsto_compute_face_pressure(const cs_mesh_t           *mesh,
                                      const cs_cdo_connect_t    *connect,
                                      const cs_cdo_quantities_t *cdoq,
                                      const cs_time_step_t      *ts,
                                      const cs_navsto_param_t   *nsp,
                                      const cs_real_t           *p_cell,
                                      cs_real_t                 *p_face);

void
cs_cdofb_navsto_add_grad_div(short int        n_fc,
                             const cs_real_t  zeta,
                             const cs_real_t  div[],
                             cs_sdm_t        *mat);

void
cs_cdofb_navsto_init_pressure(const cs_navsto_param_t   *nsp,
                              const cs_cdo_quantities_t *cdoq,
                              const cs_time_step_t      *ts,
                              cs_field_t                *pr);

void
cs_cdofb_navsto_check_init(const cs_navsto_param_t   *nsp,
                           const cs_cdo_connect_t    *connect,
                           const cs_cdo_quantities_t *cdoq,
                           const cs_time_step_t      *ts,
                           const cs_field_t          *velocity,
                           const cs_real_t           *face_vel,
                           const cs_field_t          *pressure);

void
cs_cdofb_navsto_init_face_pressure(const cs_navsto_param_t *nsp,
                                   const cs_cdo_connect_t  *connect,
                                   const cs_time_step_t    *ts,
                                   cs_real_t               *pr_f);

void
cs_cdofb_navsto_rescale_pressure_to_ref(const cs_navsto_param_t   *nsp,
                                        const cs_cdo_quantities_t *cdoq,
                                        cs_real_t                  values[]);

void
cs_cdofb_navsto_set_zero_mean_pressure(const cs_cdo_quantities_t *cdoq,
                                       cs_real_t                  values[]);

void
cs_cdofb_navsto_extra_op(const cs_navsto_param_t           *nsp,
                         const cs_mesh_t                   *mesh,
                         const cs_cdo_quantities_t         *cdoq,
                         const cs_cdo_connect_t            *connect,
                         const cs_time_step_t              *ts,
                         cs_time_plot_t                    *time_plotter,
                         const cs_turbulence_t             *turb,
                         const cs_adv_field_t              *adv_field,
                         const cs_real_t                   *mass_flux,
                         const cs_real_t                   *p_cell,
                         const cs_real_t                   *u_cell,
                         const cs_real_t                   *u_face,
                         const cs_cdo_navsto_psteady_cvg_t &ps_cvg);

void
cs_cdofb_block_dirichlet_alge(short int                                   bf,
                              [[maybe_unused]] const cs_equation_param_t *eqp,
                              [[maybe_unused]] const cs_cell_mesh_t      *cm,
                              [[maybe_unused]] const cs_property_data_t  *pty,
                              [[maybe_unused]] cs_cell_builder_t         *cb,
                              cs_cell_sys_t                              *csys);

void
cs_cdofb_block_dirichlet_pena(short int                                   bf,
                              [[maybe_unused]] const cs_equation_param_t *eqp,
                              [[maybe_unused]] const cs_cell_mesh_t      *cm,
                              [[maybe_unused]] const cs_property_data_t  *pty,
                              [[maybe_unused]] cs_cell_builder_t         *cb,
                              cs_cell_sys_t                              *csys);

void
cs_cdofb_block_dirichlet_weak(short int                                   bf,
                              [[maybe_unused]] const cs_equation_param_t *eqp,
                              [[maybe_unused]] const cs_cell_mesh_t      *cm,
                              [[maybe_unused]] const cs_property_data_t  *pty,
                              [[maybe_unused]] cs_cell_builder_t         *cb,
                              cs_cell_sys_t                              *csys);

void
cs_cdofb_block_dirichlet_wsym(short int                                   bf,
                              [[maybe_unused]] const cs_equation_param_t *eqp,
                              [[maybe_unused]] const cs_cell_mesh_t      *cm,
                              [[maybe_unused]] const cs_property_data_t  *pty,
                              [[maybe_unused]] cs_cell_builder_t         *cb,
                              cs_cell_sys_t                              *csys);
void
cs_cdofb_symmetry_weak(short int                                   bf,
                       [[maybe_unused]] const cs_equation_param_t *eqp,
                       [[maybe_unused]] const cs_cell_mesh_t      *cm,
                       [[maybe_unused]] const cs_property_data_t  *pty,
                       [[maybe_unused]] cs_cell_builder_t         *cb,
                       cs_cell_sys_t                              *csys);

void
cs_cdofb_symmetry_alge(short int                                   bf,
                       [[maybe_unused]] const cs_equation_param_t *eqp,
                       [[maybe_unused]] const cs_cell_mesh_t      *cm,
                       [[maybe_unused]] const cs_property_data_t  *pty,
                       [[maybe_unused]] cs_cell_builder_t         *cb,
                       cs_cell_sys_t                              *csys);

void
cs_cdofb_prescribed_smooth_wall_n_pena_t_robin
(
 short int                                   bf,
 [[maybe_unused]] const cs_equation_param_t *eqp,
 [[maybe_unused]] const cs_cell_mesh_t      *cm,
 [[maybe_unused]] const cs_property_data_t  *pty,
 [[maybe_unused]] cs_cell_builder_t         *cb,
 cs_cell_sys_t                              *csys
);

void
cs_cdofb_prescribed_smooth_wall_n_alge_t_robin
(
 short int                                   bf,
 [[maybe_unused]] const cs_equation_param_t *eqp,
 [[maybe_unused]] const cs_cell_mesh_t      *cm,
 [[maybe_unused]] const cs_property_data_t  *pty,
 [[maybe_unused]] cs_cell_builder_t         *cb,
 cs_cell_sys_t                              *csys
);

void
cs_cdofb_prescribed_smooth_wall_n_alge_t_neumann
(
 short int                                   fb,
 [[maybe_unused]] const cs_equation_param_t *eqp,
 [[maybe_unused]] const cs_cell_mesh_t      *cm,
 [[maybe_unused]] const cs_property_data_t  *pty,
 cs_cell_builder_t                          *cb,
 cs_cell_sys_t                              *csys);

void
cs_cdofb_prescribed_smooth_wall_n_weak_t_neumann
(
 short int                                   bf,
 [[maybe_unused]] const cs_equation_param_t *eqp,
 [[maybe_unused]] const cs_cell_mesh_t      *cm,
 [[maybe_unused]] const cs_property_data_t  *pty,
 [[maybe_unused]] cs_cell_builder_t         *cb,
 cs_cell_sys_t             *csys
);

cs_sles_convergence_state_t
cs_cdofb_navsto_nl_algo_cvg(const cs_navsto_param_t *nsp,
                            const cs_real_t         *pre_iterate,
                            cs_real_t               *cur_iterate,
                            cs_iter_algo_t          *algo);

void
cs_cdofb_navsto_set_gravity_func(const cs_navsto_param_t   *nsp,
                                 cs_cdofb_navsto_source_t **p_func);

void
cs_cdofb_navsto_gravity_term
(
 const cs_navsto_param_t                          *nsp,
 const cs_cell_mesh_t                             *cm,
 [[maybe_unused]] const cs_cdofb_navsto_builder_t *nsb,
 cs_cell_sys_t                                    *csys
);

void
cs_cdofb_navsto_boussinesq_at_cell
(
 const cs_navsto_param_t                          *nsp,
 const cs_cell_mesh_t                             *cm,
 [[maybe_unused]] const cs_cdofb_navsto_builder_t *nsb,
 cs_cell_sys_t                                    *csys
);

void
cs_cdofb_navsto_boussinesq_at_face
(
 const cs_navsto_param_t                          *nsp,
 const cs_cell_mesh_t                             *cm,
 [[maybe_unused]] const cs_cdofb_navsto_builder_t *nsb,
 cs_cell_sys_t                                    *csys
);

void
cs_cdofb_navsto_stream_source_term(cs_lnum_t        n_elts,
                                   const cs_lnum_t *elt_ids,
                                   bool             dense_output,
                                   void            *input,
                                   cs_real_t       *retval);

cs_cdo_balance_t *
cs_cdofb_navsto_balance(const cs_navsto_param_t   *nsp,
                        const cs_cdo_connect_t    *connect,
                        const cs_cdo_quantities_t *cdoq,
                        const cs_time_step_t      *ts,
                        const cs_cdo_bc_face_t    *pr_bc,
                        const cs_boundary_type_t  *bf_type,
                        const cs_real_t           *pr_c);

void
cs_cdofb_navsto_check_convergence(const cs_navsto_param_t     *nsp,
                                  const cs_cdo_quantities_t   *cdoq,
                                  const cs_time_step_t        *ts,
                                  const cs::cdo_navsto_ctx_t  *sc,
                                  const cs_turbulence_t       *tbs,
                                  cs_cdo_navsto_psteady_cvg_t &ps_cvg);

/*----------------------------------------------------------------------------*/

#endif /* CS_CDOFB_NAVSTO_H */

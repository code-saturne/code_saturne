#ifndef CS_EQUATION_H
#define CS_EQUATION_H

/*============================================================================
 * Functions to handle the cs_equation_t structure and its related structures
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_field.h"
#include "base/cs_restart.h"
#include "base/cs_time_step.h"
#include "cdo/cs_cdo_connect.h"
#include "cdo/cs_cdo_quantities.h"
#include "cdo/cs_cdo_toolbox.h"
#include "cdo/cs_equation_builder.h"
#include "cdo/cs_equation_param.h"
#include "mesh/cs_mesh.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_equation_t cs_equation_t;

/*!
 * \struct cs_equation_core_t
 * \brief Main structures on which an equation structure relies
 *
 * Intermediate structure useful to manipulate an array of (sub-)structures.
 * Especially, the scheme context relies on the space discretization and it is
 * not easy to manipulate void ** object. This is a work-around to this
 * operation.
 *
 * These three structures allow one to use nearly all operations related to
 * an equation without having to build an equation structure. This is useful
 * when handling extra-diagonal block in systems of equations.
 *
 * \var cs_equation_core_t::param
 *      Set of parameters to specifiy the settings
 *
 * \var cs_equation_core_t::builder
 *      Part of the quantities useful to build/manipulate an equation. All
 *      quantities that are shared among all discretizations are in this
 *      structure.
 *
 * \var cs_equation_core_t::scheme_context
 *      Part of the quantities useful to build/manipulate an equation. All
 *      quantities that are specific to the discrization are in this structure.
 */

struct cs_equation_core_t {

  cs_equation_param_t   *param;
  cs_equation_builder_t *builder;
  void                  *scheme_context;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

int
cs_equation_get_n_equations(void);

cs_equation_t *
cs_equation_by_name(const char *eqname);

cs_equation_t *
cs_equation_by_field_name(const char *field_name);

bool
cs_equation_has_field_name(const cs_equation_t *eq,
                           const char          *fld_name);

cs_equation_param_t *
cs_equation_param_by_name(const char *eqname);

cs_equation_param_t *
cs_equation_param_by_field_name(const char *field_name);

cs_equation_param_t *
cs_equation_get_param(const cs_equation_t *eq);

cs_equation_t *
cs_equation_by_id(int eq_id);

const char *
cs_equation_get_name(const cs_equation_t *eq);

int
cs_equation_get_id(const cs_equation_t *eq);

cs_field_t *
cs_equation_get_field(const cs_equation_t *eq);

int
cs_equation_get_field_id(const cs_equation_t *eq);

const char *
cs_equation_get_field_name(const cs_equation_t *eq);

const cs_range_set_t *
cs_equation_get_range_set(const cs_equation_t *eq);

cs_gnum_t
cs_equation_get_global_n_dofs(const cs_equation_t       *eq,
                              const cs_cdo_quantities_t *cdoq);

cs_field_t *
cs_equation_get_boundary_flux(const cs_equation_t *eq);

cs_flag_t
cs_equation_get_flag(const cs_equation_t *eq);

void
cs_equation_set_flag(cs_equation_t *eq,
                     cs_flag_t      flag);

void
cs_equation_add_build_hook(cs_equation_t            *eq,
                           void                     *context,
                           cs_equation_build_hook_t *func);

cs_equation_builder_t *
cs_equation_get_builder(const cs_equation_t *eq);

void *
cs_equation_get_scheme_context(const cs_equation_t *eq);

cs_equation_core_t
cs_equation_get_core_structure(const cs_equation_t *eq);

cs_real_t *
cs_equation_get_source_term_array(const cs_equation_t *eq);

bool
cs_equation_is_steady(const cs_equation_t *eq);

bool
cs_equation_uses_new_mechanism(const cs_equation_t *eq);

cs_equation_t *
cs_equation_add(const char         *eqname,
                const char         *varname,
                cs_equation_type_t  eqtype,
                int                 dim,
                cs_param_bc_type_t  default_bc);

cs_equation_t *
cs_equation_add_user(const char         *eqname,
                     const char         *varname,
                     int                 dim,
                     cs_param_bc_type_t  default_bc);

cs_equation_t *
cs_equation_add_user_tracer(const char         *eqname,
                            const char         *varname,
                            int                 dim,
                            cs_param_bc_type_t  default_bc,
                            cs_property_t      *time_pty,
                            cs_adv_field_t     *adv,
                            cs_property_t      *diff_pty);
void
cs_equation_destroy_all(void);

bool
cs_equation_needs_steady_state_solve(void);

void
cs_equation_get_count(int *n_equations,
                      int *n_predef_equations,
                      int *n_user_equations);

void
cs_equation_log_monitoring(const cs_time_step_t      *ts,
                           const cs_cdo_quantities_t *cdoq);

void
cs_equation_log_setup(void);

void
cs_equation_set_default_param(cs_equation_key_t  key,
                              const char        *keyval);

void
cs_equation_set_sles(void);

void
cs_equation_init_sharing(const cs_mesh_t           *mesh,
                         const cs_cdo_connect_t    *connect,
                         const cs_cdo_quantities_t *cdoq,
                         const cs_time_step_t      *time_step,
                         cs_flag_t                  cb_scheme_flag,
                         cs_flag_t                  eb_scheme_flag,
                         cs_flag_t                  fb_scheme_flag,
                         cs_flag_t                  vb_scheme_flag,
                         cs_flag_t                  vcb_scheme_flag,
                         cs_flag_t                  hho_scheme_flag,
                         cs_flag_t                  mac_scheme_flag);

void
cs_equation_finalize_sharing(cs_flag_t cb_scheme_flag,
                             cs_flag_t eb_scheme_flag,
                             cs_flag_t fb_scheme_flag,
                             cs_flag_t vb_scheme_flag,
                             cs_flag_t vcb_scheme_flag,
                             cs_flag_t hho_scheme_flag,
                             cs_flag_t mac_scheme_flag);

bool
cs_equation_set_functions(void);

void
cs_equation_lock_settings(void);

void
cs_equation_predefined_create_field(int            n_previous,
                                    cs_equation_t *eq);

void
cs_equation_user_create_fields(void);

void
cs_equation_define_builders(const cs_mesh_t *mesh);

void
cs_equation_define_context_structures(void);

void
cs_equation_define_core_structure(const cs_equation_t *eq,
                                  cs_equation_core_t **p_core);

void
cs_equation_init_field_values(const cs_mesh_t      *mesh,
                              const cs_time_step_t *ts);

void
cs_equation_solve_steady_state(const cs_mesh_t *mesh,
                               cs_equation_t   *eq);

void
cs_equation_solve(bool             cur2prev,
                  const cs_mesh_t *mesh,
                  cs_equation_t   *eq);

[[deprecated]] void
cs_equation_build_system(const cs_mesh_t *mesh,
                         cs_equation_t   *eq);

[[deprecated]] void
cs_equation_solve_deprecated(const cs_mesh_t *mesh,
                             cs_equation_t   *eq);

cs_property_t *
cs_equation_get_diffusion_property(const cs_equation_t *eq);

cs_property_t *
cs_equation_get_time_property(const cs_equation_t *eq);

cs_property_t *
cs_equation_get_reaction_property(const cs_equation_t *eq,
                                  const int            reaction_id);

cs_param_time_scheme_t
cs_equation_get_time_scheme(const cs_equation_t *eq);

cs_real_t
cs_equation_get_theta_time_val(const cs_equation_t *eq);

cs_param_space_scheme_t
cs_equation_get_space_scheme(const cs_equation_t *eq);

int
cs_equation_get_space_poly_degree(const cs_equation_t *eq);

int
cs_equation_get_var_dim(const cs_equation_t *eq);

cs_equation_type_t
cs_equation_get_type(const cs_equation_t *eq);

double
cs_equation_get_time_eval(const cs_time_step_t *ts,
                          const cs_equation_t  *eq);

void
cs_equation_current_to_previous(const cs_equation_t *eq);

void
cs_equation_get_cellwise_builders(const cs_equation_t *eq,
                                  cs_cell_sys_t      **csys,
                                  cs_cell_builder_t  **cb);

cs_real_t *
cs_equation_get_cell_values(const cs_equation_t *eq,
                            bool                 previous);

cs_real_t *
cs_equation_get_face_values(const cs_equation_t *eq,
                            bool                 previous);

cs_real_t *
cs_equation_get_edge_values(const cs_equation_t *eq,
                            bool                 previous);

cs_real_t *
cs_equation_get_vertex_values(const cs_equation_t *eq,
                              bool                 previous);

void
cs_equation_integrate_variable(const cs_cdo_connect_t    *connect,
                               const cs_cdo_quantities_t *cdoq,
                               const cs_equation_t       *eq,
                               cs_real_t                 *result);

void
cs_equation_compute_boundary_diff_flux(const cs_equation_t       *eq,
                                       const cs_equation_param_t *eqp,
                                       const cs_property_t       *diff_pty,
                                       const cs_real_t           *dof_vals,
                                       const cs_real_t           *cell_vals,
                                       cs_real_t                  t_eval,
                                       cs_real_t                 *diff_flux);

void
cs_equation_compute_flux_across_plane(const cs_equation_t *eq,
                                      const char          *ml_name,
                                      const cs_real_3_t    direction,
                                      cs_real_t           *diff_flux,
                                      cs_real_t           *conv_flux);

void
cs_equation_compute_diffusive_flux(const cs_equation_t       *eq,
                                   const cs_equation_param_t *eqp,
                                   const cs_property_t       *diff_pty,
                                   const cs_real_t           *dof_vals,
                                   const cs_real_t           *cell_vals,
                                   cs_flag_t                  location,
                                   cs_real_t                  t_eval,
                                   cs_real_t                 *diff_flux);

void
cs_equation_compute_vtx_field_gradient(const cs_equation_t *eq,
                                       cs_real_t           *v_gradient);

void
cs_equation_compute_peclet(const cs_equation_t  *eq,
                           const cs_time_step_t *ts,
                           cs_real_t             peclet[]);

void
cs_equation_read_extra_restart(cs_restart_t *restart);

void
cs_equation_write_extra_restart(cs_restart_t *restart);

void
cs_equation_post_balance(const cs_mesh_t           *mesh,
                         const cs_cdo_connect_t    *connect,
                         const cs_cdo_quantities_t *cdoq,
                         const cs_time_step_t      *ts);

void
cs_equation_apply_stiffness(cs_equation_t       *eq,
                            const cs_property_t *property,
                            const cs_real_t     *pot,
                            cs_flag_t            loc_res,
                            cs_real_t           *res);

void
cs_equation_extra_post(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_EQUATION_H */

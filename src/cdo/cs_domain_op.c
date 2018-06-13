/*============================================================================
 * Manage specific post-processing related to a computational domain and
 * restart files
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_advection_field.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_gwf.h"
#include "cs_log_iteration.h"
#include "cs_post.h"
#include "cs_property.h"
#include "cs_prototypes.h"
#include "cs_navsto_system.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain_op.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for advection fields.
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure
 * \param[in]  quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 * \param[in]  dt_cur      value of the current time step
 */
/*----------------------------------------------------------------------------*/

static void
_post_advection_field(const cs_adv_field_t       *adv,
                      const cs_cdo_quantities_t  *quant,
                      const cs_time_step_t       *time_step,
                      double                      dt_cur)
{
  if (adv == NULL)
    return;
  if (adv->flag == 0)
    return;

  const bool post_courant =
    (adv->flag & CS_ADVECTION_FIELD_POST_COURANT) ? true : false;

  if (post_courant) { /* Compute and postprocess the Courant number */

    double  hc;
    cs_nvec3_t  adv_c;

    char  *label = NULL;
    double  *courant = NULL;
    int  len = strlen(adv->name) + 8 + 1;

    BFT_MALLOC(courant, quant->n_cells, double);
    BFT_MALLOC(label, len, char);
    sprintf(label, "%s.Courant", adv->name);

    if (adv->cell_field_id > -1) { /* field is defined at cell centers */

      cs_field_t  *fld = cs_field_by_id(adv->cell_field_id);

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_nvec3(fld->val + 3*c_id, &adv_c);
        hc = cbrt(quant->cell_vol[c_id]);
        courant[c_id] = dt_cur * adv_c.meas / hc;
      }

    }
    else {

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_advection_field_get_cell_vector(c_id, adv, &adv_c);
        hc = cbrt(quant->cell_vol[c_id]);
        courant[c_id] = dt_cur * adv_c.meas / hc;
      }

    }

    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      label,
                      1,
                      true,           // interlace
                      true,           // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      courant,        // values on cells
                      NULL,           // values at internal faces
                      NULL,           // values at border faces
                      time_step);     // time step management struct.

    BFT_FREE(label);
    BFT_FREE(courant);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the computational domain
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_domain_post(void                      *input,
             int                        mesh_id,
             int                        cat_id,
             int                        ent_flag[5],
             cs_lnum_t                  n_cells,
             cs_lnum_t                  n_i_faces,
             cs_lnum_t                  n_b_faces,
             const cs_lnum_t            cell_ids[],
             const cs_lnum_t            i_face_ids[],
             const cs_lnum_t            b_face_ids[],
             const cs_time_step_t      *time_step)
{
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (input == NULL)
    return;

  if (mesh_id != -1) /* Post-processing only on the generic volume mesh */
    return;

  cs_domain_t  *d = (cs_domain_t *)input;

  if (cs_domain_get_cdo_mode(d) == CS_DOMAIN_CDO_MODE_OFF)
    return;

  assert(time_step == d->time_step);

  /* Post-processing related to advection fields */
  int n_adv_fields = cs_advection_field_get_n_fields();
  for (int adv_id = 0; adv_id < n_adv_fields; adv_id++)
    _post_advection_field(cs_advection_field_by_id(adv_id),
                          d->cdo_quantities,
                          time_step,
                          d->dt_cur);

  /* Post-processing related to equations */
  cs_equation_extra_post_all(d->mesh,
                             d->connect,
                             d->cdo_quantities,
                             time_step,
                             d->dt_cur);

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the generic post-processing related to a domain
 *
 * \param[in]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post_init(cs_domain_t   *domain)
{
  /* Set pointers of function if additional postprocessing is requested */
  cs_post_add_time_mesh_dep_output(_domain_post, domain);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Process the computational domain after the resolution
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post(cs_domain_t  *domain)
{
  cs_timer_t  t0 = cs_timer_time();

  /* Pre-stage for post-processing for the current time step */
  cs_post_time_step_begin(domain->time_step);

  /* Extra-operations */
  /* ================ */

  /* Predefined extra-operations related to advection fields */
  assert(domain->cdo_context != NULL);

  if (domain->cdo_context->force_advfield_update)
    cs_advection_field_update(domain->time_step->t_cur, true);

  /* User-defined extra operations */
  cs_user_cdo_extra_op(domain);

  /* Log output */
  if (cs_domain_needs_log(domain))
    cs_log_iteration();

  /* Post-processing */
  /* =============== */

  /* Predefined extra-operations related to
     - the domain (advection fields and properties),
     - equations
     - groundwater flows
     are also handled during the call of this function thanks to
     cs_post_add_time_mesh_dep_output() function pointer
  */
  cs_post_time_step_output(domain->time_step);

  cs_post_time_step_end();

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(domain->tcp), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a restart file for the CDO/HHO module
 *
 * \param[in, out]  domain     pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_read_restart(cs_domain_t  *domain)
{
  if (cs_restart_present() == false)
    return;

  cs_restart_t  *restart = cs_restart_create("main", /* restart file name */
                                             NULL,   /* directory name */
                                             CS_RESTART_MODE_READ);

  const char err_i_val[] = N_("Restart mismatch for: %s\n"
                              "read: %d\n"
                              "expected: %d.");
  int i_val;
  int retval;

  /* Read a new section: version */
  int  version = 400000;
  retval = cs_restart_read_section(restart,
                          "code_saturne:checkpoint:main:version", // secname
                          CS_MESH_LOCATION_NONE,                  // ml_id
                          1,                                      // nb. values
                          CS_TYPE_cs_int_t,                       // val. type
                          &i_val);                                // value(s)

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);

  if (i_val != version)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "code_saturne:checkpoint:main:version", version, i_val);

  /* Read a new section: field information */
  cs_map_name_to_id_t  *old_field_map = NULL;

  cs_restart_read_field_info(restart, &old_field_map);

  /* Read a new section */
  int  n_equations = cs_equation_get_n_equations();
  retval = cs_restart_read_section(restart,
                                   "cdo:n_equations",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &i_val);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);
  if (i_val != n_equations)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "cdo:n_equations", n_equations, i_val);

  /* Read a new section */
  int  n_properties = cs_property_get_n_properties();
  retval = cs_restart_read_section(restart,
                                   "cdo:n_properties",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &i_val);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);
  if (i_val != n_properties)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "cdo:n_properties", n_properties, i_val);

  /* Read a new section */
  int  n_adv_fields = cs_advection_field_get_n_fields();
  retval = cs_restart_read_section(restart,
                                   "cdo:n_adv_fields",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &i_val);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);
  if (i_val != n_adv_fields)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "cdo:n_adv_fields", n_adv_fields, i_val);

  /* Read a new section: activation or not of the groundwater flow module */
  int  igwf = 0; /* not activated by default */
  if (cs_gwf_is_activated()) igwf = 1;
  retval = cs_restart_read_section(restart,
                                   "groundwater_flow_module",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &i_val);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);
  if (i_val != igwf)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "groundwater_flow_module", igwf, i_val);

  /* Read a new section: activation or not of the Navier-Stokes system */
  int  inss = 0; /* not activated by default */
  if (cs_navsto_system_is_activated()) inss = 1;
  retval = cs_restart_read_section(restart,
                                   "navier_stokes_system",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &i_val);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);
  if (i_val != inss)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "navier_stokes_system", inss, i_val);

  /* Read a new section: computation or not of the wall distance */
  int  iwall = 0;
  if (cs_walldistance_is_activated()) iwall = 1;
  retval = cs_restart_read_section(restart,
                                   "wall_distance",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &i_val);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);
  if (i_val != iwall)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "wall_distance", iwall, i_val);

  /* Read a new section: number of computed time steps */
  int  nt_cur = 0;
  retval = cs_restart_read_section(restart,
                                   "cur_time_step",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_int_t,
                                   &nt_cur);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);

  /* Read a new section: number of computed time steps */
  cs_real_t  t_cur = 0;
  retval = cs_restart_read_section(restart,
                                   "cur_time",
                                   CS_MESH_LOCATION_NONE,
                                   1,
                                   CS_TYPE_cs_real_t,
                                   &t_cur);

  if (retval != CS_RESTART_SUCCESS)
    bft_error(__FILE__, __LINE__, 0, " %s: error %d while reading restart",
              __func__, retval);

  domain->time_step->nt_cur = nt_cur;
  domain->time_step->t_cur = t_cur;
  cs_time_step_redefine_cur(nt_cur, t_cur);
  cs_time_step_define_prev(nt_cur, t_cur);

  /* Main variables */
  int  t_id_flag = 0; /* Only current values */
  cs_restart_read_variables(restart, old_field_map, t_id_flag, NULL);

  cs_map_name_to_id_destroy(&old_field_map);

  /* Read fields related to the main restart file */
  cs_restart_read_fields(restart, CS_RESTART_MAIN);

  /* TODO: read field values for previous time step if needed */

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    cs_field_current_to_previous(f);
  }

  /* Read additional arrays of values according to the type of
     equation and the discretization scheme */
  cs_equation_read_extra_restart(restart);

  /* Finalize restart process */
  cs_restart_destroy(&restart);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a restart file for the CDO/HHO module
 *
 * \param[in]  domain     pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_write_restart(const cs_domain_t  *domain)
{
  if (cs_restart_checkpoint_required(domain->time_step) == false)
    return;

  cs_restart_t  *restart = cs_restart_create("main", /* restart file name */
                                             NULL,   /* directory name */
                                             CS_RESTART_MODE_WRITE);

  /* Write a new section: version */
  int  version = 400000;
  cs_restart_write_section(restart,
                           "code_saturne:checkpoint:main:version", // secname
                           CS_MESH_LOCATION_NONE,                  // ml_id
                           1,                                      // nb. values
                           CS_TYPE_cs_int_t,                       // val. type
                           &version);                              // value(s)

  /* Write a new section: field information */
  cs_restart_write_field_info(restart);

  /* Write a new section */
  int  n_equations = cs_equation_get_n_equations();
  cs_restart_write_section(restart,
                           "cdo:n_equations",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &n_equations);

  /* Write a new section */
  int  n_properties = cs_property_get_n_properties();
  cs_restart_write_section(restart,
                           "cdo:n_properties",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &n_properties);

  /* Write a new section */
  int  n_adv_fields = cs_advection_field_get_n_fields();
  cs_restart_write_section(restart,
                           "cdo:n_adv_fields",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &n_adv_fields);

  /* Write a new section: activation or not of the groundwater flow module */
  int  igwf = 0; /* not activated by default */
  if (cs_gwf_is_activated()) igwf = 1;
  cs_restart_write_section(restart,
                           "groundwater_flow_module",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &igwf);

  /* Write a new section: activation or not of the Navier-Stokes system */
  int  inss = 0; /* not activated by default */
  if (cs_navsto_system_is_activated()) inss = 1;
  cs_restart_write_section(restart,
                           "navier_stokes_system",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &inss);

  /* Write a new section: computation or not of the wall distance */
  int  iwall = 0;
  if (cs_walldistance_is_activated()) iwall = 1;
  cs_restart_write_section(restart,
                           "wall_distance",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &iwall);

  /* Write a new section: number of computed time steps */
  int  ntcabs = domain->time_step->nt_cur;
  cs_restart_write_section(restart,
                           "cur_time_step",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &ntcabs);

  /* Read a new section: number of computed time steps */
  cs_real_t  ttcabs = domain->time_step->t_cur;
  cs_restart_write_section(restart,
                           "cur_time",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_real_t,
                           &ttcabs);

  /* Main variables */
  int  t_id_flag = 0; /* Only current values */
  cs_restart_write_variables(restart, t_id_flag, NULL);

  /* Write fields related to the main restart file */
  cs_restart_write_fields(restart, CS_RESTART_MAIN);

  /* TODO: write field values for previous time step if needed */

  /* Write additional arrays of values according to the type of
     equation and the discretization scheme */
  cs_equation_write_extra_restart(restart);

  /* Indicate that a chechpoint has been done */
  cs_restart_checkpoint_done(domain->time_step);

  /* Finalize restart process */
  cs_restart_destroy(&restart);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

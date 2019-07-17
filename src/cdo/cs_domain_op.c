/*============================================================================
 * Manage specific post-processing related to a computational domain and
 * restart files
 *============================================================================*/

/* VERS */

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

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_advection_field.h"
#include "cs_array_reduce.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_equation_param.h"
#include "cs_gwf.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_maxwell.h"
#include "cs_navsto_system.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_property.h"
#include "cs_prototypes.h"
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
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_cdo_post_domain(void);

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to compute at least one adimensional number
 */
/*----------------------------------------------------------------------------*/

static bool
_needs_adimensional_numbers(void)
{
  bool is_needed = false;

  /* Check for the Courant number */
  int n_adv_fields = cs_advection_field_get_n_fields();
  for (int adv_id = 0; adv_id < n_adv_fields; adv_id++) {
    const cs_adv_field_t  *adv = cs_advection_field_by_id(adv_id);
    if (adv->flag & CS_ADVECTION_FIELD_POST_COURANT)
      return true;
  }

  /* Check for the Peclet number */
  int n_equations = cs_equation_get_n_equations();
  for (int i = 0; i < n_equations; i++) {
    cs_equation_t  *eq = cs_equation_by_id(i);
    cs_equation_param_t  *eqp = cs_equation_get_param(eq);
    if (eqp->process_flag & CS_EQUATION_POST_PECLET)
      return true;
  }

  /* Check for the Fourier number */
  int  n_properties = cs_property_get_n_properties();
  for (int i = 0; i < n_properties; i++) {
    const cs_property_t  *pty = cs_property_by_id(i);
    if (pty->process_flag & CS_PROPERTY_POST_FOURIER)
      return true;
  }

  return is_needed;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute very simple statistic for an array with values located
 *         at cells (scalar-valued)
 *
 * \param[in]  cdoq       pointer to a cs_cdo_quantities_t struct.
 * \param[in]  basename   label for output in the log
 * \param[in]  array      pointer to the array to analyze
 */
/*----------------------------------------------------------------------------*/

static void
_analyze_cell_array(const cs_cdo_quantities_t   *cdoq,
                    const char                   basename[],
                    const cs_real_t              array[])
{
  cs_real_t  _min = array[0], _max = array[0], _sum = 0.;

  cs_array_reduce_simple_stats_l(cdoq->n_cells, 1, NULL, array,
                                 &_min, &_max, &_sum);

  /* Parallel treatment */
  cs_real_t  min, max, sum;
  if (cs_glob_n_ranks > 1) {

    cs_real_t  minmax[2] = {-_min, _max};

    cs_parall_max(2, CS_REAL_TYPE, minmax);
    min = -minmax[0], max = minmax[1];

    sum = _sum;
    cs_parall_sum(1, CS_REAL_TYPE, &sum);

  }
  else
    min = _min, max = _max, sum = _sum;


  cs_log_printf(CS_LOG_DEFAULT, "s- %20s  % -6.4e % -6.4e % -6.4e\n",
                basename, min, max, sum/cdoq->n_g_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined Courant number post-processing for advection fields.
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure
 * \param[in]  cdoq        pointer to a cs_cdo_quantities_t struct.
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_post_courant_number(const cs_adv_field_t       *adv,
                     const cs_cdo_quantities_t  *cdoq,
                     const cs_time_step_t       *time_step)
{
  if (adv == NULL)
    return;
  if (!(adv->flag & CS_ADVECTION_FIELD_POST_COURANT))
    return;

  int  len = strlen(adv->name) + 8 + 1;
  char  *label = NULL;
  BFT_MALLOC(label, len, char);
  sprintf(label, "%s.Courant", adv->name);

  cs_real_t  *courant = cs_equation_get_tmpbuf();

  cs_advection_get_courant(adv, time_step->dt[0], courant);

  /* Brief output for the log */
  _analyze_cell_array(cdoq, label, courant);

  /* Postprocessing */
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    label,
                    1,
                    true,       /* interlaced? */
                    true,       /* true = original mesh */
                    CS_POST_TYPE_cs_real_t,
                    courant,    /* values on cells */
                    NULL, NULL, /* values at internal,border faces */
                    time_step); /* time step management struct. */

  BFT_FREE(label);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined Peclet number post-processing for an equation
 *
 * \param[in]  eq          pointer to a cs_equation_t structure
 * \param[in]  cdoq        pointer to a cs_cdo_quantities_t struct.
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_post_peclet_number(const cs_equation_t        *eq,
                    const cs_cdo_quantities_t  *cdoq,
                    const cs_time_step_t       *time_step)
{
  if (eq == NULL)
    return;

  const cs_equation_param_t  *eqp = cs_equation_get_param(eq);

  assert(eqp != NULL);
  if (!(eqp->process_flag & CS_EQUATION_POST_PECLET))
    return;

  int  len = strlen(eqp->name) + 8 + 1;
  char  *label = NULL;
  BFT_MALLOC(label, len, char);
  sprintf(label, "%s.Peclet", eqp->name);

  cs_real_t  *peclet = cs_equation_get_tmpbuf();
  cs_equation_compute_peclet(eq, time_step, peclet);

  /* Brief output for the log */
  _analyze_cell_array(cdoq, label, peclet);

  /* Postprocessing */
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    label,
                    1,
                    true,       /* interlaced? */
                    true,       /* true = original mesh */
                    CS_POST_TYPE_cs_real_t,
                    peclet,    /* values on cells */
                    NULL, NULL, /* values at internal,border faces */
                    time_step); /* time step management struct. */

  BFT_FREE(label);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined Fourier number post-processing for a property
 *
 * \param[in]  pty         pointer to a cs_property_t structure
 * \param[in]  cdoq        pointer to a cs_cdo_quantities_t struct.
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_post_fourier_number(const cs_property_t        *pty,
                     const cs_cdo_quantities_t  *cdoq,
                     const cs_time_step_t       *time_step)
{
  if (pty == NULL)
    return;
  if (!(pty->process_flag & CS_PROPERTY_POST_FOURIER))
    return;

  cs_real_t  *fourier = cs_equation_get_tmpbuf();

  cs_property_get_fourier(pty, time_step->t_cur, time_step->dt[0], fourier);

  int  len = strlen(pty->name) + 8 + 1;
  char  *label = NULL;
  BFT_MALLOC(label, len, char);
  sprintf(label, "%s.Fourier", pty->name);

  /* Brief output for the log */
  _analyze_cell_array(cdoq, label, fourier);

  /* Postprocessing */
  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    label,
                    1,
                    true,       /* interlaced? */
                    true,       /* true = original mesh */
                    CS_POST_TYPE_cs_real_t,
                    fourier,    /* values on cells */
                    NULL, NULL, /* values at internal,border faces */
                    time_step); /* time step management struct. */

  BFT_FREE(label);
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
  CS_UNUSED(time_step);

  if (input == NULL)
    return;

  cs_domain_t  *d = (cs_domain_t *)input;

  if (cs_domain_get_cdo_mode(d) == CS_DOMAIN_CDO_MODE_OFF)
    return;

  if (mesh_id != -1) /* Post-processing only on the generic volume mesh */
    return;

  /* Additional extra-operation(s) specific to a numerical scheme */
  cs_equation_extra_post();
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve only steady-state equations
 */
/*----------------------------------------------------------------------------*/

void
cs_f_cdo_post_domain(void)
{
  cs_domain_post(cs_glob_domain);
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
 * \brief  Post-processing of the computational domain after a resolution
 *
 * \param[in]  domain            pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_post(cs_domain_t  *domain)
{
  cs_timer_t  t0 = cs_timer_time();

  assert(domain->cdo_context != NULL);

  /* Extra-operations */
  /* ================ */

  /* Predefined extra-operations related to advection fields */
  cs_advection_field_update(domain->time_step->t_cur, true);

  /* Log output */
  if (cs_domain_needs_log(domain)) {

    /* Basic statistic related to variables */
    if (domain->cdo_context->mode == CS_DOMAIN_CDO_MODE_ONLY)
      cs_log_iteration(); /* Otherwise called from the FORTRAN part */

    /* Post-processing */
    /* =============== */

    /* Post-processing of adimensional numbers */
    if (_needs_adimensional_numbers()) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         " ------------------------------------------------------------\n");
      cs_log_printf(CS_LOG_DEFAULT, "s- %20s %10s %10s %10s\n",
                    "Adim. number", "min", "max", "mean");

      /* 1. Courant numbers */
      int n_adv_fields = cs_advection_field_get_n_fields();
      for (int adv_id = 0; adv_id < n_adv_fields; adv_id++)
        _post_courant_number(cs_advection_field_by_id(adv_id),
                             domain->cdo_quantities,
                             domain->time_step);

      /* 2. Peclet numbers */
      int n_equations = cs_equation_get_n_equations();
      for (int i = 0; i < n_equations; i++)
        _post_peclet_number(cs_equation_by_id(i),
                            domain->cdo_quantities,
                            domain->time_step);

      /* 3. Fourier numbers */
      int  n_properties = cs_property_get_n_properties();
      for (int i = 0; i < n_properties; i++)
        _post_fourier_number(cs_property_by_id(i),
                             domain->cdo_quantities,
                             domain->time_step);

      cs_log_printf
        (CS_LOG_DEFAULT,
         " ------------------------------------------------------------\n");

    } /* Needs to compute adimensional numbers */

    /* 4. Equation balance */
    cs_equation_post_balance(domain->mesh,
                             domain->connect,
                             domain->cdo_quantities,
                             domain->time_step);

    /* 5.a Specific operations for the GWF module */
    if (cs_gwf_is_activated())
      cs_gwf_extra_op(domain->connect, domain->cdo_quantities);

    /* 5.b Specific operations for the Maxwell module */
    if (cs_maxwell_is_activated())
      cs_maxwell_extra_op(domain->connect, domain->cdo_quantities);

    /* 5.c Specific operations for the Navier-Stokes module */
    if (cs_navsto_system_is_activated())
      cs_navsto_system_extra_op(domain->connect, domain->cdo_quantities);

  } /* Needs a new log */

  /* Predefined extra-operations related to
     - the domain (advection fields and properties),
     - equations
     - groundwater flows
     - Maxwell module
     are also handled during the call of this function thanks to
     cs_post_add_time_mesh_dep_output() function pointer
  */
  cs_post_time_step_output(domain->time_step);

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
  if (cs_restart_present() == false) {

    /* Initialize time step if a restart frequency is set */
    cs_restart_checkpoint_set_last_ts(domain->time_step->t_cur);

    return;
  }

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

  /* Initialize time step if a restart frequency is set */
  cs_restart_checkpoint_set_last_ts(nt_cur);

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

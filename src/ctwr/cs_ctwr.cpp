/*============================================================================
 * Cooling towers related functions
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr.h"
#include "cs_ctwr_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

static cs_ctwr_option_t  _ctwr_option = {
  .evap_model = CS_CTWR_NONE,
  .has_rain = false,
  .solve_rain_velocity = false,
  .air_rain_friction = false,
  .rain_to_packing = false};

const cs_ctwr_option_t *cs_glob_ctwr_option = &_ctwr_option;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* array of exchanges area */

static int               _n_ct_zones_max = 0;
static int               _n_ct_zones     = 0;
static cs_ctwr_zone_t  **_ct_zone   = NULL;

/* Restart file */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Additional output for cooling towers
 *
 * parameters:
 *   input       <-> pointer to optional (untyped) value or structure;
 *                   here, we should point to _default_input.
 *   mesh_id     <-- id of the output mesh for the current call
 *   cat_id      <-- category id of the output mesh for the current call
 *   ent_flag    <-- indicate global presence of cells (ent_flag[0]), interior
 *                   faces (ent_flag[1]), boundary faces (ent_flag[2]),
 *                   particles (ent_flag[3]) or probes (ent_flag[4])
 *   n_cells     <-- local number of cells of post_mesh
 *   n_i_faces   <-- local number of interior faces of post_mesh
 *   n_b_faces   <-- local number of boundary faces of post_mesh
 *   cell_ids    <-- list of cells (0 to n-1) of post-processing mesh
 *   i_face_ids  <-- list of interior faces (0 to n-1) of post-processing mesh
 *   b_face_ids  <-- list of boundary faces (0 to n-1) of post-processing mesh
 *   ts          <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_write_liquid_vars(void                  *input,
                   int                    mesh_id,
                   int                    cat_id,
                   int                    ent_flag[5],
                   cs_lnum_t              n_cells,
                   cs_lnum_t              n_i_faces,
                   cs_lnum_t              n_b_faces,
                   const cs_lnum_t        cell_ids[],
                   const cs_lnum_t        i_face_ids[],
                   const cs_lnum_t        b_face_ids[],
                   const cs_time_step_t  *ts)
{
  CS_UNUSED(input);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (cat_id == CS_POST_MESH_VOLUME) {

    const cs_mesh_t *mesh = cs_glob_mesh;

    /* Liquid fraction enthalpy */

    cs_real_t *yh_l_p = (cs_real_t *)CS_F_(yh_l_pack)->val;   /* Liquid enthalpy */
    cs_real_t *y_l_p = (cs_real_t *)CS_F_(y_l_pack)->val;  /* Liquid mass per unit
                                                            cell volume */

    cs_real_t *val;
    BFT_MALLOC(val, mesh->n_cells, cs_real_t);

    /* Value on all cells */

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      val[i] = 0;

    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
      for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
        cs_lnum_t cell_id = ze_cell_ids[i];
        if (y_l_p[cell_id] > 0.)
          val[cell_id] = yh_l_p[cell_id]/y_l_p[cell_id];
      }
    }

    /* Values may be restricted to selection */

    if (cell_ids != NULL) {
      cs_real_t *_val;
      BFT_MALLOC(_val, n_cells, cs_real_t);
      for (cs_lnum_t i = 0; i < n_cells; i++)
        _val[i] = val[cell_ids[i]];
      BFT_FREE(val);
      val = _val;
    }

    const char name[] = "Enthalpy liq packing";

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      _(name),
                      1,      /* dim */
                      true,   /* interlace */
                      false,  /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      val,    /* cell values */
                      NULL,   /* internal face values */
                      NULL,   /* boundary face values */
                      ts);

    BFT_FREE(val);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the list of cells attached to a packing zone
 *         Function pointer to mesh location elements selection definition.
 *
 * \param[in]   input        pointer to a structure cast on-the-fly
 * \param[in]   m            pointer to associated mesh structure.
 * \param[in]   location_id  id of associated location.
 * \param[out]  n_elts       number of selected elements
 * \param[out]  elt_list     list of selected elements.
 */
/*----------------------------------------------------------------------------*/

static void
_packing_selection(void              *input,
                   const cs_mesh_t   *m,
                   int                location_id,
                   cs_lnum_t         *n_elts,
                   cs_lnum_t        **elt_ids)
{
  CS_UNUSED(location_id);

  const cs_ctwr_zone_t **cts = (const cs_ctwr_zone_t **)input;

  bool  *is_packing = NULL;
  BFT_MALLOC(is_packing, m->n_cells, bool);

#   pragma omp parallel for if (m->n_cells> CS_THR_MIN)
  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    is_packing[i] = false;

  for (int ict = 0; ict < _n_ct_zones; ict++) {
    const cs_ctwr_zone_t *ct = cts[ict];

    const int z_id = ct->z_id;
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    /* At this stage, zone are not defined contrary to the mesh location
     * So, we retrieve the mesh location information
     */
    const int  ml_id = z->location_id;
    const cs_lnum_t  _n_elts = cs_mesh_location_get_n_elts(ml_id)[0];
    const cs_lnum_t  *_elt_ids = cs_mesh_location_get_elt_ids(ml_id);

    if (_elt_ids == NULL)
      for (cs_lnum_t j = 0; j < _n_elts; j++) is_packing[j] = true;
    else
      for (cs_lnum_t j = 0; j < _n_elts; j++) is_packing[_elt_ids[j]] = true;

  }

  /* Count the number of cells attached to a packing zone */
  cs_lnum_t  n_pack_elts = 0;
  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    if (is_packing[i]) n_pack_elts++;

  cs_lnum_t *pack_elts = NULL;
  if (n_pack_elts < m->n_cells) {

    /* Fill list  */
    BFT_MALLOC(pack_elts, n_pack_elts, cs_lnum_t);

    cs_lnum_t shift = 0;
    for (cs_lnum_t i = 0; i < m->n_cells; i++)
      if (is_packing[i]) pack_elts[shift++] = i;

    assert(shift == n_pack_elts);

  } /* Build elt_ids */

  BFT_FREE(is_packing);

  /* Return pointers */
  *n_elts = n_pack_elts;
  *elt_ids = pack_elts;
}


/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide access to cs_ctwr_option
 *----------------------------------------------------------------------------*/

cs_ctwr_option_t *
cs_get_glob_ctwr_option(void)
{
  return &_ctwr_option;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_ctwr_zones
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t **
cs_get_glob_ctwr_zone(void)
{
  return _ct_zone;
}

/*----------------------------------------------------------------------------
 * Provide access to number of ct zones
 *----------------------------------------------------------------------------*/

int *
cs_get_glob_ctwr_n_zones(void)
{
  return &_n_ct_zones;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cooling tower exchange zone
 *
 * \param[in]  zone_criteria  zone selection criteria (or NULL)
 * \param[in]  z_id           z_id if zone already created (-1 otherwise)
 * \param[in]  zone_type      exchange zone type
 * \param[in]  delta_t        imposed delta temperature delta between inlet
 *                            and outlet of the zone
 * \param[in]  relax          relaxation of the imposed delta temperature
 * \param[in]  t_l_bc         liquid water temperature at the inlet
 * \param[in]  q_l_bc         mass flow rate at the inlet
 * \param[in]  xap            beta_x_0 of the exchange law
 * \param[in]  xnp            exponent n of the exchange law
 * \param[in]  surface        total Surface of ingoing water
 * \param[in]  xleak_fac      leakage factor (ratio of outlet/inlet flow rate)
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define(const char           zone_criteria[],
               int                  z_id,
               cs_ctwr_zone_type_t  zone_type,
               cs_real_t            delta_t,
               cs_real_t            relax,
               cs_real_t            t_l_bc,
               cs_real_t            q_l_bc,
               cs_real_t            xap,
               cs_real_t            xnp,
               cs_real_t            surface,
               cs_real_t            xleak_fac)
{
  cs_ctwr_zone_t  *ct;
  int length;

  /* Verify input parameters */
  bool valid = true;

  if (   zone_type != CS_CTWR_COUNTER_CURRENT
      && zone_type != CS_CTWR_CROSS_CURRENT
      && zone_type != CS_CTWR_INJECTION) {
    /* Error message */
    bft_printf("Unrecognised packing zone type. The zone type must be either: \n"
               "CS_CTWR_COUNTER_CURRENT or CS_CTWR_CROSS_CURRENT\n");
    valid = false;
  }

  if (xleak_fac > 1.0) {
    /* Error message */
    bft_printf("Out of range leak factor.  The leak factor is a percentage and"
               "must be either: \n"
               "Negative, to indicate that the packing zone does not leak, or\n"
               "Between 0 and 1 to specify the fraction of liquid mass flow rate"
               "leaking out of the zone\n");
    valid = false;
  }

  if (!valid) {
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid packing zone specification\n"
                "Verify parameters\n"));
  }

  /* Define  a new exchange zone */

  BFT_MALLOC(ct, 1, cs_ctwr_zone_t);

  ct->criteria = NULL;
  if (zone_criteria != NULL) {
    BFT_MALLOC(ct->criteria, strlen(zone_criteria)+1, char);
    strcpy(ct->criteria, zone_criteria);
  }
  ct->num = _n_ct_zones + 1;
  ct->z_id = z_id;

  ct->type = zone_type;

  ct->name = NULL;
  const cs_zone_t *z = NULL;
  if (z_id > -1) {
    z = cs_volume_zone_by_id(z_id);
    length = strlen(z->name) + 1;
    BFT_MALLOC(ct->name, length, char);
    strcpy(ct->name, z->name);
  }
  else {
    length = strlen("cooling_towers_") + 3;
    BFT_MALLOC(ct->name, length, char);
    sprintf(ct->name, "cooling_towers_%02d", ct->num);
  }
  ct->file_name = NULL;

  if (ct->type != CS_CTWR_INJECTION)
    ct->delta_t = delta_t;
  else if (ct->type == CS_CTWR_INJECTION && delta_t > 0.){
    bft_printf("WARNING: imposed temperature difference is not possible\n"
               "for injection zone. Value will not be considered.\n\n");
    ct->delta_t = -1;
  }

  ct->relax   = relax;
  ct->t_l_bc  = t_l_bc;
  ct->q_l_bc  = q_l_bc;

  ct->xap = xap;
  ct->xnp = xnp;

  ct->surface_in  = 0.;
  ct->surface_out = 0.;
  ct->surface = surface;

  ct->xleak_fac = xleak_fac;
  ct->v_liq_pack = 0.1; /* Usual value of liquid film velocity in packing,
                           see Jourdan et al. 2022 IEC research */

  ct->n_cells = 0;

  ct->up_ct_id = -1;

  ct->n_inlet_faces = 0;
  ct->n_outlet_faces = 0;
  ct->inlet_faces_ids = NULL;
  ct->outlet_faces_ids = NULL;

 /* Different from number of faces if split faces on cells
    Can not allow non-conformal or there could be a mix up between leaking and
    non-leaking zones. */

  ct->n_outlet_cells = 0;
  ct->outlet_cells_ids = NULL;

  ct->q_l_in = 0.0;
  ct->q_l_out = 0.0;
  ct->t_l_in = 0.0;
  ct->t_l_out = 0.0;
  ct->h_l_in = 0.0;
  ct->h_l_out = 0.0;
  ct->t_h_in = 0.0;
  ct->t_h_out = 0.0;
  ct->xair_e = 0.0; //FIXME useless ?
  ct->xair_s = 0.0;
  ct->h_h_in = 0.0;
  ct->h_h_out = 0.0;
  ct->q_h_in = 0.0;
  ct->q_h_out = 0.0;

  if (_n_ct_zones >= _n_ct_zones_max) {
    _n_ct_zones_max = (_n_ct_zones_max + 1);
    BFT_REALLOC(_ct_zone, _n_ct_zones_max, cs_ctwr_zone_t *);
  }

  /* Add it to exchange zones array */

  _ct_zone[_n_ct_zones] = ct;
  _n_ct_zones += 1;

  if (cs_glob_rank_id <= 0) {
    length = strlen("cooling_towers_balance.") + 2 + 1;
    for (int _num = ct->num; _num > 99; _num /= 10)
      length += 1;
    BFT_MALLOC(ct->file_name, length, char);
    sprintf(ct->file_name, "cooling_towers_balance.%02d", ct->num);

    FILE *f = fopen(ct->file_name, "a");

    fprintf(f, "# Balance for the exchange zone %02d\n", ct->num);
    fprintf(f, "# ================================\n");
    fprintf(f, "# Time  Flux air/liq");
    fprintf(f, "\tTemp liq in");
    fprintf(f, "\tTemp liq out");
    fprintf(f, "\tTemp air in");
    fprintf(f, "\tTemp air out");
    fprintf(f, "\tFlow liq in\tFlow liq out");
    fprintf(f, "\tFlow air in\tFlow air out");
    fprintf(f, "\tPressure in\tPressure out\n");
    fclose(f);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define cooling tower zones.
 *
 * TODO rename this: definition (at setup stage) and build (instanciation on
 *      actual mesh are not the same).
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define_zones(void)
{
  /* Check if there are any leaking packing zones, if yes, there is rain */
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
  for (int ict = 0; ict < _n_ct_zones && !(ct_opt->has_rain); ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];
    if (ct->xleak_fac > 0.0)
      ct_opt->has_rain = true;
  }

  /* Define the zones with source terms */
  if (ct_opt->has_rain) {
    /* Phase change may take place in the entire computational domain
     * so activate mass source term to zone 0 */
    cs_volume_zone_set_type(0, CS_VOLUME_ZONE_MASS_SOURCE_TERM);

    /* Identify cooling towers zones for cs_ctwr_build_all
       but don't redeclare the cells as mass_source_term
       to avoid double counting */
    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      int z_id = ct->z_id;
      if (z_id > -1)
        cs_volume_zone_set_type(z_id, CS_VOLUME_ZONE_INITIALIZATION);
      else {
        z_id = cs_volume_zone_define(ct->name,
                                     ct->criteria,
                                     CS_VOLUME_ZONE_INITIALIZATION);
        ct->z_id = z_id;
      }
    }
  }
  else {
    /* Phase change will take place only in the packing zones */
    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      int z_id = ct->z_id;
      if (z_id > -1)
        cs_volume_zone_set_type(z_id, CS_VOLUME_ZONE_MASS_SOURCE_TERM);
      else {
        z_id = cs_volume_zone_define(ct->name,
                                     ct->criteria,
                                     CS_VOLUME_ZONE_MASS_SOURCE_TERM);
        ct->z_id = z_id;
      }
    }
  }

  /* Define the packing zone (union of all packings), "auto:packings" */
  if (_n_ct_zones > 0){
    const char  zone_name[] = "auto:packings";
    int z_id = cs_volume_zone_define_by_func(zone_name,
                                             _packing_selection,
                                             _ct_zone, /* input */
                                             0); /* flag */

    cs_volume_zone_set_overlay(z_id, true);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different packing zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_all(void)
{
  /* Loop over exchange zones: set number of cells */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    /* Set number of cells */
    const cs_zone_t *z = cs_volume_zone_by_name(ct->name);
    ct->n_cells = z->n_elts;
    ct->vol_f = z->f_measure;

    void *input = (void *) ct;
    cs_equation_param_t *eqp =
        cs_field_get_equation_param(CS_F_(p));

    /* Set the bulk mass source term function (associated to pressure field) */
    cs_equation_add_volume_mass_injection_by_dof_func
      (eqp,
       z->name,
       cs_flag_primal_cell,
       cs_ctwr_volume_mass_injection_dof_func,
       input);

    eqp = cs_field_get_equation_param(CS_F_(ym_w));

    /* Set value of ingoing water */
    cs_real_t yw_in = 1.;
    cs_equation_add_volume_mass_injection_by_value(eqp, z->name, &yw_in);

    /* Injection zone */
    if (ct->xleak_fac > 0.0 && ct->type == CS_CTWR_INJECTION) {

      /* Rain mass fraction */
      cs_field_t *f_yp = cs_field_by_name("ym_l_r");
      eqp = cs_field_get_equation_param(f_yp);

      /* Set value of ingoing rain */
      cs_real_t y_in = 1.;
      cs_equation_add_volume_mass_injection_by_value(eqp, z->name, &y_in);

      /* Rain enthalpy */
      cs_field_t *f_yh_rain = cs_field_by_name("ymh_l_r"); /* Yp times Tp */
      eqp = cs_field_get_equation_param(f_yh_rain);
      cs_real_t t_in = ct->t_l_bc;

      // FIXME: There should be a y_p factor in there so that
      // mass and enthalpy are compatible
      /* The transported variable is y_rain * h_rain */
      cs_real_t h_in = cs_liq_t_to_h(t_in);
      cs_equation_add_volume_mass_injection_by_value(eqp, z->name, &h_in);

    }


  }
  /* Define the zones with source terms */
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
  if (ct_opt->has_rain) {
    /* Set number of cells */
    const cs_zone_t *z = cs_volume_zone_by_id(0);

    cs_equation_param_t *eqp =
        cs_field_get_equation_param(CS_F_(p));

    cs_equation_add_volume_mass_injection_by_dof_func
      (eqp,
       z->name,
       cs_flag_primal_cell,
       cs_ctwr_volume_mass_injection_rain_dof_func,
       NULL);

    /* Rain mass fraction */
    cs_field_t *f_yp = cs_field_by_name("ym_l_r");
    eqp = cs_field_get_equation_param(f_yp);

    /* Set value of ingoing rain */
    cs_real_t y_in = 1.;
    cs_equation_add_volume_mass_injection_by_value(eqp, z->name, &y_in);

    /* Rain enthalpy (yp.hp) */
    cs_field_t *f_yphp = cs_field_by_name("ymh_l_r");
    eqp = cs_field_get_equation_param(f_yphp);
    cs_equation_add_volume_mass_injection_by_dof_func(eqp,
                                                      z->name,
                                                      cs_flag_primal_cell,
                                                      cs_ctwr_volume_mass_injection_yh_rain_dof_func,
                                                      NULL);

  }

  /* Post-processing: multiply enthalpy by fraction */

  cs_field_t *f = cs_field_by_name_try("enthalpy_liquid");
  if (f != NULL) {
    const int vis_key_id = cs_field_key_id("post_vis");
    if (cs_field_get_key_int(f, vis_key_id) & CS_POST_ON_LOCATION) {
      cs_post_add_time_mesh_dep_output(_write_liquid_vars, NULL);
      cs_field_clear_key_int_bits(f, vis_key_id, CS_POST_ON_LOCATION);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy cs_ctwr_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_all_destroy(void)
{
  for (int id = 0; id < _n_ct_zones; id++) {

    cs_ctwr_zone_t  *ct = _ct_zone[id];
    BFT_FREE(ct->criteria);
    BFT_FREE(ct->name);
    BFT_FREE(ct->file_name);
    BFT_FREE(ct->inlet_faces_ids);
    BFT_FREE(ct->outlet_faces_ids);
    BFT_FREE(ct->outlet_cells_ids);
    BFT_FREE(ct);

  }

  _n_ct_zones_max = 0;
  _n_ct_zones = 0;

  BFT_FREE(_ct_zone);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log Packing zone definition setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_setup(void)
{
  if (_n_ct_zones < 1)
    return;

  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

  /* Verify the input parameters */
  if (ct_opt->evap_model != CS_CTWR_NONE
      && ct_opt->evap_model != CS_CTWR_POPPE
      && ct_opt->evap_model != CS_CTWR_MERKEL) {

    bft_printf("Unrecognised evaporation model. "
               "The evaporation model must be either:\n"
               "CS_CTWR_NONE or CS_CTWR_POPPE or CS_CTWR_MERKEL\n");
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid evaporation model specification\n"
                "Verify parameters\n"));
  }

  const char *model_type_name[] = {"None", "Poppe", "Merkel"};

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Cooling towers\n"
                  "--------------\n"
                  "  Droplet diameter: %f\n"
                  "  Evaporation model: %s\n"),
                cs_glob_air_props->droplet_diam,
                model_type_name[ct_opt->evap_model]);

  for (int i = 0; i < _n_ct_zones; i++) {
    cs_ctwr_zone_t *ct = _ct_zone[i];

    if (ct->criteria != NULL)
      cs_log_printf
        (CS_LOG_SETUP,
         _("  Cooling tower zone num: %d\n"
           "    zone id: %d\n"
           "    criterion: ""%s""\n"
           "    Parameters:\n"
           "      Lambda of the exchange law: %f\n"
           "      Exponent n of the exchange law: %f\n"
           "      Type: %d\n"
           "      Delta Temperature: %f\n"
           "        Relaxation: %f\n"
           "      Injected water temperature: %f\n"
           "      Injected mass flow rate: %f\n"
           "      Total surface of ingoing water: %f\n"),
         ct->num,
         ct->z_id,
         ct->criteria,
         ct->xap,
         ct->xnp,
         ct->type,
         ct->delta_t,
         ct->relax,
         ct->t_l_bc,
         ct->q_l_bc,
         ct->surface);
    else
      cs_log_printf
        (CS_LOG_SETUP,
         _("  Cooling tower num: %d\n"
           "    zone id: %d\n"
           "    Parameters:\n"
           "      Lambda of the exchange law: %f\n"
           "      Exponent n of the exchange law: %f\n"
           "      Type: %d\n"
           "      Delta Temperature: %f\n"
           "        Relaxation: %f\n"
           "      Injected water temperature: %f\n"
           "      Injected mass flow rate: %f\n"
           "      Total surface of ingoing water: %f\n"),
         ct->num,
         ct->z_id,
         ct->xap,
         ct->xnp,
         ct->type,
         ct->delta_t,
         ct->relax,
         ct->t_l_bc,
         ct->q_l_bc,
         ct->surface);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform balances in packing zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_log_balance(void)
{
  //TODO : Separate log depending on zone type (exchange or injection)
  if (_n_ct_zones < 1)
    return;

  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->i_face_normal;
  cs_real_t *p = (cs_real_t *)CS_F_(p)->val;        /* Pressure */
  cs_real_t *t_h = NULL;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID)
    t_h = cs_field_by_name("real_temperature")->val; /* Humid air temp */
  else
    t_h = cs_field_by_name("temperature")->val; /* Humid air temp */

  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;      /* Humid air enthalpy */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l_pack)->val;    /* Liquid temperature */
  cs_real_t *yh_l_p = (cs_real_t *)CS_F_(yh_l_pack)->val;    /* Liquid enthalpy */
  cs_real_t *y_l_p = (cs_real_t *)CS_F_(y_l_pack)->val;   /* Liquid mass per unit
                                                           cell volume */

  // FIXME take the good one... for y_p
  cs_real_t *liq_mass_flow
    = cs_field_by_name("inner_mass_flux_y_l_packing")->val;
  cs_real_t *mass_flow = cs_field_by_name("inner_mass_flux")->val;

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    ct->p_in = 0.0;
    ct->p_out = 0.0;
    ct->q_l_in = 0.0;
    ct->q_l_out = 0.0;
    ct->t_l_in = 0.0;
    ct->h_l_out = 0.0;
    ct->h_l_in = 0.0;
    ct->t_l_out = 0.0;
    ct->t_h_in = 0.0;
    ct->t_h_out = 0.0;
    ct->xair_e = 0.0;
    ct->xair_s = 0.0;
    ct->h_h_in = 0.0;
    ct->h_h_out = 0.0;
    ct->q_h_in = 0.0;
    ct->q_h_out = 0.0;

    /* Compute liquid water quantities
     * And humid air quantities at liquid inlet */
    for (cs_lnum_t i = 0; i < ct->n_inlet_faces; i++) {

      cs_lnum_t face_id = ct->inlet_faces_ids[i];
      cs_lnum_t cell_id_l, cell_id_h;
      cs_real_t face_surf = cs_math_3_norm(i_face_normal[face_id]);

      /* Convention: inlet is negative mass flux
       * Then upwind cell for liquid is i_face_cells[][1] */
      int sign = 1;
      if (liq_mass_flow[face_id] > 0) {
        sign = -1;
        cell_id_l = i_face_cells[face_id][0];
        cell_id_h = i_face_cells[face_id][1];
      }
      else {
        cell_id_l = i_face_cells[face_id][1];
        cell_id_h = i_face_cells[face_id][0];
      }

      /* Liquid inlet = air outlet -> outlet pressure */
      ct->p_out += p[cell_id_h] * face_surf;
      /* (y_l. h_l) is transported with (rho u_l)
       * so h_l is transported with (y_l rho u_l) */
      ct->t_l_in += sign * t_l[cell_id_l]
        * y_l_p[cell_id_l] * liq_mass_flow[face_id];
      ct->h_l_in += sign * yh_l_p[cell_id_l] * liq_mass_flow[face_id];
      ct->q_l_in += sign * y_l_p[cell_id_l] * liq_mass_flow[face_id];

      ct->t_h_out += sign * t_h[cell_id_h] * mass_flow[face_id];
      ct->h_h_out += sign * h_h[cell_id_h] * mass_flow[face_id];
      ct->q_h_out += sign * mass_flow[face_id];

      //ct->xair_s  += debit*xa[icel];
    }
    double stmp[7] = {ct->t_l_in, ct->h_l_in, ct->q_l_in,
      ct->t_h_out, ct->h_h_out, ct->q_h_out, ct->p_out};

    cs_parall_sum(7, CS_DOUBLE, stmp);

    ct->t_l_in = stmp[0]; ct->h_l_in = stmp[1]; ct->q_l_in = stmp[2];
    ct->t_h_out = stmp[3]; ct->h_h_out = stmp[4]; ct->q_h_out = stmp[5];
    ct->p_out = stmp[6];

    ct->t_l_in /= ct->q_l_in;
    ct->h_l_in /= ct->q_l_in;
    ct->q_l_in /= ct->surface_in;
    ct->p_out /= ct->surface_in;

    if (CS_ABS(ct->q_h_out) > 1e-10) {
      ct->t_h_out /= ct->q_h_out;
      ct->h_h_out /= ct->q_h_out;
    }
    ct->q_h_out /= ct->surface_in;

    /* Compute liquid water quantities
     * And humid air quantities at liquid packing outlet  */
    for (cs_lnum_t i = 0; i < ct->n_outlet_faces; i++) {

      cs_lnum_t face_id = ct->outlet_faces_ids[i];
      cs_lnum_t cell_id_l, cell_id_h;
      cs_real_t face_surf = cs_math_3_norm(i_face_normal[face_id]);

      /* Convention: outlet is positive mass flux
       * Then upwind cell for liquid is i_face_cells[][0] */
      int sign = 1;
      if (liq_mass_flow[face_id] < 0) {
        sign = -1;
        cell_id_l = i_face_cells[face_id][1];
        cell_id_h = i_face_cells[face_id][0];
      }
      else {
        cell_id_l = i_face_cells[face_id][0];
        cell_id_h = i_face_cells[face_id][1];
      }

      /* Liquid outlet = air inlet -> inlet pressure */
      ct->p_in += p[cell_id_h] * face_surf;
      /* h_l is in fact (y_l. h_l),
       * and the transport field is (y_l*liq_mass_flow) */
      ct->t_l_out += sign * t_l[cell_id_l]
        * y_l_p[cell_id_l] * liq_mass_flow[face_id];
      ct->q_l_out += sign * y_l_p[cell_id_l] * liq_mass_flow[face_id];
      ct->h_l_out += sign * yh_l_p[cell_id_l] * liq_mass_flow[face_id];

      // FIXME: Sign coming from liq_mass_flow
      // and applied to mass_flow - correct?
      ct->t_h_in  += sign * t_h[cell_id_h] * mass_flow[face_id];
      ct->h_h_in  += sign * h_h[cell_id_h] * mass_flow[face_id];
      ct->q_h_in  += sign * mass_flow[face_id];
    }

    cs_parall_sum(1, CS_REAL_TYPE, &(ct->t_l_out));
    cs_parall_sum(1, CS_REAL_TYPE, &(ct->q_l_out));
    cs_parall_sum(1, CS_REAL_TYPE, &(ct->h_l_out));

    cs_parall_sum(1, CS_REAL_TYPE, &(ct->t_h_in));
    cs_parall_sum(1, CS_REAL_TYPE, &(ct->h_h_in));
    cs_parall_sum(1, CS_REAL_TYPE, &(ct->q_h_in));

    cs_parall_sum(1, CS_REAL_TYPE, &(ct->p_in));

    ct->t_l_out /= ct->q_l_out;
    ct->h_l_out /= ct->q_l_out;
    ct->q_l_out /= ct->surface_out;
    ct->p_in /= ct->surface_out;

    if (CS_ABS(ct->q_h_in) > 1e-10) {
      ct->t_h_in /= ct->q_h_in;
      ct->h_h_in /= ct->q_h_in;
    }
    ct->q_h_in /= ct->surface_out;

    /* Writings */
    if (cs_glob_rank_id <= 0) {
      if (CS_ABS(ct->h_l_in - ct->h_l_out)> 1.e-6) {
        FILE *f = fopen(ct->file_name, "a");
        cs_real_t aux = cs_math_fabs(  (ct->h_h_out - ct->h_h_in)
            / (ct->h_l_in - ct->h_l_out));
        fprintf(f,
            "%10f\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t"
            "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n",
            cs_glob_time_step->t_cur,
            aux,
            ct->t_l_in,
            ct->t_l_out,
            ct->t_h_in,
            ct->t_h_out,
            ct->q_l_in,
            ct->q_l_out,
            ct->q_h_in,
            ct->q_h_out,
            ct->p_in,
            ct->p_out);
        fclose(f);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

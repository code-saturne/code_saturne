/*============================================================================
 *
 * Definitions, Global variables variables, and functions associated with the
 * exchange zones
 *============================================================================*/

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_ctwr_air_props.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr.h"

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
  .has_rain = false};

/* Cooling tower exchange zone structure definition */
/*--------------------------------------------------*/

struct _cs_ctwr_zone_t {

  int                  num;        /* Exchange zone number */
  char                *criteria;   /* Exchange zone selcection criteria */
  char                *name;       /* Exchange zone name */
  char                *file_name;  /* Exchange zone budget file name */
  cs_ctwr_zone_type_t  type;       /* Zone type */

  cs_real_t  hmin;               /* Minimum vertical height of exchange zone */
  cs_real_t  hmax;               /* Maximum height of exchange zone */
  cs_real_t  delta_t;            /* Temperature delta required for exchange zone
                                    if positive */
  cs_real_t  relax;              /* Relaxation of the imposed temperature */

  cs_real_t  t_l_bc;             /* Water entry temperature */
  cs_real_t  q_l_bc;             /* Water flow */
  cs_real_t  y_l_bc;             /* Mass fraction of water */

  cs_real_t  xap;                /* Exchange law a_0 coefficient */
  cs_real_t  xnp;                /* Exchange law n exponent */

  cs_real_t  surface_in;         /* Water inlet surface */
  cs_real_t  surface_out;        /* Water outlet surface */
  cs_real_t  surface;            /* Total surface */

  cs_real_t  xleak_fac;          /* Leakage factor (ratio of outlet/inlet
                                    flow rate) */

  cs_int_t   n_cells;            /* Number of air cells belonging to the zone */

  cs_int_t   up_ct_id;           /* Id of upstream exchange zone (if any) */

  cs_lnum_t  n_inlet_faces;      /* Number of inlet faces */
  cs_lnum_t  n_outlet_faces;     /* Number of outlet faces */
  cs_lnum_t *inlet_faces_ids;    /* List of inlet faces */
  cs_lnum_t *outlet_faces_ids;   /* List of outlet faces */

  cs_lnum_t  n_outlet_cells;     /* Number of outlet cells */
  cs_lnum_t *outlet_cells_ids;   /* List of outlet cells */

  cs_real_t  q_l_in;          /* Water entry flow */
  cs_real_t  q_l_out;         /* Water exit flow */
  cs_real_t  t_l_in;          /* Mean water entry temperature */
  cs_real_t  t_l_out;         /* Mean water exit temperature */
  cs_real_t  h_l_in;          /* Mean water entry enthalpy */
  cs_real_t  h_l_out;         /* Mean water exit enthalpy */
  cs_real_t  t_h_in;          /* Mean air entry temperature */
  cs_real_t  t_h_out;         /* Mean air exit temperature */
  cs_real_t  xair_e;          /* Mean air entry humidity */
  cs_real_t  xair_s;          /* Mean air exit humidity */
  cs_real_t  h_h_in;          /* Mean air entry enthalpy */
  cs_real_t  h_h_out;         /* Mean air exit enthalpy */
  cs_real_t  q_h_in;          /* Air entry flow */
  cs_real_t  q_h_out;         /* Air exit flow */

};

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

    cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;   /* liquid enthalpy */
    cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;  /* liquid mass per unit
                                                        cell volume*/

    cs_real_t *val;
    BFT_MALLOC(val, mesh->n_cells, cs_real_t);

    /* Value on all cells */

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      val[i] = 0;

    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->cell_ids;
      for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
        cs_lnum_t cell_id = ze_cell_ids[i];
        if (y_l[cell_id] > 0.)
          val[cell_id] = h_l[cell_id]/y_l[cell_id];
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

    const char name[] = "Liquid fraction enthalpy";

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Provide acces to cs_ctwr_option
 *----------------------------------------------------------------------------*/

cs_ctwr_option_t *
cs_get_glob_ctwr_option(void)
{
  return &_ctwr_option;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cooling tower exchange zone
 *
 * \param[in]   zone_criteria   Zone selection criteria
 * \param[in]   zone_type       exchange zone type
 * \param[in]   delta_t         Imposed delta temperature delta between inlet
 *                              and oulet of the zone
 * \param[in]   relax           Relaxation of the imposed delta temperature
 * \param[in]   t_l_bc          Liquid water temperature at the inlet
 * \param[in]   q_l_bc          Mass flow rate at the inlet
 * \param[in]   xap             Lambda of the exchange law
 * \param[in]   xnp             Exponent n of the exchange law
 * \param[in]   surface         Total Surface of ingoing water
 * \param[in]   xleak_fac       Leakage factor (ratio of outlet/inlet flow rate)
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define(const char           zone_criteria[],
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
  FILE *f;

  /* Define  a new exchange zone */

  BFT_MALLOC(ct, 1, cs_ctwr_zone_t);

  ct->criteria = NULL;
  BFT_MALLOC(ct->criteria, strlen(zone_criteria)+1, char);
  strcpy(ct->criteria, zone_criteria);

  ct->num = _n_ct_zones + 1;

  ct->type = zone_type;

  ct->name = NULL;
  length = strlen("cooling_towers_") + 3;
  BFT_MALLOC(ct->name, length, char);
  sprintf(ct->name, "cooling_towers_%02d", ct->num);

  ct->file_name = NULL;

  ct->delta_t = delta_t;
  ct->relax   = relax;
  ct->t_l_bc  = t_l_bc;
  ct->q_l_bc  = q_l_bc;
  ct->y_l_bc  = -1; /* Mass of liquid water divided by the mass of humid air
                        in packing zones.
                        Factice version, initialize after */
  ct->xap = xap;
  ct->xnp = xnp;

  ct->surface_in  = 0.;
  ct->surface_out = 0.;
  ct->surface = surface;

  ct->xleak_fac = xleak_fac;

  ct->n_cells = 0;

  ct->up_ct_id = -1;

  ct->n_inlet_faces = 0;
  ct->n_outlet_faces = 0;
  ct->inlet_faces_ids = NULL;
  ct->outlet_faces_ids = NULL;

 /* Different from number of faces if split faces on cells
    Can not allow non-conformal or there could be a mix up between leaking and
    non-leaking zones
  */
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
  ct->xair_e = 0.0;//FIXME useless ?
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

    f = fopen(ct->file_name, "a");

    fprintf(f, "# Balance for the exchange zone %02d\n", ct->num);
    fprintf(f, "# ==========================================================\n");
    fprintf(f, "\tTime\tFlux air/liq");
    fprintf(f, "\tTemp liq in");
    fprintf(f, "\tTemp liq out");
    fprintf(f, "\tTemp air in");
    fprintf(f, "\tTemp air out");
    fprintf(f, "\tDeb liq in\tDeb liq out");
    fprintf(f, "\tDeb air in\tDeb air out\n");
    fclose(f);
  }


}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map fields used by the cooling tower module to pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_field_pointer_map(void)
{
  /* No need to redefine the temperature and enthalpy for humid air as they
     have already been defined in 'cs_field_pointer_map',
     which comes after 'ctvarp' */
  cs_field_pointer_map(CS_ENUMF_(humid), cs_field_by_name_try("humidity"));
  cs_field_pointer_map(CS_ENUMF_(ym_w), cs_field_by_name_try("ym_water"));
  cs_field_pointer_map(CS_ENUMF_(t_l), cs_field_by_name_try("temperature_liquid"));
  cs_field_pointer_map(CS_ENUMF_(h_l), cs_field_by_name_try("enthalpy_liquid"));
  cs_field_pointer_map(CS_ENUMF_(y_l_pack), cs_field_by_name_try("y_l_packing"));
  cs_field_pointer_map(CS_ENUMF_(thermal_diff_h),
                       cs_field_by_name_try("thermal_conductivity"));
}

/*----------------------------------------------------------------------------*/

/*!
 * \brief  Define zones.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_zones(void)
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
    /* Phase change may take place in the entire computational domain */
    void *input;
    cs_volume_zone_define("rain_zone",
                          "all[]",
                          CS_VOLUME_ZONE_MASS_SOURCE_TERM);

    /* Identify packing zones for cs_ctwr_build_all
       but don't redeclare the cells as mass_source_term to avoid double counting */
    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      cs_volume_zone_define(ct->name,
                            ct->criteria,
                            CS_VOLUME_ZONE_INITIALIZATION);
    }
  } else {
    /* Phase change will  take place only in the packing zones */
    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      cs_volume_zone_define(ct->name,
                            ct->criteria,
                            CS_VOLUME_ZONE_MASS_SOURCE_TERM);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different packing zones.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_all(void)
{
  /* Loop over exchange zones: set number of cells */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    /* Set number of cells */
    ct->n_cells = cs_volume_zone_by_name(ct->name)->n_cells;
  }

  /* Postprocessing: multiply enthalpy by fraction */

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
  cs_ctwr_zone_t  *ct;

  for (int id = 0; id < _n_ct_zones; id++) {

    ct = _ct_zone[id];
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

  char model_type[16];
  if (ct_opt->evap_model == CS_CTWR_NONE) {
    snprintf(model_type, 15, "None");
  } else if (ct_opt->evap_model == CS_CTWR_POPPE) {
    snprintf(model_type, 15, "Poppe");
  } else if (ct_opt->evap_model == CS_CTWR_MERKEL) {
    snprintf(model_type, 15, "Merkel");
  }

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Cooling towers\n"
                  "--------------\n"
                  "  Droplet diameter: %f\n"
                  "  Evaporation model: %s\n"),
                cs_glob_ctwr_props->droplet_diam,
                model_type
                );

  for (int i = 0; i < _n_ct_zones; i++) {
    cs_ctwr_zone_t *ct = _ct_zone[i];

    cs_log_printf
      (CS_LOG_SETUP,
       _("  Cooling tower zone id: %d\n"
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
       ct->criteria,
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
  if (_n_ct_zones < 1)
    return;

  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);
  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;  /* humid air bulk density */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;      /* humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;      /* humid air enthalpy */
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val;   /* Water mass fraction
                                                       in humid air */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;    /* humid air bulk humidity */

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /* liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /* liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;   /* liquid mass per unit
                                                       cell volume */

  cs_real_t *liq_mass_flow = cs_field_by_name("inner_mass_flux_y_l_packing")->val;//FIXME take the good one... for y_p
  cs_real_t *mass_flow = cs_field_by_name("inner_mass_flux")->val;

  FILE *f;

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

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

      /* Convention: inlet is negativ mass flux
       * Then upwind cell for liquid is i_face_cells[][1] */
      int sign = 1;
      if (liq_mass_flow[face_id] > 0) {
        sign = -1;
        cell_id_l = i_face_cells[face_id][0];
        cell_id_h = i_face_cells[face_id][1];
      } else {
        cell_id_l = i_face_cells[face_id][1];
        cell_id_h = i_face_cells[face_id][0];
      }

      /* (y_l. h_l) is transported with (rho u_l)
       * so h_l is transported with (y_l rho u_l) */
      ct->t_l_in += sign * t_l[cell_id_l]
                         * y_l[cell_id_l] * liq_mass_flow[face_id];
      ct->h_l_in += sign * h_l[cell_id_l] * liq_mass_flow[face_id];
      ct->q_l_in += sign * y_l[cell_id_l] * liq_mass_flow[face_id];

      ct->t_h_out += sign * t_h[cell_id_h] * mass_flow[face_id];
      ct->h_h_out += sign * h_h[cell_id_h] * mass_flow[face_id];
      ct->q_h_out += sign * mass_flow[face_id];

      //ct->xair_s  += debit*xa[icel];
    }

    double stmp[6] = {ct->t_l_in, ct->h_l_in, ct->q_l_in,
                      ct->t_h_out, ct->h_h_out, ct->q_h_out};

    cs_parall_sum(6, CS_DOUBLE, stmp);

    ct->t_l_in = stmp[0]; ct->h_l_in = stmp[1]; ct->q_l_in = stmp[2];
    ct->t_h_out = stmp[3]; ct->h_h_out = stmp[4]; ct->q_h_out = stmp[5];

    ct->t_l_in /= ct->q_l_in;
    ct->h_l_in /= ct->q_l_in;
    ct->q_l_in /= ct->surface_in;

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

      /* Convention: outlet is positive mass flux
       * Then upwind cell for liquid is i_face_cells[][0] */
      int sign = 1;
      if (liq_mass_flow[face_id] < 0) {
        sign = -1;
        cell_id_l = i_face_cells[face_id][1];
        cell_id_h = i_face_cells[face_id][0];
      } else {
        cell_id_l = i_face_cells[face_id][0];
        cell_id_h = i_face_cells[face_id][1];
      }

      /* h_l is in fact (y_l. h_l),
       * and the transport field is (y_l*liq_mass_flow) */
      ct->t_l_out += sign * t_l[cell_id_l]
        * y_l[cell_id_l] * liq_mass_flow[face_id];
      ct->q_l_out += sign * y_l[cell_id_l] * liq_mass_flow[face_id];
      ct->h_l_out += sign * h_l[cell_id_l] * liq_mass_flow[face_id];

      ct->t_h_in  += sign * t_h[cell_id_h] * mass_flow[face_id];
      ct->h_h_in  += sign * h_h[cell_id_h] * mass_flow[face_id];
      ct->q_h_in  += sign * mass_flow[face_id];
    }

    cs_parall_sum(1, CS_DOUBLE, &(ct->t_l_out));
    cs_parall_sum(1, CS_DOUBLE, &(ct->q_l_out));
    cs_parall_sum(1, CS_DOUBLE, &(ct->h_l_out));

    cs_parall_sum(1, CS_DOUBLE, &(ct->t_h_in));
    cs_parall_sum(1, CS_DOUBLE, &(ct->h_h_in));
    cs_parall_sum(1, CS_DOUBLE, &(ct->q_h_in));

    ct->t_l_out /= ct->q_l_out;
    ct->h_l_out /= ct->q_l_out;
    ct->q_l_out /= ct->surface_out;

    if (CS_ABS(ct->q_h_in) > 1e-10) {
      ct->t_h_in /= ct->q_h_in;
      ct->h_h_in /= ct->q_h_in;
    }
    ct->q_h_in /= ct->surface_out;

    /* Writings */
    if (cs_glob_rank_id <= 0) {
      if (CS_ABS(ct->h_l_in - ct->h_l_out)> 1.e-6) {
        f = fopen(ct->file_name, "a");
        cs_real_t aux = CS_ABS(  (ct->h_h_out - ct->h_h_in)
                               / (ct->h_l_in - ct->h_l_out));
        fprintf(f,
                "%10f\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t"
                "%12.5e\t%12.5e\t%12.5e\t%12.5e\n",
                cs_glob_time_step->t_cur,
                aux,
                ct->t_l_in,
                ct->t_l_out,
                ct->t_h_in,
                ct->t_h_out,
                ct->q_l_in,
                ct->q_l_out,
                ct->q_h_in,
                ct->q_h_out);
        fclose(f);
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the field variables
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_init_field_vars(cs_real_t  rho0,
                        cs_real_t  t0,
                        cs_real_t  p0,
                        cs_real_t  molmassrat)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_halo_t *halo = m->halo;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;

  // Initialise the fields - based on map
  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;   /* humid air (bulk) density */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* humid air temperature */
  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre; /* humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* humid air enthalpy */
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val;    /* water mass fraction in
                                                        humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val; /* humidity in humid air (bulk) */

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;     /* liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;     /* liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;  /* liquid mass per unit */

  /* Packing zone liquidus vertical velocity component */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;

  /* Rain drops variables */
  cs_field_t *cfld_yp = cs_field_by_name_try("y_p");         /* Rain drops mass fraction */
  cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_p");
  cs_field_t *cfld_drift_vel = cs_field_by_name_try("drift_vel_y_p");

  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

  cs_real_t *cpro_taup = NULL;
  if (cfld_taup != NULL)
    cpro_taup = cfld_taup->val;
  else
    BFT_MALLOC(cpro_taup, n_cells_with_ghosts, cs_real_t);

  const cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;
  cs_real_t rho_l = ct_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t droplet_diam = ct_prop->droplet_diam;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field in case users have updated the initial
       dry air mass fraction.
       Note: this is a bit dubious as users could also have chosen
       to reset the humidity ? */

    /* Clippings of water mass fraction */
    if (y_w[cell_id] < 0.0)
      y_w[cell_id] = 0; //TODO count it

    if (y_w[cell_id] >= 1.0)
      y_w[cell_id] = 1. - cs_math_epzero; //TODO count it

    /* Note: the drops contribute to the bulk density */
    x[cell_id] = y_w[cell_id]/(1.0-y_w[cell_id]);

    /* Bulk humid air temperature */
    t_h[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
    t_h_a[cell_id] = t_h[cell_id];

    /* Update the humid air density */
    rho_h[cell_id] = cs_ctwr_rho_humidair(x[cell_id],
                                          rho0,
                                          p0,
                                          t0,
                                          molmassrat,
                                          t_h[cell_id]);

    /* Update the humid air enthalpy */
    x_s[cell_id] = cs_ctwr_xsath(t_h[cell_id],p0);
    cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_ctwr_h_humidair(cp_h,
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    /* Initialise the liquid vertical velocity component
     * this is correct for droplet and extended for other packing zones
     * NB: this value is derived from the drag coef:
     * C_D = 24 / Re * (1 + 0.15 * Re^0.687)
     * See ZOPLU HT-31-08-06 */

    cs_real_t v_lim = pow(droplet_diam, 2.) * rho_l / (18. * visc)
                    * cs_math_3_norm(gravity);

    cs_real_t reynolds_old = 0.;
    cs_real_t reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;

    for (int sweep = 0; sweep < 100 && CS_ABS(reynolds - reynolds_old) > 0.001; sweep++) {
      reynolds_old = reynolds;
      v_lim = pow(droplet_diam, 2.) * rho_l / (18. * visc * (1 + 0.15 * pow(reynolds, 0.687)))
            * cs_math_3_norm(gravity);
      reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
    }

    cpro_taup[cell_id] = v_lim / cs_math_3_norm(gravity);

    /* Initialize rain variables (note that Yp is already set to 0) */
    if (ct_opt->has_rain) {
      cs_real_3_t *drift_vel = (cs_real_3_t *restrict)(cfld_drift_vel->val);
      drift_vel[cell_id][0] = cpro_taup[cell_id] * gravity[0];
      drift_vel[cell_id][1] = cpro_taup[cell_id] * gravity[1];
      drift_vel[cell_id][2] = cpro_taup[cell_id] * gravity[2];
    }

  }

  /* Loop over exchange zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->cell_ids;
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Initialize with the injection water temperature */
      t_l[cell_id] = ct->t_l_bc;

      /* Update the injected liquid enthalpy */
      h_l[cell_id] = cs_ctwr_h_liqwater(t_l[cell_id]);

      /* Initialise the liquid vertical velocity component
       * this is correct for droplet and extended for other packing zones */
      vel_l[cell_id] = cpro_taup[cell_id] * cs_math_3_norm(gravity);

      /* Note that rho_h * Y_l * vel_l * Stot = q_l_bc */
      ct->y_l_bc = ct->q_l_bc / (rho_h[cell_id] * vel_l[cell_id] * ct->surface);

      /* Initialise the liquid transported variables:
         liquid mass and enthalpy corrected by the density ratio */
      y_l[cell_id] = ct->y_l_bc;

      /* The transported value is (y_l.h_l) and not (h_l) */
      h_l[cell_id] *= y_l[cell_id];

    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, vel_l);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_taup);
    if (cfld_yp != NULL)
      cs_halo_sync_var(halo, CS_HALO_STANDARD, cfld_yp->val);
    if (cfld_drift_vel != NULL) {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, cfld_drift_vel->val, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, CS_HALO_STANDARD, cfld_drift_vel->val);
    }
  }

  /* Free memory */
  if (cfld_taup == NULL)
    BFT_FREE(cpro_taup);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the flow variables relevant to the cooling tower scalars
 * inside the packing zones
 *
 * \param[in,out] liq_mass_flow Liquid mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_init_flow_vars(cs_real_t  liq_mass_flow[])
{

  const cs_real_t  *rho_h = (cs_real_t *)CS_F_(rho)->val; /* humid air (bulk) density */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val; /* mass fraction of liquidus */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /*liquid enthalpy */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /*liquid temperature */

  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->i_face_normal;
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;

  const cs_halo_t *halo = cs_glob_mesh->halo;
  cs_real_t norm_g;
  cs_lnum_t *packing_cell;

  /* Normalised gravity vector */

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  norm_g = cs_math_3_norm(gravity);

  gravity[0] /= norm_g;
  gravity[1] /= norm_g;
  gravity[2] /= norm_g;

  /* Initialise the liquid mass flux to null */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
    liq_mass_flow[face_id] = 0.0;

  /* Tag and initialise the ct values in the packing zone cells */

  BFT_MALLOC(packing_cell, n_cells_with_ghosts, cs_lnum_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++)
    packing_cell[cell_id] = -1;

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_MALLOC(ct->inlet_faces_ids, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_faces_ids, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_cells_ids, n_i_faces, cs_lnum_t);
    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->cell_ids;
    for (int i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];
      packing_cell[cell_id] = ict;

    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_untyped(halo, CS_HALO_STANDARD, sizeof(int), packing_cell);
  }

  /* Initialise the liquid mass flux at packing zone faces
   * and the ghost cells for the liquid mass and enthalpy
   * Initialise the couples (inlet faces, upwind cells) and
   * (outlet faces, upwind cells) arrays */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id_1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id_2 = i_face_cells[face_id][1];

    if (packing_cell[cell_id_1] != -1 || packing_cell[cell_id_2] != -1) {

      int ct_id = CS_MAX(packing_cell[cell_id_1], packing_cell[cell_id_2]);
      cs_ctwr_zone_t *ct = _ct_zone[ct_id];

      // Vertical (align with gravity) component of the surface vector
      cs_real_t liq_surf = cs_math_3_dot_product(gravity, i_face_normal[face_id]);

      /* Face mass flux of the liquid */
      liq_mass_flow[face_id] = ct->q_l_bc / (ct->surface * ct->y_l_bc) * liq_surf;

      /* Initialise a band of ghost cells on the top side of the
         packing zone in order to impose boundary values
         Take the upwinded value for initialisation */

      if (packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] == -1) {

        /* cell_id_2 is an inlet halo */
        if (liq_mass_flow[face_id] < 0.0) {

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
          y_l[cell_id_2] = ct->y_l_bc;
          t_l[cell_id_2] = ct->t_l_bc;
          h_l[cell_id_2] = cs_ctwr_h_liqwater(ct->t_l_bc);
          /* The transported value is (y_l.h_l) and not (h_l) */
          h_l[cell_id_2] *= y_l[cell_id_2];
        }
        /* face_id is an outlet */
        else {

          /* cell_id_2 is an outlet halo */
          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_2;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;

          ct->surface_out += liq_surf;
        }
      }
      else if (packing_cell[cell_id_1] == -1 && packing_cell[cell_id_2] >= 0) {

        /* cell_id_1 is an inlet halo */
        if (liq_mass_flow[face_id] > 0.0) {

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
          y_l[cell_id_1] = ct->y_l_bc;
          t_l[cell_id_1] = ct->t_l_bc;
          h_l[cell_id_1] = cs_ctwr_h_liqwater(ct->t_l_bc);
          /* The transported value is (y_l.h_l) and not (h_l) */
          h_l[cell_id_1] *= y_l[cell_id_1];
        }
        /* cell_id_1 is an outlet */
        else {
          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_1;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;

          ct->surface_out += liq_surf;
        }

        /* Neighbouring zones, inlet for one, outlet fot the other */
      } else if (  packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] >= 0
                && packing_cell[cell_id_1] != packing_cell[cell_id_2]) {

        /* cell_id_1 is an inlet for CT2, an outlet for CT1 */
        if (liq_mass_flow[face_id] > 0.0) {
          /* CT2 */
          ct = _ct_zone[packing_cell[cell_id_2]];

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;

          /* CT1 */
          ct = _ct_zone[packing_cell[cell_id_1]];

          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_1;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;
          ct->surface_out += liq_surf;

        }
        /* cell_id_2 is an inlet for CT1, an outlet for CT2 */
        else {
          /* CT2 */
          ct = _ct_zone[packing_cell[cell_id_2]];

          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_2;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;
          ct->surface_out += liq_surf;

          /* CT1 */
          ct = _ct_zone[packing_cell[cell_id_1]];

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
        }

      }
    } else {
      liq_mass_flow[face_id] = 0.0;
    }
  }

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_REALLOC(ct->inlet_faces_ids, ct->n_inlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_faces_ids, ct->n_outlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_cells_ids, ct->n_outlet_cells, cs_lnum_t);

    cs_parall_sum(1, CS_DOUBLE, &(ct->surface_in));
    cs_parall_sum(1, CS_DOUBLE, &(ct->surface_out));
  }

  BFT_FREE(packing_cell);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset the field variables based on the restart values
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     humidity0   Reference humidity
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_restart_field_vars(cs_real_t  rho0,
                           cs_real_t  t0,
                           cs_real_t  p0,
                           cs_real_t  humidity0,
                           cs_real_t  molmassrat)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_halo_t *halo = m->halo;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;

  // Initialise the fields - based on map
  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;   /* humid air (bulk) density */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;     /* humid air (bulk) Cp */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* humid air temperature */
  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre; /* humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* humid air enthalpy */
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val;    /* water mass fraction in
                                                        humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;     /* humidity in humid air (bulk) */

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;     /* liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;     /* liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;/*liquid mass per unit cell volume*/

  /* liquid vertical velocity component */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;

  cs_field_t *cfld_yp = cs_field_by_name_try("y_p");
  cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_p");
  cs_field_t *cfld_drift_vel = cs_field_by_name_try("drift_vel_y_p");

  cs_real_t *cpro_taup = NULL;
  if (cfld_taup != NULL)
    cpro_taup = cfld_taup->val;
  else
    BFT_MALLOC(cpro_taup, n_cells_with_ghosts, cs_real_t);

  /* Check if there are any leaking packing zones, if yes, there is rain */
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
  for (int ict = 0; ict < _n_ct_zones && !(ct_opt->has_rain); ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];
    if (ct->xleak_fac > 0.0)
      ct_opt->has_rain = true;
  }

  const cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;
  cs_real_t rho_l = ct_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t droplet_diam = ct_prop->droplet_diam;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};


  /* Recompute the initial values which were used in the initialisation of
   * the calculation which is being restarted */
  cs_real_t y_w_ini = humidity0 / ( 1.0 + humidity0); // From 'ctiniv'
  if (y_w_ini < 0.0)
    y_w_ini = 0;

  if (y_w_ini >= 1.0)
    y_w_ini = 1. - cs_math_epzero;

  cs_real_t x_ini = y_w_ini/(1.0-y_w_ini);

  cs_real_t t_h_ini = t0 - cs_physical_constants_celsius_to_kelvin;

  cs_real_t rho_h_ini = cs_ctwr_rho_humidair(x_ini,
                                             rho0,
                                             p0,
                                             t0,
                                             molmassrat,
                                             t_h_ini);

  /* Initialise the cooling towers variables */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field */

    /* Clippings of water mass fraction */
    if (y_w[cell_id] < 0.0)
      y_w[cell_id] = 0; //TODO count it

    if (y_w[cell_id] >= 1.0)
      y_w[cell_id] = 1. - cs_math_epzero; //TODO count it

    x[cell_id] = y_w[cell_id]/(1.0-y_w[cell_id]);

    /* Bulk humid air temperature at the reference temperature
       This is only calculated once at the beginning so same as 'cs_ctwr_init_field_vars'
       No, this would be the value at the previous time step - At present, it is
       not stored in the restart file, so for lack of information initialise it with
       the present value of the temperatur0 */
    t_h_a[cell_id] = t_h[cell_id];

    /* Update the humid air enthalpy based on the solved value of T_h */
    //FIXME Need to use the method of 'cs_ctwr_phyvar_update'

    x_s[cell_id] = cs_ctwr_xsath(t_h[cell_id],p0);

    cp_h[cell_id] = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_ctwr_h_humidair(cp_h[cell_id],
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    /* Update the liquidus temperature based on the solved liquidus enthalpy
     * NB: May not be required as it is also done in 'cs_ctwr_phyvar_update'?
     * No, it must be done here because here we sweep over the entire computational
     * domain whereas 'cs_ctwr_phyvar_update' updates T_l only over the packing zones */
    t_l[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
    if (y_l[cell_id] > 0.) {
      cs_real_t h_liq = h_l[cell_id] / y_l[cell_id];
      t_l[cell_id] = cs_ctwr_t_liqwater(h_liq);
    }

    /* Initialise the liquid vertical velocity component
     * this is correct for droplet and extended for other packing zones
     * NB: this value is derived from the drag coef:
     * C_D = 24 / Re * (1 + 0.15 * Re^0.687)
     * See ZOPLU HT-31-08-06 */

    cs_real_t v_lim = pow(droplet_diam, 2.) * rho_l / (18. * visc)
                    * cs_math_3_norm(gravity);

    cs_real_t reynolds_old = 0.;

    /* Use the same humid air density which was used at the beginning of the
     * calculation being restarted, otherwise since rho_h changes during the
     * calculation, reynolds, v_lim and cpro_taup will end up being different
     * from the initial values used in the calculation being restarted */

    //    cs_real_t reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
    cs_real_t reynolds = rho_h_ini * v_lim * droplet_diam / visc;

    for (int sweep = 0; sweep < 100 && CS_ABS(reynolds - reynolds_old) > 0.001; sweep++) {
      reynolds_old = reynolds;
      v_lim = pow(droplet_diam, 2.) * rho_l / (18. * visc * (1 + 0.15 * pow(reynolds, 0.687)))
            * cs_math_3_norm(gravity);
      //      reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
      reynolds = rho_h_ini * v_lim * droplet_diam / visc;
    }

    cpro_taup[cell_id] = v_lim / cs_math_3_norm(gravity);

    /* Initialize rain variable */
    if (ct_opt->has_rain) {//FIXME useless
      cs_real_3_t *drift_vel = (cs_real_3_t *restrict)(cfld_drift_vel->val);
      drift_vel[cell_id][0] = cpro_taup[cell_id] * gravity[0];
      drift_vel[cell_id][1] = cpro_taup[cell_id] * gravity[1];
      drift_vel[cell_id][2] = cpro_taup[cell_id] * gravity[2];
    }

  }

  /* Loop over exchange zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->cell_ids;
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Update the liquidus temperature based on the solved liquidus enthalpy
       * NB: May not be required as it is also done in 'cs_ctwr_phyvar_update'?
       * No, it must be done here because here we sweep over the entire computational
       * domain whereas 'cs_ctwr_phyvar_update' updates T_l only over the packing zones */

      /* Initialise the liquid vertical velocity component
       * this is correct for dorplet and extended for other packing zones */
      vel_l[cell_id] = cpro_taup[cell_id] * cs_math_3_norm(gravity);

      /* Note that rho_h * Y_l * vel_l * Stot = q_l_bc */
      /* Use the same humid air density which was used at the beginning of the
       * calculation being restarted, otherwise since rho_h changes during the
       * calculation, reynolds, v_lim and cpro_taup will end up being different
       * from the initial values used in the calculation being restarted */
      //      ct->y_l_bc = ct->q_l_bc / (rho_h[cell_id] * vel_l[cell_id] * ct->surface);
      ct->y_l_bc = ct->q_l_bc / (rho_h_ini * vel_l[cell_id] * ct->surface);

      /* Initialise the liquid transported variables:
         liquid mass and enthalpy corrected by the density ratio */
      //FIXME y_l[cell_id] = ct->y_l_bc;

      /* The transported value is (y_l.h_l) and not (h_l) */
      /* No need for that, the value read from the restart file is already
       * the value multiplied by the mass fraction of liquid */
      //FIXME h_l[cell_id] *= y_l[cell_id];

    }
  }

  //Check enthalpies
  cs_real_t h_max = -1.e12;
  cs_real_t h_min = 1.e12;
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    h_min = CS_MIN(h_min,h_h[cell_id]);
    h_max = CS_MAX(h_max,h_h[cell_id]);
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, vel_l);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_taup);
    if (cfld_yp != NULL)
      cs_halo_sync_var(halo, CS_HALO_STANDARD, cfld_yp->val);
    if (cfld_drift_vel != NULL) {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, cfld_drift_vel->val, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, CS_HALO_STANDARD, cfld_drift_vel->val);
    }
  }

  /* Free memory */
  if (cfld_taup == NULL)
    BFT_FREE(cpro_taup);
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the thermo physical properties fields for the humid air and
 *        the liquid
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_phyvar_update(cs_real_t  rho0,
                      cs_real_t  t0,
                      cs_real_t  p0,
                      cs_real_t  molmassrat)
{
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);
  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)(cs_glob_mesh->b_face_cells);
  const cs_halo_t *halo = cs_glob_mesh->halo;

  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;      /* humid air (bulk) density */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;        /* humid air (bulk) Cp */

  // Fields based on maps
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* humid air temperature */
  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre;  /* humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* humid air enthalpy */
  cs_real_t *therm_diff_h = cs_field_by_name_try("thermal_conductivity")->val;
  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;
  cs_real_t *bpro_x1 = cs_field_by_name("b_x_c")->val;
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val; /* Water mass fraction
                                                     in humid air */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val; /* humidity in humid air (bulk) */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /*liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /*liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;      /*liquid mass per unit cell volume*/

  cs_real_t *liq_mass_flow = cs_field_by_name("inner_mass_flux_y_l_packing")->val;//FIXME

  /* Variable and properties for rain zones */
  cs_field_t *cfld_yp = cs_field_by_name_try("y_p");

  cs_real_t *y_p = NULL;
  if (cfld_yp != NULL)
    y_p = cfld_yp->val;

  cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t lambda_h = cs_glob_ctwr_props->lambda_h;
  cs_real_t cp_l = cs_glob_ctwr_props->cp_l;
  cs_real_t lambda_l = cs_glob_ctwr_props->lambda_l;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Clippings of water mass fraction */
    if (y_w[cell_id] < 0.0)
      y_w[cell_id] = 0; //TODO count it

    if (y_w[cell_id] >= 1.0)
      y_w[cell_id] = 1. - cs_math_epzero; //TODO count it

    if (y_p != NULL) {
      if (y_p[cell_id] < 0.0)
        y_p[cell_id] = 0; //TODO count it

      if ((y_p[cell_id] + y_w[cell_id]) >= 1.0)
        y_p[cell_id] = 1. - y_w[cell_id] - cs_math_epzero; //TODO count it

      /* Continuous phase mass fraction */
      cpro_x1[cell_id] = 1. - y_p[cell_id];//TODO not one for rain zones - Why not?
      //If it represents the humid air, then it should be one?  If it represents
      //the dry air, then it should account for both y_p and y_w
    }

    /* Update humidity field */
    x[cell_id] = y_w[cell_id]/(1.0-y_w[cell_id]); //FIXME for drops - This should be the
    //the proportion of 'gaseous' water (dissolved and condensate) in the humid air:
    // Y(dry air)+ Y(gasesous water) + Y(drops) = 1 in all computational cells

    /* Saturated humidity */
    x_s[cell_id] = cs_ctwr_xsath(t_h[cell_id], p0);

    /* Update the humid air temperature using new enthalpy but old
     * Specific heat */

    cp_h[cell_id] = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);

    //FIXME - What is the formula below - Inconsistent with taking into
    //account the saturated phase in the enthalpy in 'cs_ctwr_h_humidair'
    h_h[cell_id] += (t_h[cell_id] - t_h_a[cell_id]) * cp_h[cell_id];

    // Udate the humid air enthalpy diffusivity lambda_h if solve for T_h?
    // Need to update since a_0 is variable as a function of T and humidity
    therm_diff_h[cell_id] = lambda_h;

    /* Update the humid air density */
    rho_h[cell_id] = cs_ctwr_rho_humidair(x[cell_id],
                                          rho0,
                                          p0,
                                          t0,
                                          molmassrat,
                                          t_h[cell_id]);

  }

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->cell_ids;

    /* Packing zone */
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Update the injected liquid temperature
       * NB: (y_l.h_l) is transported and not (h_l) */
      if (y_l[cell_id] > 0.) {
        cs_real_t h_liq = h_l[cell_id] / y_l[cell_id];
        t_l[cell_id] = cs_ctwr_t_liqwater(h_liq);
      }
    }

    /* Update Inlet packing zone temperature if imposed */
    if (ct->delta_t > 0) {
      /* Recompute outgoing temperature */
      ct->t_l_out = 0.0;

      /* Compute liquid water quantities
       * And humid air quantities at liquid outlet */
      for (cs_lnum_t i = 0; i < ct->n_outlet_faces; i++) {

        cs_lnum_t face_id = ct->outlet_faces_ids[i];
        cs_lnum_t cell_id_l, cell_id_h;

        /* Convention: outlet is positive mass flux
         * Then upwind cell for liquid is i_face_cells[][0] */
        int sign = 1;
        if (liq_mass_flow[face_id] < 0) {
          sign = -1;
          cell_id_l = i_face_cells[face_id][1];
          cell_id_h = i_face_cells[face_id][0];
        } else {
          cell_id_l = i_face_cells[face_id][0];
          cell_id_h = i_face_cells[face_id][1];
        }

        /* h_l is in fact (y_l. h_l),
         * and the transport field is (y_l*liq_mass_flow) */
        ct->t_l_out += sign * t_l[cell_id_l]
          * y_l[cell_id_l] * liq_mass_flow[face_id];
        ct->q_l_out += sign * y_l[cell_id_l] * liq_mass_flow[face_id];
      }

      cs_parall_sum(1, CS_DOUBLE, &(ct->t_l_out));
      cs_parall_sum(1, CS_DOUBLE, &(ct->q_l_out));

      ct->t_l_out /= ct->q_l_out;

      /* Relaxation of ct->t_l_bc */
      ct->t_l_bc = (1. - ct->relax) * ct->t_l_bc
                 + ct->relax * (ct->t_l_out + ct->delta_t);

      /* Clippling between 0 and 100 */
      ct->t_l_bc = CS_MAX(CS_MIN(ct->t_l_bc, 100.), 0.);

    }

  }


  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, x);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, x_s);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_x1);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cp_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, h_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, rho_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, t_l);
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    bpro_x1[face_id] = cpro_x1[b_face_cells[face_id]];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in]     p0            Reference pressure
 * \param[in]     molmassrat    dry air to water vapor molecular mass ratio
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_source_term(int              f_id,
                    const cs_real_t  p0,
                    const cs_real_t  molmassrat,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[])
{

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)(m->i_face_cells);
  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)(m->b_face_cells);

  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->b_face_normal;

  cs_real_t  *rho_h = (cs_real_t *)CS_F_(rho)->val; /* humid air (bulk) density */
  cs_real_3_t *vel_h = (cs_real_3_t *)CS_F_(u)->val;   /* humid air (bulk) */

  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val; /* Water mass fraction
                                                     in humid air */

  cs_real_t *t_h = cs_field_by_name("temperature")->val; /* humid air temperature */
  cs_real_t *t_l = cs_field_by_name("temperature_liquid")->val;      /*liquid temperature */
  cs_real_t *x = cs_field_by_name("humidity")->val; /* humidity in humid air (bulk) */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;  /*liquid vertical velocity component */
  cs_real_t *y_l = CS_F_(y_l_pack)->val;

  /* Variable and properties for rain */
  cs_field_t *cfld_yp = cs_field_by_name_try("y_p");         /* Rain drops mass fraction */
  cs_field_t *cfld_tp = cs_field_by_name_try("y_p_t_l");     /* Rain drops temperature */
  cs_field_t *cfld_drift_vel = cs_field_by_name_try("drift_vel_y_p");
  cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_p");

  cs_real_t vertical[3], horizontal[3], norme_g;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

  int evap_model = ct_opt->evap_model;

  /* Need to cook up the cell value of the liquid mass flux
     In the old code, it seems to be taken as the value of the
     face mass flux upstream of the cell */

  cs_ctwr_fluid_props_t *ct_prop = cs_glob_ctwr_props;

  cs_real_t v_air;

  cs_real_t mass_flux_h = 0.; // Highly suspicious for rain zones - not recomputed

  /* Identify the source term formulation for the required field */

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_real_t *f_var = f->val;  /* field variable */

  /* Compute the source terms */

  vertical[0] = -gravity[0];
  vertical[1] = -gravity[1];
  vertical[2] = -gravity[2];

  norme_g = cs_math_3_norm(vertical);

  vertical[0] /= norme_g;
  vertical[1] /= norme_g;
  vertical[2] /= norme_g;
  horizontal[0] = vertical[0] -1.;
  horizontal[1] = vertical[1] -1.;
  horizontal[2] = vertical[2] -1.;

  cs_real_t cp_v = ct_prop->cp_v;
  cs_real_t cp_l = ct_prop->cp_l;
  cs_real_t hv0 = ct_prop->hv0;
  cs_real_t rho_l = ct_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t lambda_h = ct_prop->lambda_h;
  cs_real_t droplet_diam  = ct_prop->droplet_diam;

  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    /* Packing zone characteristics */
    cs_real_t a_0 = ct->xap;
    cs_real_t xnp = ct->xnp;
    int zone_type = ct->type;

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->cell_ids;


    /* Phase change source terms
       ========================= */
    if (evap_model != CS_CTWR_NONE) {

      for (cs_lnum_t j = 0; j < ct->n_cells; j++) {

        cs_lnum_t cell_id = ze_cell_ids[j];

        /* For correlations, T_h cannot be greter than T_l */
        cs_real_t temp_h = CS_MIN(t_h[cell_id], t_l[cell_id]);

        /* saturation humidity at humid air temperature */
        cs_real_t x_s_th = cs_ctwr_xsath(temp_h, p0);

        /* saturation humidity at injected liquid temperature */
        cs_real_t x_s_tl = cs_ctwr_xsath(t_l[cell_id], p0);

        cs_real_t beta_x_a = 0.;
        cs_real_t xlew = 1.;

        /*--------------------------------------------*
         * Counter or cross flow packing zone         *
         *--------------------------------------------*/

        if (   zone_type == CS_CTWR_COUNTER_CURRENT
            || zone_type == CS_CTWR_CROSS_CURRENT) {

          if (zone_type == CS_CTWR_COUNTER_CURRENT) {
            /* Counter flow packing */
            v_air = CS_ABS(cs_math_3_dot_product(vel_h[cell_id], vertical));
          }
          else {
            /* Cross flow packing */
            v_air = CS_ABS(cs_math_3_dot_product(vel_h[cell_id], horizontal));
          }

          /* Dry air flux */
          mass_flux_h = rho_h[cell_id] * v_air * (1. - y_w[cell_id]);

          /* Liquid mass flux */
          cs_real_t mass_flux_l = rho_h[cell_id] * y_l[cell_id] * vel_l[cell_id];

          /* Evaporation coefficient 'Beta_x' times exchange surface 'a' */
          beta_x_a = a_0*mass_flux_l*pow((mass_flux_h/mass_flux_l), xnp);

        }

        //TODO make a private function
        /* Poppe Model */
        if (evap_model == CS_CTWR_POPPE) {
          /* Compute evaporation source terms using Bosnjakovic hypothesis
           * NB: clippings ensuring xi > 1 and xlew > 0 */
          cs_real_t xi = (molmassrat + x_s_tl)/(molmassrat + CS_MIN(x[cell_id], x_s_tl));
          if ((xi - 1.) < 1.e-15)
            xlew = pow(0.866,(2./3.));
          else
            xlew = pow(0.866,(2./3.))*(xi-1.)/log(xi);

        }
        /* Merkel Model */
        else if (evap_model == CS_CTWR_MERKEL) {
          /* Hypothes of Lewis */
          xlew = 1.;
        }

        /* Source terms for the different equations */

        /* Humid air mass source term */
        cs_real_t mass_source = 0.0;
        if (x[cell_id] <= x_s_th) {
          mass_source = beta_x_a*(x_s_tl - x[cell_id]);
        } else {
          mass_source = beta_x_a*(x_s_tl - x_s_th);
        }
        mass_source = CS_MAX(mass_source, 0.);

        cs_real_t vol_mass_source = mass_source * cell_f_vol[cell_id];
        cs_real_t vol_beta_x_a = beta_x_a * cell_f_vol[cell_id];

        /* Global mass source term for continuity (pressure) equation
         * Note that rain is already considered in the bulk, so inner
         * mass transfer between liquid and vapor disappears */
        if (f_id == (CS_F_(p)->id)) {
          /* Warning: not multiplied by Cell volume! no addition neither */
          exp_st[cell_id] = mass_source;
        }

        /* Air mass fraction equation except rain */
        else if (f_id == (CS_F_(ym_w)->id)) {
          exp_st[cell_id] += vol_mass_source*(1. - f_var[cell_id]); //TODO add mass_from_rain
          imp_st[cell_id] += vol_mass_source;
        }

        /* Injected liquid mass equation (solve in drift model form) */
        else if (f_id == (CS_F_(y_l_pack)->id)) {
          exp_st[cell_id] -= vol_mass_source * y_l[cell_id];
          imp_st[cell_id] += vol_mass_source;
        }

        /* Humid air temperature equation */
        else if (f_id == (CS_F_(t)->id)) {
          /* Because the writing is in a non-conservtiv form */
          cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);
          cs_real_t l_imp_st = vol_mass_source * cp_h;
          if (x[cell_id] <= x_s_th) {
            /* Implicit term */
            l_imp_st += vol_beta_x_a * ( xlew * cp_h
                                   + (x_s_tl - x[cell_id]) * cp_v
                                   / (1. + x[cell_id]));
            exp_st[cell_id] += l_imp_st * (t_l[cell_id] - f_var[cell_id]);
          } else {
            cs_real_t coeft = xlew * cp_h;
            /* Implicit term */
            l_imp_st += vol_beta_x_a * ( coeft
                                   + (x_s_tl - x_s_th) * cp_l / (1. + x[cell_id]));
            exp_st[cell_id] += vol_beta_x_a * ( coeft * t_l[cell_id]
                                    + (x_s_tl - x_s_th) * (cp_v * t_l[cell_id] + hv0)
                                    / (1. + x[cell_id])
                                    )
                       - l_imp_st * f_var[cell_id];
          }
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
        }

        /* Injected liquid enthalpy equation (solve in drift model form)
         * NB: it is in fact "y_l x h_l" */
        else if (f_id == (CS_F_(h_l)->id)) {
          /* Implicit term */
          cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);
          cs_real_t l_imp_st = vol_mass_source;
          if (x[cell_id] <= x_s_th) {
            cs_real_t coefh = vol_beta_x_a * ( xlew * cp_h
                                        + (x_s_tl - x[cell_id]) * cp_v
                                        / (1. + x[cell_id]));
            exp_st[cell_id] += coefh * (t_h[cell_id] - t_l[cell_id]);
          } else {
            cs_real_t coefh = xlew * cp_h;
            exp_st[cell_id] += vol_beta_x_a * ( coefh * t_h[cell_id]
                                    - coefh * t_l[cell_id]
                                    + (x_s_tl - x_s_th) / (1. + x[cell_id])
                                      * (  cp_l * t_h[cell_id]
                                        - (cp_v * t_l[cell_id] + hv0)
                                        )
                                    );
          }
          /* Because we deal with an increment */
          exp_st[cell_id] -= l_imp_st * f_var[cell_id];
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);

        }

        /* Injected liquid mass equation for rain zones (solve in drift model form) */
        else if (cfld_yp != NULL) {
          if (f_id == cfld_yp->id) {
            /* Because we deal with an increment */
            exp_st[cell_id] -= vol_mass_source * f_var[cell_id];
            imp_st[cell_id] += vol_mass_source;
            //FIXME other terms ???
          }
        }

      } /* end loop over cells */
    } /* end evaporation model */


    /* Leaking packing zones source terms
       ================================== */
    if (ct->xleak_fac > 0.0) {

      cs_real_t *liq_mass_frac = CS_F_(y_l_pack)->val;   /* liquid mass fraction */
      cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;   /* liquid enthalpy x liquid mass fraction */
      cs_real_t *leak_mass_flow
        = cs_field_by_name("inner_mass_flux_y_l_packing")->val; /* Inner mass flux of liquidus (in the packing) */
      cs_real_t *y_rain = (cs_real_t *)cfld_yp->val;
      cs_real_t *t_rain = (cs_real_t *)cfld_tp->val;

      for (cs_lnum_t i = 0; i < ct->n_outlet_cells; i++) {

        cs_lnum_t cell_id = ct->outlet_cells_ids[i];
        cs_lnum_t face_id = ct->outlet_faces_ids[i];

        cs_real_t inj_mass_flow = ct->xleak_fac * leak_mass_flow[face_id];

        cs_real_t q_inj = liq_mass_frac[cell_id] * inj_mass_flow;
        cs_real_t t_inj = ( cs_ctwr_t_liqwater(h_l[cell_id]/liq_mass_frac[cell_id]) + //FIXME - Careful about T units
            cs_physical_constants_celsius_to_kelvin ) * q_inj;

        /* Global bulk mass - continuity */  //FIXME - Ignore because drops are not in the bulk?
        //        if (f_id == (CS_F_(p)->id)) {
        //          /* Warning: not multiplied by Cell volume! */
        //          exp_st[cell_id] = q_inj;
        //        }

        if (f_id == (CS_F_(p)->id)) {
          /* Global bulk mass - continuity */
          /* Warning: not multiplied by Cell volume! no addition either (why?)*/
          exp_st[cell_id] = q_inj;
        }
        else if (f_id == (CS_F_(ym_w)->id)) {
          //Opposite of rain mass source - make room for drops //FIXME - Check bounds: there should be water in the first place
          cs_real_t vol_mass_source = -q_inj * cell_f_vol[cell_id];

          exp_st[cell_id] += vol_mass_source*(1. - f_var[cell_id]);
          imp_st[cell_id] += vol_mass_source;
        }
        else if (f_id == cfld_yp->id) {
          /* Rain drops mass (particle phase 'p') */
          cs_real_t vol_mass_source = q_inj * cell_f_vol[cell_id];
          //          exp_st[cell_id] -= vol_mass_source * y_rain[cell_id];
          //          imp_st[cell_id] += vol_mass_source;
          exp_st[cell_id] += vol_mass_source;
          imp_st[cell_id] += 0.0;
        }
        else if (f_id == cfld_tp->id) {
          /* Rain drops temperature (particle phase 'p') */
          cs_real_t vol_temp_source = t_inj * cell_f_vol[cell_id];
          //FIXME - Careful about T units
          //FIXME - Cp?  T form not h

          //          exp_st[cell_id] -= vol_temp_source * t_rain[cell_id];
          //          imp_st[cell_id] += vol_temp_source;

          /* Skip for debugging
             exp_st[cell_id] += vol_temp_source;
             imp_st[cell_id] += 0.0;
             */
        }

      }

    } /* End leaking packing zone */

  } /* end packing zone */


  /*--------------------------------------------*/
  /* Rain treatment                             */
  /*--------------------------------------------*/

  if (ct_opt->has_rain) {
    cs_real_3_t *drift_vel = (cs_real_3_t *restrict)(cfld_drift_vel->val);

    for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

      /* For correlations, T_h cannot be greter than T_l */
      cs_real_t temp_h = CS_MIN(t_h[cell_id], t_l[cell_id]);//FIXME is T_l the temperatur of the rain?

      /* saturation humidity at humid air temperature */
      cs_real_t x_s_th = cs_ctwr_xsath(temp_h, p0);

      /* saturation humidity at injected liquid temperature */
      cs_real_t x_s_tl = cs_ctwr_xsath(t_l[cell_id], p0);//FIXME is T_p the temperatur of the rain?

      cs_real_t beta_x_a = 0.;
      cs_real_t xlew = 1.;


      /* Drift velocity of droplets relative to humid air (bulk) */
      cs_real_t dv_p = cs_math_3_norm(drift_vel[cell_id]);

      /* Liquid mass flux divided by v_p */
      cs_real_t mass_flux_l_vp = rho_h[cell_id] * y_l[cell_id];

      cs_real_t cp_h;

      /* Poppe Model */
      if (evap_model == CS_CTWR_POPPE)
        cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s_th);
      /* Merkel Model: humidity is at saturation */
      else
        cp_h = cs_ctwr_cp_humidair(x_s_th, x_s_th);

      /* Reynolds number p. 32 */
      cs_real_t rre = dv_p * rho_h[cell_id]*(1. + x_s_th)*droplet_diam / visc; //FIXME I don't understand 1 + x_s_th

      /* Prandtl number p. 31 */
      cs_real_t rpr = cp_h * visc / lambda_h; //TODO check was cpx

      /* Nusselt number p. 31 */
      cs_real_t anu = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));

      beta_x_a = (6. * lambda_h * anu * mass_flux_l_vp) //FIXME TODO NAT check
        / (0.92 * rho_l * pow(droplet_diam, 2.) * cp_h);//FIXME cpx of the bulk?

      //TODO make a private function
      /* Poppe Model */
      if (evap_model == CS_CTWR_POPPE) {
        /* Compute evaporation source terms using Bosnjakovic hypothesis
         * NB: clippings ensuring xi > 1 and xlew > 0 */
        cs_real_t xi = (molmassrat + x_s_tl)/(molmassrat + CS_MIN(x[cell_id], x_s_tl));
        if ((xi - 1.) < 1.e-15)
          xlew = pow(0.866,(2./3.));
        else
          xlew = pow(0.866,(2./3.))*(xi-1.)/log(xi);

      }
      /* Merkel Model */
      else if (evap_model == CS_CTWR_MERKEL) {
        /* Hypothes of Lewis */
        xlew = 1.;
      }

      /* Source terms for the different equations */

      /* Humid air mass source term */
      cs_real_t mass_source = 0.0;
      if (x[cell_id] <= x_s_th) {
        mass_source = beta_x_a*(x_s_tl - x[cell_id]);
      } else {
        mass_source = beta_x_a*(x_s_tl - x_s_th);
      }
      mass_source = CS_MAX(mass_source, 0.);

      cs_real_t vol_mass_source = mass_source * cell_f_vol[cell_id];
      cs_real_t vol_beta_x_a = beta_x_a * cell_f_vol[cell_id];

      /* Global mass source term for continuity (pressure) equation
       * Note that rain is already considered in the bulk, so inner
       * mass transfer between liquid and vapor disappears */
      if (f_id == (CS_F_(p)->id)) {
        //      /* Warning: not multiplied by Cell volume! no addition neither */
        //      exp_st[cell_id] = mass_source;
      }//FIXME rm because rain is part of the bulk

      /* Air mass fraction equation containing rain */
      else if (f_id == (CS_F_(ym_w)->id)) {
        exp_st[cell_id] += vol_mass_source*(1. - f_var[cell_id]); //TODO NAT Check
        imp_st[cell_id] += vol_mass_source;
      }

      /* Humid air temperature equation */
      else if (f_id == (CS_F_(t)->id)) {
        /* Because the writing is in a non-conservtiv form */
        cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);
        cs_real_t l_imp_st = vol_mass_source * cp_h;
        if (x[cell_id] <= x_s_th) {
          /* Implicit term */
          l_imp_st += vol_beta_x_a * ( xlew * cp_h
              + (x_s_tl - x[cell_id]) * cp_v
              / (1. + x[cell_id]));
          exp_st[cell_id] += l_imp_st * (t_l[cell_id] - f_var[cell_id]);
        } else {
          cs_real_t coeft = xlew * cp_h;
          /* Implicit term */
          l_imp_st += vol_beta_x_a * ( coeft
              + (x_s_tl - x_s_th) * cp_l / (1. + x[cell_id]));
          exp_st[cell_id] += vol_beta_x_a * ( coeft * t_l[cell_id]
              + (x_s_tl - x_s_th) * (cp_v * t_l[cell_id] + hv0)
              / (1. + x[cell_id])
              )
            - l_imp_st * f_var[cell_id];
        }
        imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
      }

      /* Injected liquid mass equation for rain zones (solve in drift model form) */
      if (f_id == cfld_yp->id) {
        /* Because we deal with an increment */
        exp_st[cell_id] -= vol_mass_source * f_var[cell_id];
        imp_st[cell_id] += vol_mass_source;
        //FIXME other terms ???
        //TODO NAT check if it is conservative...
      }

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change mass source term from the evaporating liquid to the
 *        bulk, humid air.
 *
 * Careful, this is different from an injection source term, which would
 * normally be handled with 'cs_user_mass_source_term'
 *
 * \param[in]   p0              Reference pressure
 * \param[in]   molmassrat      Dry air to water vapor molecular mass ratio
 * \param[in]   mass_source     Mass source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_bulk_mass_source_term(const cs_real_t   p0,
                              const cs_real_t   molmassrat,
                              cs_real_t         mass_source[])
{

  cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  /* Compute the mass exchange term */
  cs_real_t *imp_st;

  BFT_MALLOC(imp_st, n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    imp_st[cell_id] = 0.0;
  }

  cs_ctwr_source_term(CS_F_(p)->id, /* Bulk mass source term is
                                       stored for pressure */
      p0,
      molmassrat,
      mass_source,
      imp_st);

  BFT_FREE(imp_st);
}

/*----------------------------------------------------------------------------
 * Get pointer to exchange area.
 *
 * parameters:
 *   ct_id  -->  Id (0 to n-1) of exchange area
 *
 * returns:
 *   pointer to exchange area structure
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t *
cs_ctwr_by_id(int ct_id)
{
  cs_ctwr_zone_t  *retval = NULL;

  if (ct_id > -1 && ct_id <  _n_ct_zones)
    retval = _ct_zone[ct_id];

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

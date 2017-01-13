/*============================================================================
 *
 * Definitions, Global variables variables, and functions associated with the
 * exchange zones
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/* Cooling tower exchange zone structure definition */
/*--------------------------------------------------*/

struct _cs_ctwr_zone_t {

  int                  num;        /* Exchange zone number */
  char                *criteria;   /* Exchange zone elements name */
  cs_ctwr_model_t      model;      /* Exchange model type */
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

  cs_int_t   n_cells;            /* Number of air cells belonging to the zone */

  cs_int_t   up_ct_id;           /* Id of upstream exchange zone (if any) */

  cs_int_t   *ze_cell_list;      /* List of cells of ct criteria */
  cs_int_t   *mark_ze;           /* Cell marker for ct */

  cs_lnum_t  n_inlet_faces;      /* Number of inlet faces */
  cs_lnum_t  n_outlet_faces;     /* Number of outlet faces */
  cs_lnum_t *inlet_faces_list;   /* List of inlet faces */
  cs_lnum_t *outlet_faces_list;  /* List of outlet faces */

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

  cs_real_t  droplet_diam;    /* Drop diameter for rain zones */

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
    cs_real_t *y_l = (cs_real_t *)CS_F_(ym_l)->val;  /* liquid mass per unit
                                                        cell volume*/

    cs_real_t *val;
    BFT_MALLOC(val, mesh->n_cells, cs_real_t);

    /* Value on all cells */

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      val[i] = 0;

    for (int ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
        cs_lnum_t cell_id = ct->ze_cell_list[i];
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cooling tower exchange zone
 *
 * \param[in]   zone_criteria   Zone selection criteria
 * \param[in]   model           model type
 * \param[in]   zone_type       exchange zone type
 * \param[in]   delta_t         Imposed delta temperature delta between inlet
 *                              and oulet of the zone
 * \param[in]   relax           Relaxation of the imposed delta temperature
 * \param[in]   t_l_bc          Liquid water temperature at the inlet
 * \param[in]   q_l_bc          Mass flow rate at the inlet
 * \param[in]   xap             Lambda of the exchange law
 * \param[in]   xnp             Exponent n of the exchange law
 * \param[in]   surface         Total Surface of ingoing water
 * \param[in]   droplet_diam    Droplet diameter in rain zones
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_define(const char           zone_criteria[],
               cs_ctwr_model_t      model,
               cs_ctwr_zone_type_t  zone_type,
               cs_real_t            delta_t,
               cs_real_t            relax,
               cs_real_t            t_l_bc,
               cs_real_t            q_l_bc,
               cs_real_t            xap,
               cs_real_t            xnp,
               cs_real_t            surface,
               cs_real_t            droplet_diam)
{
  cs_ctwr_zone_t  *ct;
  int length;
  FILE *f;
  char  *file_name = NULL;

  /* Definition d'une nouvelle zone d'echange */

  BFT_MALLOC(ct, 1, cs_ctwr_zone_t);

  ct->criteria = NULL;
  BFT_MALLOC(ct->criteria, strlen(zone_criteria)+1, char);
  strcpy(ct->criteria, zone_criteria);

  ct->num = _n_ct_zones + 1;

  ct->model = model;
  ct->type = zone_type;

  ct->delta_t   = delta_t;
  ct->relax = relax;
  ct->t_l_bc   = t_l_bc;
  ct->q_l_bc   = q_l_bc;
  ct->y_l_bc   = 100.; /* Mass of liquid water divided by the mass of humid air
                          in packing zones.
                        Factice version, will be used for rain zones */
  ct->xap = xap;
  ct->xnp = xnp;

  ct->surface_in  = 0.;
  ct->surface_out = 0.;
  ct->surface = surface;

  ct->n_cells = 0;

  ct->up_ct_id = -1;

  ct->n_inlet_faces = 0;
  ct->n_outlet_faces = 0;
  ct->inlet_faces_list = NULL;
  ct->outlet_faces_list = NULL;

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

  ct->droplet_diam = droplet_diam;

  /* Cells selection */
  ct->ze_cell_list = NULL;

  if (_n_ct_zones >= _n_ct_zones_max) {
    _n_ct_zones_max = (_n_ct_zones_max + 1);
    BFT_REALLOC(_ct_zone, _n_ct_zones_max, cs_ctwr_zone_t *);
  }

  /* Add it to exchange zones array */

  _ct_zone[_n_ct_zones] = ct;
  _n_ct_zones += 1;

  if (cs_glob_rank_id <= 0) {
    length = strlen("cooling_towers_balance.") + 3;
    BFT_MALLOC(file_name, length, char);
    sprintf(file_name, "cooling_towers_balance.%02d", ct->num);

    f = fopen(file_name, "a");

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
    BFT_FREE(file_name);
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
  cs_field_pointer_map(CS_ENUMF_(ym_a), cs_field_by_name_try("ym_dry_air"));
  cs_field_pointer_map(CS_ENUMF_(t_l), cs_field_by_name_try("temperature_liquid"));
  cs_field_pointer_map(CS_ENUMF_(h_l), cs_field_by_name_try("enthalpy_liquid"));
  cs_field_pointer_map(CS_ENUMF_(ym_l), cs_field_by_name_try("ym_liquid"));
  cs_field_pointer_map(CS_ENUMF_(thermal_diff_h),
                       cs_field_by_name_try("thermal_conductivity"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cells belonging to the different packing zones.
 *
 * \param[in]   mesh             associated mesh structure
 * \param[in]   mesh_quantities  mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_build_all(const cs_mesh_t              *mesh,
                  const cs_mesh_quantities_t   *mesh_quantities)
{
  /* Create an array for cells flagging */
  /*------------------------------------*/

  cs_ctwr_zone_t  *ct;

  for (int id = 0; id < _n_ct_zones; id++) {

    ct = _ct_zone[id];
    /* Cells selection */
    BFT_MALLOC(ct->ze_cell_list, mesh->n_cells_with_ghosts, cs_lnum_t);

    cs_selector_get_cell_list(ct->criteria, &(ct->n_cells), ct->ze_cell_list);

    BFT_REALLOC(ct->ze_cell_list, ct->n_cells, cs_lnum_t);
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
    BFT_FREE(ct->ze_cell_list);
    BFT_FREE(ct->inlet_faces_list);
    BFT_FREE(ct->outlet_faces_list);
    BFT_FREE(ct);

  }

  _n_ct_zones_max = 0;
  _n_ct_zones = 0;

  BFT_FREE(cs_glob_ctwr_props);

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

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Cooling towers\n"
                  "--------------\n"));

  for (int i = 0; i < _n_ct_zones; i++) {
    cs_ctwr_zone_t *ct = _ct_zone[i];

    char model_type[16];
    if (ct->model == CS_CTWR_NONE) {
      snprintf(model_type, 15, "None");
    } else if (ct->model == CS_CTWR_POPPE) {
      snprintf(model_type, 15, "Poppe");
    } else if (ct->model == CS_CTWR_MERKEL) {
      snprintf(model_type, 15, "Merkel");
    }

    cs_log_printf
      (CS_LOG_SETUP,
       _("  Cooling tower zone id: %d\n"
         "    criterion: ""%s""\n"
         "    Parameters:\n"
         "      Model: %s\n"
         "      Lambda of the exchange law: %f\n"
         "      Exponent n of the exchange law: %f\n"
         "      Type: %d\n"
         "        Droplet diameter: %f\n"
         "      Delta Temperature: %f\n"
         "        Relaxation: %f\n"
         "      Injected water temperature: %f\n"
         "      Injected mass flow rate: %f\n"
         "      Total surface of ingoing water: %f\n"),
       ct->num,
       ct->criteria,
       model_type,
       ct->xap,
       ct->xnp,
       ct->type,
       ct->droplet_diam,
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
  cs_real_t *y_a = (cs_real_t *)CS_F_(ym_a)->val;   /* dry air mass fraction
                                                       in humid air */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;    /* humid air bulk humidity */

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /* liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /* liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(ym_l)->val;   /* liquid mass per unit
                                                       cell volume */

  cs_real_t *liq_mass_flow = cs_field_by_name("inner_mass_flux_ym_liquid")->val;
  cs_real_t *mass_flow = cs_field_by_name("inner_mass_flux")->val;

  int length;
  FILE *f;
  char  *file_name = NULL;
  cs_real_t cp_l = cs_glob_ctwr_props->cp_l;

  /* Loop over Cooling tower zones */
  for (int ict=0; ict < _n_ct_zones; ict++) {

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

      cs_lnum_t face_id = ct->inlet_faces_list[i];
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
     * And humid air quantities at liquid outlet */
    for (cs_lnum_t i = 0; i < ct->n_outlet_faces; i++) {

      cs_lnum_t face_id = ct->outlet_faces_list[i];
      cs_lnum_t cell_id_l, cell_id_h;

      /* Convention: outlet is positiv mass flux
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
      length = strlen("cooling_towers_balance.") + 3;
      BFT_MALLOC(file_name, length, char);
      sprintf(file_name, "cooling_towers_balance.%02d", ct->num);

      if (CS_ABS(ct->h_l_in - ct->h_l_out)> 1.e-6) {
        f = fopen(file_name, "a");
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

    BFT_FREE(file_name);
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
  // Initialise the fields - based on map
  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;   /* humid air (bulk) density */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* humid air temperature */
  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre; /* humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* humid air enthalpy */
  cs_real_t *y_a = (cs_real_t *)CS_F_(ym_a)->val;    /* dry air mass fraction in
                                                        humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val; /* humidity in humid air (bulk) */

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;     /* liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;     /* liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(ym_l)->val;    /* liquid mass per unit
                                                        cell volume */

  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field in case users have updated the initial
       dry air mass fraction.
       Note: this is a bit dubious as users could also have chosen
       to reset the humidity ? */

    if (y_a[cell_id] > 0.0 && y_a[cell_id] <= 1.0) {
      x[cell_id] = (1.0-y_a[cell_id])/y_a[cell_id];
    } else {
      //Put stop signal here - !!FIXME
    }

    /* Bulk humid air temperature */
    t_h[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
    t_h_a[cell_id] = t_h[cell_id];

    // Update the humid air density
    rho_h[cell_id] = cs_ctwr_rho_humidair(x[cell_id],
                                          rho0,
                                          p0,
                                          t0,
                                          molmassrat,
                                          t_h[cell_id]);

    // Update the humid air enthalpy
    x_s[cell_id] = cs_ctwr_xsath(t_h[cell_id],p0);
    cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_ctwr_h_humidair(cp_h,
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

  }

  /* Loop over exchange zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ct->ze_cell_list[i];

      /* Initialize with the injection water temperature */
      t_l[cell_id] = ct->t_l_bc;

      /* Update the injected liquid enthalpy */
      h_l[cell_id] = cs_ctwr_h_liqwater(t_l[cell_id]);

      /* Initialise the liquid transported variables:
         liquid mass and enthalpy corrected by the density ratio */
      y_l[cell_id] = ct->y_l_bc;

      /* The transported value is (y_l.h_l) and not (h_l) */
      h_l[cell_id] *= y_l[cell_id];

    }
  }
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
  /* humid air (bulk) density */
  const cs_real_t  *rho_h = (cs_real_t *)CS_F_(rho)->val;
  cs_real_t *y_l = cs_field_by_name("ym_liquid")->val;
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /*liquid enthalpy */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /*liquid temperature */

  /* liquid vertical velocity component */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;

  const cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->i_face_normal;

  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;

  const cs_halo_t *halo = cs_glob_mesh->halo;

  cs_real_t gravity[3], norm_g;
  cs_lnum_t *packing_cell;

  cs_lnum_t cell_id;

  /* Normalised gravity vector */

  gravity[0] = ct_prop->gravx;
  gravity[1] = ct_prop->gravy;
  gravity[2] = ct_prop->gravz;

  norm_g = cs_math_3_norm(gravity);

  gravity[0] /= norm_g;
  gravity[1] /= norm_g;
  gravity[2] /= norm_g;

  /* Tag and initialise the ct values in the packing zone cells */

  BFT_MALLOC(packing_cell, n_cells_with_ghosts, cs_lnum_t);

  for (cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++)
    packing_cell[cell_id] = -1;

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_MALLOC(ct->inlet_faces_list, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_faces_list, n_i_faces, cs_lnum_t);
    for (int i = 0; i < ct->n_cells; i++) {
      cell_id = ct->ze_cell_list[i];
      packing_cell[cell_id] = ict;
      /* Initialise the liquid vertical velocity component */
      vel_l[cell_id] = ct->q_l_bc / (rho_h[cell_id] * ct->y_l_bc * ct->surface);
    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, vel_l);
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
      cs_real_t y_l_bc = ct->y_l_bc;
      liq_mass_flow[face_id] = ct->q_l_bc / (ct->surface * ct->y_l_bc) * liq_surf;

      /* Initialise a band of ghost cells on the top side of the
         packing zone in order to impose boundary values
         Take the upwinded value for initialisation */

      if (packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] == -1) {

        /* cell_id_2 is an inlet halo */
        if (liq_mass_flow[face_id] < 0.0) {

          ct->inlet_faces_list[ct->n_inlet_faces] = face_id;

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
          ct->outlet_faces_list[ct->n_outlet_faces] = face_id;

          ct->n_outlet_faces ++;
          ct->surface_out += liq_surf;
        }
      }
      else if (packing_cell[cell_id_1] == -1 && packing_cell[cell_id_2] >= 0) {

        /* cell_id_1 is an inlet halo */
        if (liq_mass_flow[face_id] > 0.0) {

          ct->inlet_faces_list[ct->n_inlet_faces] = face_id;

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
          ct->outlet_faces_list[ct->n_outlet_faces] = face_id;

          ct->n_outlet_faces ++;
          ct->surface_out += liq_surf;
        }

        /* Neighbouring zones, inlet for one, outlet fot the other */
      } else if (  packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] >= 0
                && packing_cell[cell_id_1] != packing_cell[cell_id_2]) {

        /* cell_id_1 is an inlet for CT2, an outlet for CT1 */
        if (liq_mass_flow[face_id] > 0.0) {
          /* CT2 */
          ct = _ct_zone[packing_cell[cell_id_2]];

          ct->inlet_faces_list[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;

          /* CT1 */
          ct = _ct_zone[packing_cell[cell_id_1]];

          ct->outlet_faces_list[ct->n_outlet_faces] = face_id;

          ct->n_outlet_faces ++;
          ct->surface_out += liq_surf;

        }
        /* cell_id_2 is an inlet for CT1, an outlet for CT2 */
        else {
          /* CT2 */
          ct = _ct_zone[packing_cell[cell_id_2]];

          ct->outlet_faces_list[ct->n_outlet_faces] = face_id;

          ct->n_outlet_faces ++;
          ct->surface_out += liq_surf;

          /* CT1 */
          ct = _ct_zone[packing_cell[cell_id_1]];

          ct->inlet_faces_list[ct->n_inlet_faces] = face_id;

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

    BFT_REALLOC(ct->inlet_faces_list, ct->n_inlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_faces_list, ct->n_outlet_faces, cs_lnum_t);

    cs_parall_sum(1, CS_DOUBLE, &(ct->surface_in));
    cs_parall_sum(1, CS_DOUBLE, &(ct->surface_out));
  }

  BFT_FREE(packing_cell);
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
  const cs_halo_t *halo = cs_glob_mesh->halo;

  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;      /* humid air (bulk) density */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;        /* humid air (bulk) Cp */

  // Fields based on maps
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* humid air temperature */
  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre;  /* humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* humid air enthalpy */
  cs_real_t *therm_diff_h = cs_field_by_name_try("thermal_conductivity")->val;
  cs_real_t *y_a = (cs_real_t *)CS_F_(ym_a)->val;       /* dry air mass fraction in humid air */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val; /* humidity in humid air (bulk) */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /*liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /*liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(ym_l)->val;      /*liquid mass per unit cell volume*/

  cs_real_t *liq_mass_flow = cs_field_by_name("inner_mass_flux_ym_liquid")->val;

  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_t lambda_h = cs_glob_ctwr_props->cond_h;
  cs_real_t cp_l = cs_glob_ctwr_props->cp_l;
  cs_real_t lambda_l = cs_glob_ctwr_props->cond_l;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field */
    if (y_a[cell_id] > 0.0 && y_a[cell_id] <= 1.0) {
      x[cell_id] = (1.0-y_a[cell_id])/y_a[cell_id];
    } else {
      //Put stop signal here - !!FIXME
    }

    // Update the humid air Cp - Not completely right
    // ultimately, should have a coupled solution on t_h and humidity
    x_s[cell_id] = cs_ctwr_xsath(t_h[cell_id], p0);

    /* Update the humid air temperature using new enthalpy but old
     * Specific heat */

    cp_h[cell_id] = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] += (t_h[cell_id] - t_h_a[cell_id]) * cp_h[cell_id];

    // Udate the umid air enthalpy diffusivity lambda_h if solve for T_h?
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

    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ct->ze_cell_list[i];

      /* Update the injected liquid temperature
       * NB: (y_l.h_l) is transported and not (h_l) */
      if (y_l[cell_id] > 0.) {
        cs_real_t h_liq = h_l[cell_id]/y_l[cell_id];
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

        cs_lnum_t face_id = ct->outlet_faces_list[i];
        cs_lnum_t cell_id_l, cell_id_h;

        /* Convention: outlet is positiv mass flux
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
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cp_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, h_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, rho_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, t_l);
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
  cs_lnum_t  iloc = 0;

  cs_real_t  *rho_h = (cs_real_t *)CS_F_(rho)->val; /* humid air (bulk) density */
  cs_real_3_t *u_air = (cs_real_3_t *)CS_F_(u)->val;   /* humid air (bulk) */

  cs_real_t *y_a = (cs_real_t *)CS_F_(ym_a)->val;       /* dry air mass fraction in humid air */

  cs_real_t *t_h = cs_field_by_name("temperature")->val; /* humid air temperature */
  cs_real_t *h_h   = cs_field_by_name("enthalpy")->val;    /* humid air enthalpy */
  cs_real_t *t_l = cs_field_by_name("temperature_liquid")->val;      /*liquid temperature */
  cs_real_t *x = cs_field_by_name("humidity")->val; /* humidity in humid air (bulk) */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;  /*liquid vertical velocity component */
  cs_real_t *y_l = cs_field_by_name("ym_liquid")->val;

  cs_real_t vertical[3], horizontal[3], norme_g;
  cs_real_t vvai, vhai;
  cs_real_t dvg, cpx, rre, rpr, anu;
  cs_real_t cp_l, cp_v, cp_a, visc, conduc;

  /* Need to cook up the cell value of the liquid mass flux
     In the old code, it seems to be taken as the value of the
     face mass flux upstream of the cell */
  cs_real_t mass_flux_l;      /* injected liquid mass flux */

  cs_ctwr_fluid_props_t *ct_prop = cs_glob_ctwr_props;

  cs_real_t v_air, xi;

  cs_real_t mass_flux_h = 0.; // Highly suspicious for rain zones - not recomputed

  /* Identify the source term formulation for the required field */

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_real_t *f_var = f->val;  /* field variable */

  /* Compute the source terms */

  vertical[0] = -ct_prop->gravx;
  vertical[1] = -ct_prop->gravy;
  vertical[2] = -ct_prop->gravz;

  norme_g = cs_math_3_norm(vertical);

  vertical[0] /= norme_g;
  vertical[1] /= norme_g;
  vertical[2] /= norme_g;
  horizontal[0] = vertical[0] -1.;
  horizontal[1] = vertical[1] -1.;
  horizontal[2] = vertical[2] -1.;

  cp_a = ct_prop->cp_a;
  cp_v = ct_prop->cp_v;
  cp_l = ct_prop->cp_l;
  cs_real_t hv0 = ct_prop->hv0;
  cs_real_t rho_l = ct_prop->rho_l;
  visc = ct_prop->visc;
  conduc = ct_prop->cond_h ;

  cs_lnum_t i = 0;

  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    /* Packing zone characteristics */
    cs_real_t drop_diam  = ct->droplet_diam;
    cs_real_t a_0 = ct->xap;
    cs_real_t xnp = ct->xnp;
    int zone_type = ct->type;
    int evap_model = ct->model;

    if (evap_model > CS_CTWR_NONE) {

      for (cs_lnum_t j = 0; j < ct->n_cells; j++) {

        cs_lnum_t cell_id = ct->ze_cell_list[j];

        /* For correlations, T_h cannot be greter than T_l */
        cs_real_t temp_h = CS_MIN(t_h[cell_id], t_l[cell_id]);

        /* saturation humidity at humid air temperature */
        cs_real_t x_s_th = cs_ctwr_xsath(temp_h, p0);

        /* saturation humidity at injected liquid temperature */
        cs_real_t x_s_tl = cs_ctwr_xsath(t_l[cell_id], p0);

        cs_real_t beta_x_a, xlew;

        if (evap_model == CS_CTWR_POPPE) {

          /*--------------------------------------------*
           * Poppe Model                                *
           *--------------------------------------------*/

          if (   zone_type == CS_CTWR_COUNTER_CURRENT
              || zone_type == CS_CTWR_CROSS_CURRENT) {

            /*--------------------------------------------*
             * Counter or cross flow packing zone         *
             *--------------------------------------------*/

            if (zone_type == CS_CTWR_COUNTER_CURRENT) {
              /* Counter flow packing */
              v_air = CS_ABS(cs_math_3_dot_product(u_air[cell_id], vertical));
            }
            else {
              /* Cross flow packing */
              v_air = CS_ABS(cs_math_3_dot_product(u_air[cell_id], horizontal));
            }

            /* Dry air flux */
            mass_flux_h = rho_h[cell_id] * v_air * y_a[cell_id];

            /* Liquid mass flux */
            mass_flux_l = rho_h[cell_id] * y_l[cell_id] * vel_l[cell_id];

            /* Evaporation coefficient 'Beta_x' times exchange surface 'a' */
            beta_x_a = a_0*mass_flux_l*pow((mass_flux_h/mass_flux_l), xnp);

            /* Compute evaporation source terms using Bosnjakovic hypothesis
             * NB: clippings ensuring xi > 1 and xlew > 0 */
            xi = (molmassrat + x_s_tl)/(molmassrat + CS_MIN(x[cell_id], x_s_tl));
            if ((xi - 1.) < 1.e-15)
              xlew = pow(0.866,(2./3.));
            else
              xlew = pow(0.866,(2./3.))*(xi-1.)/log(xi);

          }

          else if (zone_type == CS_CTWR_RAIN) {//FIXME

            /*--------------------------------------------*/
            /* Rain zone                                  */
            /*--------------------------------------------*/

            cs_real_t *vgin = NULL;
            if (CS_ABS(vgin[iloc]) >= 0.1) { /* vgin looks like the drop velocity */
              /* Is it the modulus ? */
              vvai = CS_ABS(cs_math_3_dot_product(u_air[cell_id], vertical));
              vhai = CS_ABS(cs_math_3_dot_product(u_air[cell_id], horizontal));

              dvg = sqrt(pow((vvai+vgin[iloc]),2.)+pow(vhai,2.)); /* This looks wrong - should be
                                                                     the difference: relative
                                                                     velocity p. 32 */

              //Looks wrong too: mass_flux_h = 0. from initialisation at the top
              //Why is it here anyway since it is recalculated below? - old copy/paste error?
              beta_x_a = a_0*mass_flux_l*pow((mass_flux_h/mass_flux_l),xnp);

              if (x[cell_id] <= x_s_th) {
                cpx = cp_a + x[cell_id]*cp_v;
              }
              else {
                cpx = cp_a + x_s_th*cp_v + (x[cell_id] - x_s_th)*cp_l;
              }

              rre = dvg*rho_h[cell_id]*(1. + x_s_th)*drop_diam/visc; /* Reynolds number p. 32 */
              rpr = cpx*visc/conduc; /* Prandtl number p. 31 */
              anu = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.)); /* Nusselt number p. 31 */

              beta_x_a = (6.*conduc*anu*mass_flux_l)/(0.92*rho_l*vgin[iloc]*pow(drop_diam,2.)*cpx);

              /* Compute evaporation source terms using Bosnjakovic hypothesis
               * NB: clippings ensuring xi > 1 and xlew > 0 */ //FIXME xi not computed
              xi = (molmassrat + x_s_tl)/(molmassrat + CS_MIN(x[cell_id], x_s_tl));
              if ((xi - 1.) < 1.e-15)
                xlew = pow(0.866,(2./3.));
              else
                xlew = pow(0.866,(2./3.))*(xi-1.)/log(xi);

            }
          }
        }
        else if (evap_model == CS_CTWR_MERKEL) {

          /*--------------------------------------------*
           * Merkel Model                               *
           *--------------------------------------------*/

          if (zone_type <= 2) {

            /*--------------------------------------------*
             * Counter or cross flow packing zone         *
             *--------------------------------------------*/

            /* Hypothes of Lewis */
            xlew = 1.;

            /* Liquid mass flux - Fe in reference */
            mass_flux_l = rho_h[cell_id] * y_l[cell_id] * vel_l[cell_id];

            if (mass_flux_l > 1.e-6) {

              /* Counter flow packing */
              if (zone_type == CS_CTWR_COUNTER_CURRENT) {
                v_air = CS_ABS(cs_math_3_dot_product(u_air[cell_id], vertical));
              }
              /* Cross flow packing */
              else {
                v_air = CS_ABS(cs_math_3_dot_product(u_air[cell_id], horizontal));
              }

              /* Dry air flux - Fa in reference */
              mass_flux_h = rho_h[cell_id] * v_air * y_a[cell_id];

              /* Evaporation coefficient Beta_x times exchange surface 's' */
              beta_x_a = a_0*mass_flux_l*pow((mass_flux_h/mass_flux_l),xnp);

            }
          }

          else if (zone_type == CS_CTWR_RAIN) {  //FIXME

            /*--------------------------------------------*/
            /* Rain zone                                  */
            /*--------------------------------------------*/

            cs_real_t *vgin = NULL;
            if (CS_ABS(vgin[iloc])>=0.1) {

              vvai = CS_ABS(cs_math_3_dot_product(u_air[cell_id], vertical));
              vhai = CS_ABS(cs_math_3_dot_product(u_air[cell_id], horizontal));
              dvg = sqrt(pow((vvai+vgin[iloc]),2.)+pow(vhai,2.)); //FIXME "vvai+vgin" should be "vvai-vgin"

              cpx = cp_a + x_s_th*cp_v;
              rre = dvg*rho_h[cell_id]*(1. + x_s_th)*drop_diam/visc;
              rpr = cpx*visc/conduc;
              anu = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));

              beta_x_a = (6.*conduc*anu*mass_flux_l)/(0.92*rho_l*vgin[iloc]*pow(drop_diam,2.)*cpx);

            }
          }
        } /* end evaporation model */

        /* Source terms for the different equations */

        // Humid air mass source term
        cs_real_t mass_source = 0.0;
        if (x[cell_id] <= x_s_th) {
          mass_source = beta_x_a*(x_s_tl - x[cell_id]);
        } else {
          mass_source = beta_x_a*(x_s_tl - x_s_th);
        }
        mass_source = CS_MAX(mass_source, 0.);

        /* Global continuity (pressure) equation */
        if (f_id == (CS_F_(p)->id)) {
          exp_st[i] = mass_source;
          imp_st[i] = 0.0;
        }

        /* Dry air mass fraction equation */
        else if (f_id == (CS_F_(ym_a)->id)) {
          exp_st[i] = -mass_source*f_var[cell_id];
          imp_st[i] = CS_MAX(mass_source, 0.);
        }

        /* Injected liquid mass equation (solve in drift model form) */
        else if (f_id == (CS_F_(ym_l)->id)) {
          exp_st[i] = -mass_source * y_l[cell_id];
          imp_st[i] = CS_MAX(mass_source, 0.);
        }

        /* Humid air temperature equation */
        else if (f_id == (CS_F_(t)->id)) {
          /* Because the writing is in a non-conservtiv form */
          cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);
          imp_st[i] = CS_MAX(mass_source, 0.)*cp_h; //TODO check was working before fix...
          if (x[cell_id] <= x_s_th) {
            /* Implicit term */
            imp_st[i] += beta_x_a * ( xlew * cp_h
                                   + (x_s_tl - x[cell_id]) * cp_v
                                   / (1. + x[cell_id]));
            exp_st[i] += imp_st[i] * (t_l[cell_id] - f_var[cell_id]);
          } else {
            cs_real_t coeft = xlew * cp_h;
            /* Implicit term */
            imp_st[i] += beta_x_a * ( coeft
                                    + (x_s_tl - x_s_th) * cp_l / (1. + x[cell_id]));
            exp_st[i] += beta_x_a * ( coeft * t_l[cell_id]
                                    + (x_s_tl - x_s_th) * (cp_v * t_l[cell_id] + hv0)
                                    / (1. + x[cell_id])
                                    )
                       - imp_st[i] * f_var[cell_id];
          }
          imp_st[i] = CS_MAX(imp_st[i], 0.);
        }

        /* Injected liquid enthalpy equation (solve in drift model form)
         * NB: it is in fact "y_l x h_l" */
        else if (f_id == (CS_F_(h_l)->id)) {
          /* Implicit term */
          cs_real_t cp_h = cs_ctwr_cp_humidair(x[cell_id], x_s[cell_id]);
          imp_st[i] = CS_MAX(mass_source, 0.);
          if (x[cell_id] <= x_s_th) {
            cs_real_t coefh = beta_x_a * ( xlew * cp_h
                                        + (x_s_tl - x[cell_id]) * cp_v
                                        / (1. + x[cell_id]));
            exp_st[i] = coefh * (t_h[cell_id] - t_l[cell_id]);
          } else {
            cs_real_t coefh = xlew * cp_h;
            exp_st[i] += beta_x_a * ( coefh * t_h[cell_id]
                                    - coefh * t_l[cell_id]
                                    + (x_s_tl - x_s_th) / (1. + x[cell_id])
                                      * (  cp_l * t_h[cell_id]
                                        - (cp_v * t_l[cell_id] + hv0)
                                        )
                                    );
          }
          /* Because we deal with an increment */
          exp_st[i] -= imp_st[i] * f_var[cell_id];

        }

        i++;

      } /* end loop over cells */
    } /* end evaporation model */
  } /* end packing zone */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change mass source term from the evaporating liquid to the
 *        bulk, humid air.
 *
 * Careful, this is different from an injection source term, which would
 * normally be handled with 'cs_user_mass_source_term'
 *
 * \param[in]   call_id          Calling sequence flag
 * \param[in]   p0              Reference pressure
 * \param[in]   molmassrat      Dry air to water vapor molecular mass ratio
 * \param[in]   n_tot           Pointer to the total number
 *                              of cells in the packing zones
 * \param[in]   packing_cell    Packing cell ids
 * \param[in]   mass_source     Mass source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_bulk_mass_source_term(int               call_id,
                              const cs_real_t   p0,
                              const cs_real_t   molmassrat,
                              int              *n_tot,
                              cs_lnum_t         packing_cell[],
                              cs_real_t         mass_source[])
{
  if (call_id == 1) {
    // Count the total number of cells in which the source term will be applied
    // This is the total number of cells in the packing regions

    if (*n_tot != 0) {
      //Error message - This would indicate that 'cs_user_mass_source_term' is also
      //in use, which would be inconsistent with activating the cooling towers model

    }
    for (cs_lnum_t ict = 0; ict < _n_ct_zones; ict++) {
      cs_ctwr_zone_t *ct = _ct_zone[ict];
      *n_tot = *n_tot + (ct->n_cells);
    }

  } else if (call_id == 2) {

    // Fill in the array of cells in which the source term will be applied
    // These are the cells located in the packing regions

    cs_lnum_t i = 0;

    for (cs_lnum_t ict = 0; ict < _n_ct_zones; ict++) {

      cs_ctwr_zone_t *ct = _ct_zone[ict];

      for (cs_lnum_t j = 0; j < ct->n_cells; j++) {
        cs_lnum_t cell_id = ct->ze_cell_list[j];
        /* Careful, cell number and not cell id because used in Fortran */
        packing_cell[i] = cell_id + 1;
        i++;
      }

    }

  } else if (call_id == 3) {

    // Compute the mass exchange term
    cs_real_t *exp_st;
    cs_real_t *imp_st;

    BFT_MALLOC(exp_st, *n_tot, cs_real_t);
    BFT_MALLOC(imp_st, *n_tot, cs_real_t);

    for (cs_lnum_t i = 0; i < *n_tot; i++) {
      exp_st[i] = 0.0;
      imp_st[i] = 0.0;
    }


    cs_ctwr_source_term(CS_F_(p)->id, /* Bulk mass source term is
                                         stored for pressure */
                        p0,
                        molmassrat,
                        exp_st,
                        imp_st);

    for (cs_lnum_t i = 0; i < *n_tot; i++) {
      mass_source[i] = mass_source[i] + exp_st[i];
    }

    BFT_FREE(exp_st);
    BFT_FREE(imp_st);
  }
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

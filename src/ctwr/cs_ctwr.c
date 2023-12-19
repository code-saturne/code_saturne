/*============================================================================
 * Cooling towers related functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
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
  .has_rain = false,
  .solve_rain_velocity = false};

const cs_ctwr_option_t *cs_glob_ctwr_option = &_ctwr_option;

/* Cooling tower exchange zone structure definition */
/*--------------------------------------------------*/

struct _cs_ctwr_zone_t {

  int                  num;        /* Exchange zone number */
  char                *criteria;   /* Exchange zone selection criteria */
  int                  z_id;       /* id of the volume zone */
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

  cs_real_t  xap;                /* Exchange law a_0 coefficient */
  cs_real_t  xnp;                /* Exchange law n exponent */

  cs_real_t  surface_in;         /* Water inlet surface */
  cs_real_t  surface_out;        /* Water outlet surface */
  cs_real_t  surface;            /* Total surface */

  cs_real_t  xleak_fac;          /* Leakage factor (ratio of outlet/inlet
                                    flow rate) */
  cs_real_t  v_liq_pack;         /* Vertical liquid film velocity in packing */

  cs_lnum_t  n_cells;            /* Number of air cells belonging to the zone */
  cs_real_t  vol_f;              /* Cooling tower zone total volume */

  int        up_ct_id;           /* Id of upstream exchange zone (if any) */

  cs_lnum_t  n_inlet_faces;      /* Number of inlet faces */
  cs_lnum_t  n_outlet_faces;     /* Number of outlet faces */
  cs_lnum_t *inlet_faces_ids;    /* List of inlet faces */
  cs_lnum_t *outlet_faces_ids;   /* List of outlet faces */

  cs_lnum_t  n_outlet_cells;     /* Number of outlet cells */
  cs_lnum_t *outlet_cells_ids;   /* List of outlet cells */

  cs_real_t  p_in;            /* Average inlet pressure */
  cs_real_t  p_out;           /* Average outlet pressure */
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

    cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;   /* Liquid enthalpy */
    cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;  /* Liquid mass per unit
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

/*----------------------------------------------------------------------------
 * Compute the Lewis factor used for the evaluation of the heat transfer
 * phase change source terms
 *
 * parameters:
 *   evap_model  <-- Evaporation model: CS_CTWR_POPPE or CS_CTWR_MERKEL
 *   molmassrat  <-- Dry air to water vapor molecular mass ratio
 *   x           <-- Humidity
 *   x_s_tl      <-- Saturation humidity at the temperature of the liquid
 *
 * returns:
 *   xlew        --> Lewis factor
 *----------------------------------------------------------------------------*/

static cs_real_t
_lewis_factor(const int        evap_model,
              const cs_real_t  molmassrat,
              const cs_real_t  x,
              const cs_real_t  x_s_tl)
{
  /* Merkel Model
     Hypothesis of unity Lewis factor */
  cs_real_t xlew = 1.;

  if (evap_model == CS_CTWR_POPPE) {
    /* Poppe evaporation model
       Compute Lewis factor using Bosnjakovic hypothesis
       NB: clippings ensuring xi > 1 and xlew > 0 */
    cs_real_t xi = (molmassrat + x_s_tl)/(molmassrat + CS_MIN(x, x_s_tl));
    if ((xi - 1.) < 1.e-15)
      xlew = pow(0.866,(2./3.));
    else
      xlew = pow(0.866,(2./3.))*(xi-1.)/log(xi);
  }

  return xlew;
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
 * Add variables fields
 *----------------------------------------------------------------------------*/

void
cs_ctwr_add_variable_fields(void)
{
  /* Key id of the scalar class */
  const int keyccl = cs_field_key_id("scalar_class");

  /* Key id for drift scalar */
  const int keydri = cs_field_key_id("drift_scalar_model");

  /* Key ids for clipping */
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  /* Key id for the diffusivity */
  const int kivisl = cs_field_key_id("diffusivity_id");

  /* Fluid properties and physical variables */
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

  /* Set fluid properties parameters */

  /* Variable density */
  fp->irovar = 1;
  /* Activate compressibility */
  cs_velocity_pressure_model_t *vp_model =
    cs_get_glob_velocity_pressure_model();
  vp_model->idilat = 2;
  /* Constant molecular viscosity */
  fp->ivivar = 0;

  /* 1. Definition of fields
   * --------------------------------------------------------------------------
   *  Bulk definition - For cooling towers, the bulk is the humid air.
   *  By definition, humid air is composed of two species: dry air and water
   *  vapor (whether in gas or condensate form)
   *  -------------------------------------------------------------------------
   */

  cs_field_t *f;

  {
    /* Thermal model - Set parameters of calculations (module optcal) */

    cs_thermal_model_t *thermal = cs_get_glob_thermal_model();

    /* Solve for temperature of bulk humid air */
    thermal->thermal_variable = CS_THERMAL_MODEL_TEMPERATURE;

    /* Temperature treated in Celsius */
    thermal->itpscl = CS_TEMPERATURE_SCALE_CELSIUS;

    /* Variable cp (0 = variable, -1 = constant) since it changed with humidity
     * Needs to be specified here because the automated creation and
     * initialization of the cell array for cp in 'iniva0' depends on its value
     * (unlike the cell arrays for the density and viscosity which are
     * initialized irrespective of the values of irovar and ivivar) */
    fp->icp = 0;

    /* The thermal transported scalar is the temperature of the bulk.
     * If the atmospheric module is switched off (i.e., iatmos!= 2)
     * , we create the field. */
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != CS_ATMO_HUMID){
      int f_id = cs_variable_field_create("temperature",
                                          "Temperature humid air",
                                          CS_MESH_LOCATION_CELLS,
                                          1);

      f = cs_field_by_id(f_id);

      /* Set variable diffusivity for the humid air enthalpy.
       * The diffusivity used in the transport equation will be the cell value
       * of the viscls array for f_id.
       * This value is updated at the top of each time step in 'ctphyv' along
       * with the other variable properties */
      int ifcvsl = 0;
      cs_field_set_key_int(f, kivisl, ifcvsl);
      cs_add_model_thermal_field_indexes(f->id);
    }
  }

  {
    /* Rain zone variables */

    /* Associate liquid water rain with class 1 */
    int class_id = 1;

    int f_id = cs_variable_field_create("y_p",
                                        "Yp rain",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    f = cs_field_by_id(f_id);

    /* Clipping of rain mass fraction 0 < y_p < 1 */
    cs_field_set_key_double(f, kscmin, 0.e0);
    cs_field_set_key_double(f, kscmax,  1.e0);

    /* Set the class index for the field */
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift: create additional mass flux.
     * This flux will then be reused for all scalars associated to this class
     * (here : injected liquid water variables)
     * We set the bit corresponding to drift flux computation to 1.
     * TODO (from old .f90 file) : make it optional ?*/
    int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;

    cs_field_set_key_int(f, keydri, drift);

    /* Set constant diffusivity for injected liquid mass fraction */
    int ifcvsl = -1;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    /* Set beta limiter to maintain y_p in the limits */
    eqp->isstpc = 2;
    eqp->blencv = 1.0;

    /* Transport and solve for the temperature of the liquid - with the same
     * drift as the mass fraction Y_l in the rain zones.
     * NB : Temperature of the liquid must be transported after the bulk
     * enthalpy. */

    f_id = cs_variable_field_create("y_p_t_l",
                                    "Yl.Tl rain",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    f = cs_field_by_id(f_id);
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift, but do not create an additional mass flux for the
     * enthalpy (use ^= to reset the bit for drift flux calculation).
     * It reuses the mass flux already identified with the mass fraction. */
    drift = CS_DRIFT_SCALAR_ON;

    cs_field_set_key_int(f, keydri, drift);

    /* Set variable diffusivity for the injected liquid enthalpy transport.
     * The diffusivity used in the transport equation will be the cell value
     * of the viscls array for field f */

    ifcvsl = 0;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    eqp = cs_field_get_equation_param(f);
    eqp->blencv = 1.0;

    /* Variable fields creation for rain drops velocities if we want to solve
     * rain fall velocity */
    cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
    if (ct_opt->solve_rain_velocity) {
      char f_name[80];
      char f_label[80];
      /* Rain drops velocities --> treated as particles */
      sprintf(f_name, "v_p_%02d", class_id);
      sprintf(f_label, "Vp_%02d", class_id);
      f_id = cs_variable_field_create(f_name, f_label,
                                      CS_MESH_LOCATION_CELLS, 3);
      f = cs_field_by_id(f_id);
      cs_field_set_key_int(f, keyccl, class_id);
      cs_add_model_field_indexes(f_id);

      /* Scalar with drift, but do not create an additional mass flux */
      drift = CS_DRIFT_SCALAR_ON;
      cs_field_set_key_int(f, keydri, drift);

      //TODO : Check equation parameters to set for v_p_ */
    }
  }

  {
    /* Packing zone variables */

    /* Associate injected liquid water in packing with class  */
    int class_id = 2;

    /* Mass fraction of liquid */
    int f_id = cs_variable_field_create("y_l_packing",
                                        "Yl packing",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    f = cs_field_by_id(f_id);

    /* Clipping of packing liquid mass fraction 0 < y_l_packing */
    cs_field_set_key_double(f, kscmin, 0.e0);

    /* Set the class index for the field */
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift: create additional mass flux.
     * This flux will then be reused for all scalars associated to this class
     * (here : injected liquid water variables)
     * We set the bit corresponding to drift flux computation to 1.
     * TODO (from old .f90 file) : make it optional ?*/
    int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX
                + CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;

    cs_field_set_key_int(f, keydri, drift);

    /* Set constant diffusivity for injected liquid mass fraction */
    int ifcvsl = -1;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    /* Upwind schemes for scalars in packing zone */
    eqp->blencv = 0.;
    eqp->idiff  = 0;
    eqp->idifft = 0;

    /* Do not show yl_packing in post-processing */
    cs_field_set_key_int(f, cs_field_key_id("post_vis"), 0);

    /* Transport and solve for the temperature of the liquid - with the same
     * drift as the mass fraction Y_l in the rain zones.
     * NB : Temperature of the liquid must be transported after the bulk
     * enthalpy. */

    f_id = cs_variable_field_create("enthalpy_liquid",
                                    "Enthalpy liq packing",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    /* TODO (from ctvarp.f90) : x_p_h_l or y_p_h_2 */

    f = cs_field_by_id(f_id);
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift, but do not create an additional mass flux for the
     * enthalpy (use ^= to reset the bit for drift flux calculation).
     * It reuses the mass flux already identified with the mass fraction. */
    drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;

    cs_field_set_key_int(f, keydri, drift);

    /* Set variable diffusivity for the injected liquid enthalpy transport.
     * The diffusivity used in the transport equation will be the cell value
     * of the viscls array for field f */

    ifcvsl = 0;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    eqp = cs_field_get_equation_param(f);
    /* Upwind schemes for scalars in packing zone */
    eqp->blencv = 0.;
    eqp->idiff  = 0;
    eqp->idifft = 0;
  }

  {
    /* Continuous phase variables */

    int class_id = -1;

    /* NB : 'c' stands for continuous and 'p' for particles */

    /* If not using the atmospheric module, we create the fields */
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != CS_ATMO_HUMID){
      /* Total mass fraction of water in the bulk humid air */
      int f_id = cs_variable_field_create("ym_water",
                                          "Ym water bulk",
                                          CS_MESH_LOCATION_CELLS,
                                          1);

      f = cs_field_by_id(f_id);

      /* Clipping : 0 < ym < 1 */
      cs_field_set_key_double(f, kscmin, 0.e0);
      cs_field_set_key_double(f, kscmax,  1.e0);

      /* Set the class index for the field */
      cs_field_set_key_int(f, keyccl, class_id);

      /* Set constant diffusivity for the dry air mass fraction.
       * The diffusivity used in the transport equation will be the cell value
       * of the" diffusivity_ref" for field f */
      int ifcvsl = -1;
      cs_field_set_key_int(f, kivisl, ifcvsl);

      /* Activate the drift for all scalars with key "drift" > 0 */
      int drift = CS_DRIFT_SCALAR_ON;

      /* Activated drift. As it is the continuous phase class (class_id = -1),
       * the convective flux is deduced for classes > 0
       * and bulk class (class_id = 0) */
      drift |= CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;
      cs_field_set_key_int(f, keydri, drift);

      cs_add_model_field_indexes(f->id);

      /* Equation parameters */
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      eqp->blencv = 1.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Add property fields
 *----------------------------------------------------------------------------*/

void
cs_ctwr_add_property_fields(void)
{
  cs_field_t *f;
  int class_id = 1;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  bool has_previous = false;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();


  /* Humid air properties */
  {
    /* Humidity field */
    f = cs_field_create("humidity",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Humidity");
  }

  {
    /* Saturated humidity field */
    f = cs_field_create("x_s",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Humidity sat");
  }

  {
    /* Relative humidity field */
    f = cs_field_create("x_rel",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Humidity rel");
  }


  {
    /* Humid air enthalpy field */
    f = cs_field_create("enthalpy",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Enthalpy humid air");
  }

  /* Liquid properties in packing */
  {
    /* Liquid temperature in packing */
    f = cs_field_create("temperature_liquid",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Temperature liq packing");
  }
  {
    /* True liquid mass fraction in packing */
    f = cs_field_create("y_liq_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Liq mass fraction packing");
  }


  {
    /* Liquid vertical velocity in packing */
    f = cs_field_create("vertvel_l",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Velocity liq packing");
  }

  {
    /* Liquid mass flux in packing */
    f = cs_field_create("mass_flux_l",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Mass flux liq packing");
  }

  /* Liquid properties in rain */
  {
    /* Rain temperature */
    f = cs_field_create("t_rain",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Temperature rain");
  }

  {
    /* True rain mass fraction */
    f = cs_field_create("y_rain",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Rain mass fraction");
  }

  /* Properties to create for rain velocity equation solving */
  if (ct_opt->solve_rain_velocity) {
    char f_name[80];

    /* Particle limit velocity */
    sprintf(f_name, "vg_lim_p_%02d", class_id);
    f = cs_field_create(f_name,
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    // cs_field_set_key_str(f, klbl, "Terminal velocity rain");
    // FIXME: labels should also be unique, so handle class id here

    /* Drift velocity for rain drops */
    sprintf(f_name, "vd_p_%02d", class_id);
    f = cs_field_create(f_name,
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    // cs_field_set_key_str(f, klbl, "Drift velocity rain");
    // FIXME: labels should also be unique, so handle class id here
  }

  /* Continuous phase properties */
  /* NB: 'c' stands for continuous and 'p' for particles */
  {
    /* Mass fraction of the continuous phase (X1) */
    f = cs_field_create("x_c",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Gas mass fraction");
  }

  {
    /* Mass fraction of the continuous phase (X1) BOUNDARY VALUE */
    f = cs_field_create("b_x_c",
                        field_type,
                        CS_MESH_LOCATION_BOUNDARY_FACES,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Boundary gas mass fraction");
  }

  /* Properties to create for rain velocity equation solving */
  if (ct_opt->solve_rain_velocity) {
    /* Continuous phase drift velocity */
    f = cs_field_create("vd_c",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Drift velocity gas phase");
  }

  /* Mass and energy exchange terms in packing and rain for evaporation rate
   * and thermal power post-processing, they are updated
   * in cs_ctwr_source_term */
  {
    /* Evaporation rate in packing */
    f = cs_field_create("evaporation_rate_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Evaporation packing");
  }

  {
    /* Evaporation rate in rain */
    f = cs_field_create("evaporation_rate_rain",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Evaporation rain");
  }

  {
    /* Thermal power in packing */
    f = cs_field_create("thermal_power_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Thermal power packing");
  }

  {
    /* Thermal power in rain */
    f = cs_field_create("thermal_power_rain",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Thermal power rain");
  }
}

/*----------------------------------------------------------------------------
 * Automatic boundary condition for cooling towers
 *----------------------------------------------------------------------------*/

void
cs_ctwr_bcond(void)
{
  /* Mesh-related data */
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const int *bc_type = cs_glob_bc_type;

  /* Fluid properties and physical variables */
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;

  cs_real_t *vel_rcodcl1 = CS_F_(vel)->bc_coeffs->rcodcl1;
  cs_field_t *y_rain= cs_field_by_name("y_p");
  cs_field_t *yt_rain= cs_field_by_name_try("y_p_t_l");
  cs_field_t *hlp= cs_field_by_name("enthalpy_liquid");
  cs_field_t *ylp = cs_field_by_name("y_l_packing");
  cs_field_t *yw = cs_field_by_name("ym_water");
  cs_field_t *t_h = cs_field_by_name("temperature");
  cs_real_t tkelvin = cs_physical_constants_celsius_to_kelvin;

  const cs_real_t xhum = air_prop->humidity0;
  const cs_real_t t0 =  cs_glob_fluid_properties->t0;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (bc_type[face_id] == CS_INLET || bc_type[face_id] == CS_FREE_INLET) {

      /* The turbulence BC values are calculated upstream using the base
         mechanism, so nothing specific is needed here. */

      /* Boundary conditions for the transported temperature of the humid air and
       * of the liquid water injected in the packing zones
       * --> Bulk values if not specified by the user
       * Assuming humid air is at conditions '0' */

      /* For humid air temperature */
      if (t_h->bc_coeffs->icodcl[face_id] == 0){
        t_h->bc_coeffs->icodcl[face_id] = 1;
        t_h->bc_coeffs->rcodcl1[face_id] = t0 - tkelvin;
      }

      /* For water mass fraction */
      if (yw->bc_coeffs->icodcl[face_id] == 0){
        yw->bc_coeffs->icodcl[face_id] = 1;
        yw->bc_coeffs->rcodcl1[face_id] = xhum / (1 + xhum);
      }

      /* For injected liquid in the packing*/
      if (ylp->bc_coeffs->icodcl[face_id] == 0){
        ylp->bc_coeffs->icodcl[face_id] = 1;
        ylp->bc_coeffs->rcodcl1[face_id] = 0.;
      }

      /* For injected liquid enthalpy in the packing*/
      if (hlp->bc_coeffs->icodcl[face_id] == 0){
        cs_real_t t_l = cs_glob_fluid_properties->t0 - tkelvin;
        cs_real_t h_l = cs_liq_t_to_h(t_l);

        /* Y_l . h_l is transported (not only h_l) */
        h_l *= ylp->bc_coeffs->rcodcl1[face_id];

        hlp->bc_coeffs->icodcl[face_id] = 1;
        hlp->bc_coeffs->rcodcl1[face_id] = h_l;
      }

    }

    /* For walls -> 0 flux for previous variables
     * Dirichlet condition y_rain = 0 to mimic water basin drain and avoid rain
     * accumulation on the floor */

    else if (   bc_type[face_id] == CS_SMOOTHWALL
             || bc_type[face_id] == CS_ROUGHWALL) {

      t_h->bc_coeffs->icodcl[face_id] = 3;
      t_h->bc_coeffs->rcodcl3[face_id] = 0.;
      yw->bc_coeffs->icodcl[face_id] = 3;
      yw->bc_coeffs->rcodcl3[face_id] = 0.;

      hlp->bc_coeffs->icodcl[face_id] = 3;
      hlp->bc_coeffs->rcodcl3[face_id] = 0.;
      ylp->bc_coeffs->icodcl[face_id] = 3;
      ylp->bc_coeffs->rcodcl3[face_id] = 0.;

      y_rain->bc_coeffs->icodcl[face_id] = 1;
      y_rain->bc_coeffs->rcodcl1[face_id] = 0.;
      if (yt_rain != NULL) {
        yt_rain->bc_coeffs->icodcl[face_id] = 1;
        yt_rain->bc_coeffs->rcodcl1[face_id] = 0.;
      }

    }
  }

  /* Extra variables to load if we solve rain velocity */

  const cs_ctwr_option_t *ct_opt = cs_glob_ctwr_option;
  if (ct_opt->solve_rain_velocity) {
    char f_name[80];
    int class_id = 1;
    sprintf(f_name, "v_p_%02d", class_id);
    cs_field_t *vp = cs_field_by_name(f_name);

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      for (cs_lnum_t i = 0; i < 3; i++){
        if (   bc_type[face_id] == CS_INLET
            || bc_type[face_id] == CS_FREE_INLET) {
          vp->bc_coeffs->icodcl[face_id] = 1;
          vp->bc_coeffs->rcodcl1[n_b_faces*i + face_id]
            = vel_rcodcl1[n_b_faces * i + face_id];
        }

        else if (   bc_type[face_id] == CS_SMOOTHWALL
                 || bc_type[face_id] == CS_ROUGHWALL) {
          vp->bc_coeffs->icodcl[face_id] = 1;
          vp->bc_coeffs->rcodcl1[n_b_faces*i + face_id] = 0.;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Initialize cooling towers fields, stage 0
 *----------------------------------------------------------------------------*/

void
cs_ctwr_fields_init0(void)
{
  int has_restart = cs_restart_present();
  cs_halo_t *halo = cs_glob_mesh->halo;

  /* Fluid properties and physical variables */
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;

  cs_field_t *t_h = cs_field_by_name("temperature");
  cs_field_t *ylp = cs_field_by_name("y_l_packing");
  cs_field_t *yw  = cs_field_by_name("ym_water");
  cs_field_t *hlp = cs_field_by_name("enthalpy_liquid");
  cs_field_t *tlp = CS_F_(t_l);

  cs_real_t tkelvin = cs_physical_constants_celsius_to_kelvin;
  const cs_real_t xhum = air_prop->humidity0;

  /* Only if the simulation is not a restart from another one */
  if (has_restart == 0){
    for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++){
      /* Humid air */
      t_h->val[cell_id] = fp->t0 - tkelvin;
      yw->val[cell_id] = xhum / (1. + xhum);

      /* Liquid in packing */
      tlp->val[cell_id] = t_h->val[cell_id];
      ylp->val[cell_id] = 0.;

    }
    if (halo != NULL) {
      cs_halo_sync_var(halo, CS_HALO_STANDARD, t_h->val);
      cs_halo_sync_var(halo, CS_HALO_STANDARD, yw->val);
      cs_halo_sync_var(halo, CS_HALO_STANDARD, tlp->val);
      cs_halo_sync_var(halo, CS_HALO_STANDARD, ylp->val);
    }

    /* Diffusivities of the dry air and the injected liquid
     * TODO : check if overwrites what users have specified */
    const int kvisl0 = cs_field_key_id("diffusivity_ref");

    cs_field_set_key_double(yw, kvisl0, 1.e-12);
    cs_field_set_key_double(ylp, kvisl0, 1.e-12);

    /* Initialize :
     * - the enthalpies, which are the solution variables
     * - the humidity, which users might have modified if they changed the
     *   mass fraction of the dry air in the humid air */

    cs_ctwr_init_field_vars(fp->ro0, fp->t0, fp->p0, air_prop->molmass_rat);

    if (air_prop->cp_l <= 0 || air_prop->lambda_l <= 0){
      bft_error(__FILE__,__LINE__, 0, _("Negative lambda or cp for liquid"));
    }

    else {
      cs_field_set_key_double(hlp, kvisl0, air_prop->lambda_l / air_prop->cp_l);
    }
  }

  else {
    /* TODO (from old ctiniv0 subroutine) Add restarts */
    const int kvisl0 = cs_field_key_id("diffusivity_ref");

    /* Diffusivities of the dry air and the injected liquid */
    cs_field_set_key_double(yw, kvisl0, 1.e-12);
    cs_field_set_key_double(ylp, kvisl0, 1.e-12);

    /* Restarts - recompute the required properties based on the saved solution
     * variables. For example : the humidity, liquid vertical velocity, etc. */
    cs_ctwr_restart_field_vars(fp->ro0, fp->t0, fp->p0, air_prop->humidity0,
                               air_prop->molmass_rat);
  }
}

/*----------------------------------------------------------------------------
 * Initialize cooling towers fields, stage 1
 *----------------------------------------------------------------------------*/

void
cs_ctwr_fields_init1(void)
{
  cs_halo_t *halo = cs_glob_mesh->halo;

  cs_field_t *t_h = cs_field_by_name("temperature");
  cs_field_t *ylp = cs_field_by_name("y_l_packing");
  cs_field_t *yw  = cs_field_by_name("ym_water");
  cs_field_t *tlp = CS_F_(t_l);

  /* Liquid inner mass flux */
  cs_lnum_t iflmas =
    cs_field_get_key_int(ylp, cs_field_key_id("inner_mass_flux_id"));
  cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  /* Liquid boundary mass flux */
  cs_lnum_t iflmab =
    cs_field_get_key_int(ylp, cs_field_key_id("boundary_mass_flux_id"));
  cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  cs_ctwr_init_flow_vars(i_mass_flux);

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, t_h->val);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, yw->val);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, tlp->val);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, ylp->val);
  }

  for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_b_faces; face_id++) {
    b_mass_flux[face_id] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Provide access to cs_ctwr_option
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
 * \param[in]  zone_criteria  zone selection criteria (or NULL)
 * \param[in]  z_id           z_id if zone already created (-1 otherwise)
 * \param[in]  zone_type      exchange zone type
 * \param[in]  delta_t        imposed delta temperature delta between inlet
 *                            and oulet of the zone
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
  cs_field_pointer_map(CS_ENUMF_(ym_w),
                       cs_field_by_name_try("ym_water"));
  cs_field_pointer_map(CS_ENUMF_(t_l),
                       cs_field_by_name_try("temperature_liquid"));
  cs_field_pointer_map(CS_ENUMF_(h_l),
                       cs_field_by_name_try("enthalpy_liquid"));
  cs_field_pointer_map(CS_ENUMF_(y_l_pack),
                       cs_field_by_name_try("y_l_packing"));
  cs_field_pointer_map(CS_ENUMF_(thermal_diff_h),
                       cs_field_by_name_try("thermal_conductivity"));
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
    ct->n_cells = cs_volume_zone_by_name(ct->name)->n_elts;
    ct->vol_f = cs_volume_zone_by_name(ct->name)->f_measure;
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
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;      /* Humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;      /* Humid air enthalpy */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;    /* Liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;    /* Liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val;   /* Liquid mass per unit
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
                         * y_l[cell_id_l] * liq_mass_flow[face_id];
      ct->h_l_in += sign * h_l[cell_id_l] * liq_mass_flow[face_id];
      ct->q_l_in += sign * y_l[cell_id_l] * liq_mass_flow[face_id];

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
        * y_l[cell_id_l] * liq_mass_flow[face_id];
      ct->q_l_out += sign * y_l[cell_id_l] * liq_mass_flow[face_id];
      ct->h_l_out += sign * h_l[cell_id_l] * liq_mass_flow[face_id];

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

  /* Initialize the fields - based on map */
  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val; /* Humid air (bulk)
                                                      density */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;     /* Humid air temperature */
  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre; /* Humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;     /* Humid air enthalpy */
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val;  /* Water mass fraction in
                                                      humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val; /* Humidity in air (bulk) */

  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;   /* Liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;   /* Liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass per unit */

  /* Packing zone liquid vertical velocity component */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;

  /* Rain drops variables */
  cs_field_t *cfld_yp = cs_field_by_name_try("y_p"); /* Rain mass fraction */
  cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_p");
  cs_field_t *cfld_drift_vel = cs_field_by_name_try("drift_vel_y_p");

  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

  cs_real_t *cpro_taup = NULL;
  if (cfld_taup != NULL)
    cpro_taup = cfld_taup->val;
  else
    BFT_MALLOC(cpro_taup, n_cells_with_ghosts, cs_real_t);

  const cs_air_fluid_props_t  *air_prop = cs_glob_air_props;
  cs_real_t rho_l = air_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t droplet_diam = air_prop->droplet_diam;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  /* Count clippings for rain / humidity variables */
  cs_gnum_t nclip_yw_min = 0;
  cs_gnum_t nclip_yw_max = 0;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field in case users have updated the initial
       dry air mass fraction.
       Note: this is a bit dubious as users could also have chosen
       to reset the humidity ? */

    /* Clippings of water mass fraction */
    if (y_w[cell_id] < 0.0){
      y_w[cell_id] = 0;
      nclip_yw_min += 1; //TODO : print it
    }

    if (y_w[cell_id] >= 1.0){
      y_w[cell_id] = 1. - cs_math_epzero;
      nclip_yw_max += 1; //TODO : print it
    }
    /* Note: the drops don't contribute to the bulk density yet */
    x[cell_id] = y_w[cell_id]/(1.0-y_w[cell_id]);

    /* Bulk humid air temperature */
    t_h[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
    t_h_a[cell_id] = t_h[cell_id];

    /* Update the humid air density */
    rho_h[cell_id] = cs_air_rho_humidair(x[cell_id],
                                         rho0,
                                         p0,
                                         t0,
                                         molmassrat,
                                         t_h[cell_id]);

    /* Update the humid air enthalpy */
    x_s[cell_id] = cs_air_x_sat(t_h[cell_id],p0);

    cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_air_h_humidair(cp_h,
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    /* Initialize the liquid vertical velocity component
     * this is correct for droplet and extended for other packing zones
     * NB: this value is derived from the drag coefficient:
     * C_D = 24 / Re * (1 + 0.15 * Re^0.687)
     * See ZOPLU HT-31-08-06 */

    cs_real_t v_lim =   cs_math_pow2(droplet_diam) * rho_l / (18. * visc)
                      * cs_math_3_norm(gravity);

    cs_real_t reynolds_old = 0.;
    cs_real_t reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;

    // FIXME make it global for the zone as restart...
    for (int sweep = 0;
         sweep < 100 && CS_ABS(reynolds - reynolds_old) > 0.001;
         sweep++) {
      reynolds_old = reynolds;
      v_lim =   pow(droplet_diam, 2.) * rho_l
              / (18. * visc * (1. + 0.15 * pow(reynolds, 0.687)))
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
  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, rho_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_taup);
  }

  /* Loop over exchange zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Initialize with the injection water temperature */
      t_l[cell_id] = ct->t_l_bc;

      /* Update the injected liquid enthalpy */
      h_l[cell_id] = cs_liq_t_to_h(t_l[cell_id]);

      /* Initialize the liquid vertical velocity component
       * this is correct for droplet and extended for other packing zones */
      vel_l[cell_id] = ct->v_liq_pack;

      /* Note that rho_h * Y_l * vel_l * Stot = q_l_bc */
      cs_real_t y_l_bc =   ct->q_l_bc
                         / (rho_h[cell_id] * vel_l[cell_id] * ct->surface);

      /* Initialize the liquid transported variables:
         liquid mass and enthalpy corrected by the density ratio */
      y_l[cell_id] = y_l_bc;

      /* The transported value is (y_l.h_l) and not (h_l) */
      h_l[cell_id] *= y_l[cell_id];
    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, vel_l);
    if (cfld_yp != NULL)
      cs_halo_sync_var(halo, CS_HALO_STANDARD, cfld_yp->val);
    if (cfld_drift_vel != NULL) {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, cfld_drift_vel->val, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, CS_HALO_STANDARD,
                                    cfld_drift_vel->val, 3);
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
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass fraction
                                                         in packing */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;      /* Liquid enthalpy
                                                         in packing */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;      /* Liquid temperature
                                                         in packing */

  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val; /* Humid air
                                                      (bulk) density */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val; /* Liquid velocity
                                                            in packing */

  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->i_face_normal;
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);

  const cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;

  const cs_halo_t *halo = cs_glob_mesh->halo;
  cs_lnum_t *packing_cell;

  /* Normalized gravity vector */

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  cs_real_t g_dir[3];
  cs_math_3_normalize(gravity, g_dir);

  /* Initialize the liquid mass flux to null */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
    liq_mass_flow[face_id] = 0.0;

  /* Tag and initialize the ct values in the packing zone cells */

  BFT_MALLOC(packing_cell, n_cells_with_ghosts, int);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++)
    packing_cell[cell_id] = -1;

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_MALLOC(ct->inlet_faces_ids, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_faces_ids, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_cells_ids, n_i_faces, cs_lnum_t);
    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
    for (int i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];
      if (ct->type != CS_CTWR_INJECTION)
        packing_cell[cell_id] = ict;
    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_untyped(halo, CS_HALO_STANDARD, sizeof(int), packing_cell);
  }

  /* Initialize the liquid mass flux at packing zone faces
   * and the ghost cells for the liquid mass and enthalpy
   * Initialize the couples (inlet faces, upwind cells) and
   * (outlet faces, upwind cells) arrays */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id_1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id_2 = i_face_cells[face_id][1];

    /* one of neigh. cells is in packing */
    if (packing_cell[cell_id_1] != -1 || packing_cell[cell_id_2] != -1) {

      int ct_id = CS_MAX(packing_cell[cell_id_1], packing_cell[cell_id_2]);
      cs_ctwr_zone_t *ct = _ct_zone[ct_id];

      /* Vertical (align with gravity) component of the surface vector */
      cs_real_t liq_surf = cs_math_3_dot_product(g_dir,
                                                 i_face_normal[face_id]);

      /* Face mass flux of the liquid */
      cs_lnum_t cell_id;
      if (liq_surf > 0.) { /* cell_id_1 is upwind cell for liq. flow */
        if (packing_cell[cell_id_1] != -1) /* cell_id_1 in the packing */
          cell_id = cell_id_1;
        else /* cell_id_1 in HALO of the packing and outside of it */
          cell_id = cell_id_2;
      }
      else { /* cell_id_2 is upwind cell for liq. flow */
        if (packing_cell[cell_id_2] != -1) /* cell_id_2 in the packing */
          cell_id = cell_id_2;
        else /* cell_id_2 in HALO of the packing and outside of it */
          cell_id = cell_id_1;
      }

      cs_real_t y_l_bc = ct->q_l_bc / (  rho_h[cell_id] * vel_l[cell_id]
                                       * ct->surface);
      liq_mass_flow[face_id] = rho_h[cell_id] * vel_l[cell_id] * liq_surf;

      /* Initialize a band of ghost cells on the top side of the
         packing zone in order to impose boundary values
         Take the upwind value for initialization */

      /* cell_id_1 in packing and not cell_id_2 */
      if (packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] == -1) {

        /* cell_id_2 is an inlet halo */
        if (liq_mass_flow[face_id] < 0.0) {

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
          y_l[cell_id_2] = y_l_bc;
          t_l[cell_id_2] = ct->t_l_bc;
          h_l[cell_id_2] = cs_liq_t_to_h(ct->t_l_bc);
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
          y_l[cell_id_1] = y_l_bc;
          t_l[cell_id_1] = ct->t_l_bc;
          h_l[cell_id_1] = cs_liq_t_to_h(ct->t_l_bc);
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
      }
      else if (   packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] >= 0
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
    }
    else {
      liq_mass_flow[face_id] = 0.0;
    }
  }

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_REALLOC(ct->inlet_faces_ids, ct->n_inlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_faces_ids, ct->n_outlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_cells_ids, ct->n_outlet_cells, cs_lnum_t);

    cs_parall_sum(1, CS_REAL_TYPE, &(ct->surface_in));
    cs_parall_sum(1, CS_REAL_TYPE, &(ct->surface_out));
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

  /* Initialize the fields - based on map */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;     /* Humid air (bulk) Cp */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* Humid air temperature */

  cs_real_t *t_h_a = (cs_real_t *)CS_F_(t)->val_pre; /* Humid air temperature
                                                        at previous time step */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* Humid air enthalpy */
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val;    /* Water mass fraction in
                                                        humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;     /* Saturated humidity */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;     /* Absolute humidity in
                                                        humid air (bulk) */

  /* Packing liquid quantities */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;      /* Liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;      /* Liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass per
                                                         unit cell volume */

  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val; /* Liquid vertical
                                                            velocity */

  /* Rain variables */
  cs_field_t *cfld_yp = cs_field_by_name_try("y_p"); /* Rain mass fraction */
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

  const cs_air_fluid_props_t  *air_prop = cs_glob_air_props;
  cs_real_t rho_l = air_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t droplet_diam = air_prop->droplet_diam;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  /* Recompute the initial values which were used in the initialization of
   * the calculation which is being restarted */
  cs_real_t y_w_ini = humidity0 / (1.0 + humidity0); // From 'ctiniv'
  if (y_w_ini < 0.0)
    y_w_ini = 0;

  if (y_w_ini >= 1.0)
    y_w_ini = 1. - cs_math_epzero;

  cs_real_t x_ini = y_w_ini/(1.0-y_w_ini);

  cs_real_t t_h_ini = t0 - cs_physical_constants_celsius_to_kelvin;

  cs_real_t rho_h_ini = cs_air_rho_humidair(x_ini,
                                            rho0,
                                            p0,
                                            t0,
                                            molmassrat,
                                            t_h_ini);

  /* Clip counters for humidity / water variables */
  int nclip_yw_min = 0;
  int nclip_yw_max = 0;

  /* Initialize the cooling towers variables */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field */

    /* Clippings of water mass fraction */
    if (y_w[cell_id] < 0.0){
      y_w[cell_id] = 0;
      nclip_yw_min += 1;
    }

    if (y_w[cell_id] >= 1.0){
      y_w[cell_id] = 1. - cs_math_epzero;
      nclip_yw_max += 1;
    }
    x[cell_id] = y_w[cell_id]/(1.0-y_w[cell_id]);

    /* Bulk humid air temperature at the reference temperature
       This is only calculated once at the beginning so same as
       'cs_ctwr_init_field_vars'
       No, this would be the value at the previous time step -
       At present, it is not stored in the restart file, so for lack
       of information initialize it with the present value of the temperature */
    t_h_a[cell_id] = t_h[cell_id];

    /* Update the humid air enthalpy based on the solved value of T_h */
    //FIXME Need to use the method of 'cs_ctwr_phyvar_update'

    x_s[cell_id] = cs_air_x_sat(t_h[cell_id],p0);

    cp_h[cell_id] = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_air_h_humidair(cp_h[cell_id],
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    /* Update the liquid temperature based on the solved liquid enthalpy
     * NB: May not be required as it is also done in 'cs_ctwr_phyvar_update'?
     * No, it must be done here because here we sweep over the entire
     * computational domain whereas 'cs_ctwr_phyvar_update' updates
     * T_l only over the packing zones */
    t_l[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
    if (y_l[cell_id] > 0.) {
      cs_real_t h_liq = h_l[cell_id] / y_l[cell_id];
      t_l[cell_id] = cs_liq_h_to_t(h_liq);
    }

    /* Initialize the liquid vertical velocity component
     * this is correct for droplet and extended for other packing zones
     * NB: this value is derived from the drag coefficient:
     * C_D = 24 / Re * (1 + 0.15 * Re^0.687)
     * See ZOPLU HT-31-08-06 */

    cs_real_t v_lim = pow(droplet_diam, 2.) * rho_l / (18. * visc)
                    * cs_math_3_norm(gravity);

    cs_real_t reynolds_old = 0.;

    /* Use the same humid air density which was used at the beginning of the
     * calculation being restarted, otherwise since rho_h changes during the
     * calculation, reynolds, v_lim and cpro_taup will end up being different
     * from the initial values used in the calculation being restarted */

    /* Droplet Reynolds number */
    //    cs_real_t reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
    cs_real_t reynolds = rho_h_ini * v_lim * droplet_diam / visc;

    for (int sweep = 0;
         sweep < 100 && fabs(reynolds - reynolds_old) > 0.001;
         sweep++) {
      reynolds_old = reynolds;
      v_lim =   pow(droplet_diam, 2.) * rho_l
              / (18. * visc * (1 + 0.15 * pow(reynolds, 0.687)))
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

  cs_gnum_t n_g_clip_yw_min = nclip_yw_min;
  cs_gnum_t n_g_clip_yw_max = nclip_yw_max;

  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yw_min);
  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yw_max);

  /* Printing clips in listing */
  if (n_g_clip_yw_min >= 1 || n_g_clip_yw_max >= 1) {
    bft_printf("WARNING : clipping on water mass fraction at restart in"
               "cs_ctwr_restart_field_vars : min_clip = %lu, max_clip = %lu\n",
                n_g_clip_yw_min, n_g_clip_yw_max);
  }

  /* Loop over exchange zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Initialize the liquid vertical velocity component
       * this is correct for droplet and extended for other packing zones */
      vel_l[cell_id] = ct->v_liq_pack;
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
        cs_halo_perio_sync_var_vect(halo, CS_HALO_STANDARD,
                                    cfld_drift_vel->val, 3);
    }
  }

  /* Free memory */
  if (cfld_taup == NULL)
    BFT_FREE(cpro_taup);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the thermo physical properties fields for the humid air and
 *        the liquid.
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_phyvar_update(cs_real_t  rho0,
                      cs_real_t  t0,
                      cs_real_t  p0)
{
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);
  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)(cs_glob_mesh->b_face_cells);
  const cs_halo_t *halo = cs_glob_mesh->halo;

  cs_air_fluid_props_t *air_prop = cs_glob_air_props;
  /* Water / air molar mass ratio */
  const cs_real_t molmassrat = air_prop->molmass_rat;

  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val;    /* Humid air
                                                         (bulk) density */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;      /* Humid air (bulk) Cp */

  /* Fields based on maps */
  cs_real_t *t_h = (cs_real_t *)CS_F_(t)->val;       /* Humid air temperature */
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* Humid air enthalpy */
  cs_real_t *therm_diff_h = cs_field_by_name("thermal_conductivity")->val;
  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;
  cs_real_t *bpro_x1 = cs_field_by_name("b_x_c")->val;
  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val;     /* Water mass fraction
                                                         in humid air */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;      /* Absolute humidity
                                                         in bulk air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;      /* Saturated humidity */
  cs_real_t *x_rel = cs_field_by_name("x_rel")->val;  /* Relative humidity */

  /* Packing zone variables */
  cs_real_t *t_l = (cs_real_t *)CS_F_(t_l)->val;      /* Liquid temperature */
  cs_real_t *h_l = (cs_real_t *)CS_F_(h_l)->val;      /* Liquid enthalpy */
  cs_real_t *y_l = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass per unit
                                                         cell volume*/
  cs_real_t *yl_pack = cs_field_by_name("y_liq_packing")->val; /* Liquid mass
                                                                * fraction in
                                                                * packing */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val; /* Liquid vertical
                                                            velocity */
  cs_real_t *mf_l = cs_field_by_name("mass_flux_l")->val; /* Liquid mass flux */

  cs_real_t *liq_mass_flow
    = cs_field_by_name("inner_mass_flux_y_l_packing")->val; //FIXME

  /* Variable and properties for rain zones */
  cs_field_t *cfld_yp = cs_field_by_name_try("y_p");   /* Rain scalar */
  cs_real_t *y_rain = cs_field_by_name("y_rain")->val; /* Yp times Tp */

  cs_real_t *yt_rain = cs_field_by_name("y_p_t_l")->val; /* Yp times Tp */
  cs_real_t *t_rain = cs_field_by_name("t_rain")->val; /* Rain temperature */

  cs_real_t *y_p = NULL;
  if (cfld_yp != NULL)
    y_p = cfld_yp->val;

  cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t lambda_h = cs_glob_air_props->lambda_h;

  /* Clipping counters for water / humidity variables */
  int nclip_yw_min = 0;
  int nclip_yw_max = 0;
  int nclip_yp_min = 0;
  int nclip_yp_max = 0;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Clippings of water mass fraction */
    if (y_w[cell_id] < 0.0){
      y_w[cell_id] = 0;
      nclip_yw_min += 1;
    }
    if (y_w[cell_id] >= 1.0){
      y_w[cell_id] = 1. - cs_math_epzero;
      nclip_yw_max += 1;
    }
    if (y_p != NULL) {
      if (y_p[cell_id] < 0.0){
        y_p[cell_id] = 0;
        nclip_yp_min += 1;
      }
      if ((y_p[cell_id] + y_w[cell_id]) >= 1.0){
        y_p[cell_id] = 1. - y_w[cell_id] - cs_math_epzero;
        nclip_yp_max += 1;
      }

      /* Recompute real rain mass fraction from Yp */
      if (y_p[cell_id] > 0.){
        y_rain[cell_id] = y_p[cell_id] / (1 + y_p[cell_id]);
      }

      /* Recompute real rain temperature from Yp.Tp */
      if (y_p[cell_id] > 1.e-4){
        t_rain[cell_id] = yt_rain[cell_id] / y_p[cell_id];
      }
      else {
        t_rain[cell_id] = 0.;
      }
    }

    /* Continuous phase mass fraction */
    cpro_x1[cell_id] = 1. - y_rain[cell_id];
    //TODO not one for rain zones - Why not?
    //If it represents the humid air, then it should be one?  If it represents
    //the dry air, then it should account for both y_p and y_w


    /* Update humidity field */
    x[cell_id] = y_w[cell_id]/(1.0-y_w[cell_id]);
    // FIXME for drops - This should be the proportion of 'gaseous' water
    // (dissolved and condensate) in the humid air:
    //   Y(dry air)+ Y(gasesous water) + Y(drops) = 1 in all computational cells
    //   If we do that, then the density needs to be revised as well and the
    //   temperatures of both the bulk (dry air + gaseous water +drops) and the
    //   humid air must be solved for.
    // Here, the approximation is that Y(drops) is negligible
    // This is NOT generally true : on MISTRAL, we reach Y(drops) > 0.5

    /* Saturated humidity */
    x_s[cell_id] = cs_air_x_sat(t_h[cell_id], p0);

    /*Relative humidity */
    x_rel[cell_id] = x[cell_id] / x_s[cell_id];


    /* Update the humid air temperature using new enthalpy but old
     * Specific heat */

    cp_h[cell_id] = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_air_h_humidair(cp_h[cell_id],
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    // Update the humid air enthalpy diffusivity lambda_h if solve for T_h?
    // Need to update since a_0 is variable as a function of T and humidity
    therm_diff_h[cell_id] = lambda_h;

    /* Update the humid air density */
    // Again, formally this should be the
    // bulk density, including the rain drops
    rho_h[cell_id] = cs_air_rho_humidair(x[cell_id],
                                         rho0,
                                         p0,
                                         t0,
                                         molmassrat,
                                         t_h[cell_id]);

  }

  cs_gnum_t n_g_clip_yw_min = nclip_yw_min;
  cs_gnum_t n_g_clip_yw_max = nclip_yw_max;

  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yw_min);
  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yw_max);

  cs_gnum_t n_g_clip_yp_min = nclip_yp_min;
  cs_gnum_t n_g_clip_yp_max = nclip_yp_max;

  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yp_min);
  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yp_max);

  /* Printing clips in listing */
  if (n_g_clip_yp_min >= 1 || n_g_clip_yp_max >= 1) {
    bft_printf("WARNING : clipping on rain mass fraction"
                "in cs_ctwr_phyvar_update : min_clip = %lu, max_clip = %lu\n",
                n_g_clip_yp_min, n_g_clip_yp_max);
  }
  if (n_g_clip_yw_min >= 1 || n_g_clip_yw_max >= 1) {
    bft_printf("WARNING : clipping on water mass fraction"
                "in cs_ctwr_phyvar_update : min_clip = %lu, max_clip = %lu\n",
                n_g_clip_yw_min, n_g_clip_yw_max);
  }

  /* If solving rain velocity */
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
  if (ct_opt->solve_rain_velocity) {
    int class_id = 1;
    char vg_lim_name[80];
    sprintf(vg_lim_name, "vg_lim_p_%02d", class_id);

    /* Drops terminal velocity fields */
    cs_field_t *vg_lim_p = cs_field_by_name(vg_lim_name);
    cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
      cs_glob_physical_constants->gravity[1],
      cs_glob_physical_constants->gravity[2]};
    cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_p");
    cs_real_t *cpro_taup = NULL;
    if (cfld_taup != NULL)
      cpro_taup = cfld_taup->val;

    /* Continuous phase drift velocity */
    cs_field_t *vd_c = cs_field_by_name("vd_c");

    /* Rain drift velocity variables */
    char f_name[80];
    sprintf(f_name, "vd_p_%02d", class_id);
    cs_field_t *vd_p = cs_field_by_name(f_name);
    sprintf(f_name, "v_p_%02d", class_id);
    cs_field_t *vp = cs_field_by_name(f_name);
    cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t i = 0; i < 3; i++) {

      /* Drops terminal velocity */
      vg_lim_p->val[cell_id * 3 + i] = cpro_taup[cell_id] * gravity[i];

      /* Continuous phase drift velocity */
      vd_c->val[cell_id * 3 + i] = 0.;

      /* Rain drops drift velocity calculation */
      if (y_p[cell_id] > 1.e-7) {
        vd_p->val[cell_id*3 + i] = vp->val[cell_id*3 + i] - vel[cell_id][i];
      }
      else {
        vd_p->val[cell_id * 3 + i] = 0.;
      }
      vd_c->val[cell_id * 3 + i] -= y_p[cell_id]
        * vd_p->val[cell_id * 3 + i];
      vd_c->val[cell_id * 3 + i] /= cpro_x1[cell_id];
      }
    }
  }

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < _n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;

    /* Packing zone */
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Update the injected liquid temperature
       * NB: (y_l.h_l) is transported and not (h_l) */
      if (y_l[cell_id] > 0.) {
        cs_real_t h_liq = h_l[cell_id] / y_l[cell_id];
        t_l[cell_id] = cs_liq_h_to_t(h_liq);
        mf_l[cell_id] = y_l[cell_id] * rho_h[cell_id] * vel_l[cell_id];
        yl_pack[cell_id] = y_l[cell_id] / (1 + y_l[cell_id]);
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
        cs_lnum_t cell_id_l;

        /* Convention: outlet is positive mass flux
         * Then upwind cell for liquid is i_face_cells[][0] */
        int sign = 1;
        if (liq_mass_flow[face_id] < 0) {
          sign = -1;
          cell_id_l = i_face_cells[face_id][1];
        }
        else {
          cell_id_l = i_face_cells[face_id][0];
        }

        /* h_l is in fact (y_l. h_l),
         * and the transport field is (y_l*liq_mass_flow) */
        ct->t_l_out += sign * t_l[cell_id_l]
          * y_l[cell_id_l] * liq_mass_flow[face_id];
        ct->q_l_out += sign * y_l[cell_id_l] * liq_mass_flow[face_id];
      }

      cs_parall_sum(1, CS_REAL_TYPE, &(ct->t_l_out));
      cs_parall_sum(1, CS_REAL_TYPE, &(ct->q_l_out));

      ct->t_l_out /= ct->q_l_out;

      /* Relaxation of ct->t_l_bc */
      ct->t_l_bc = (1. - ct->relax) * ct->t_l_bc
                 + ct->relax * (ct->t_l_out + ct->delta_t);

      /* Clipping between 0 and 100 */
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
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_source_term(int              f_id,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)(m->i_face_cells);
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;

  /* Water / air molar mass ratio */
  const cs_real_t molmassrat = air_prop->molmass_rat;

  cs_real_t  *rho_h = (cs_real_t *)CS_F_(rho)->val; /* Humid air
                                                       (bulk) density */
  cs_real_3_t *vel_h = (cs_real_3_t *)CS_F_(vel)->val; /* Humid air (bulk)
                                                          velocity*/

  cs_real_t *y_w = (cs_real_t *)CS_F_(ym_w)->val; /* Water mass fraction
                                                     in humid air */

  cs_real_t *t_h = cs_field_by_name("temperature")->val; /* Humid air
                                                            temperature */
  cs_real_t *t_l = cs_field_by_name("temperature_liquid")->val;
  cs_real_t *x = cs_field_by_name("humidity")->val; /* Humidity in air (bulk) */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;  /* Liquid velocity
                                                             in packing */
  cs_real_t *y_l = CS_F_(y_l_pack)->val;

  /* Variable and properties for rain drops */
  cs_field_t *cfld_yp = cs_field_by_name("y_p");     /* Rain mass fraction */
  cs_field_t *cfld_yt_rain = cs_field_by_name("y_p_t_l"); /* Yp times Tp */
  cs_field_t *cfld_drift_vel = cs_field_by_name("drift_vel_y_p"); /* Rain drift
                                                                     velocity */

  /* Rain inner mass flux */
  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  cs_real_t *imasfl_r = cs_field_by_id(cs_field_get_key_int(cfld_yp, kimasf))->val;

  cs_real_t vertical[3], horizontal[3];

  const cs_ctwr_option_t *ct_opt = cs_glob_ctwr_option;

  int evap_model = ct_opt->evap_model;

  /* Need to cook up the cell value of the liquid mass flux
     In the old code, it seems to be taken as the value of the
     face mass flux upstream of the cell */

  cs_real_t v_air = 0.;

  cs_real_t mass_flux_h = 0.; /* Highly suspicious for rain
                                 zones - not recomputed */

  /* Identify the source term formulation for the required field */

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_real_t *f_var = f->val;  /* Field variable */

  /* Compute the source terms */

  /* Fields for source terms post-processing */
  cs_real_t *evap_rate_pack = NULL;
  cs_real_t *evap_rate_rain = NULL;
  evap_rate_pack = cs_field_by_name("evaporation_rate_packing")->val;
  evap_rate_rain = cs_field_by_name("evaporation_rate_rain")->val;

  cs_real_t *thermal_power_pack = NULL;
  cs_real_t *thermal_power_rain = NULL;
  thermal_power_pack = cs_field_by_name("thermal_power_packing")->val;
  thermal_power_rain = cs_field_by_name("thermal_power_rain")->val;

  /* Table to track cells belonging to packing zones */
  cs_lnum_t  *packing_cell;
  BFT_MALLOC(packing_cell, m->n_cells_with_ghosts, int);
#   pragma omp parallel for if (m->n_cells> CS_THR_MIN)
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    packing_cell[cell_id] = -1;

  /* Air / fluid properties */
  cs_real_t cp_v = air_prop->cp_v;
  cs_real_t cp_l = air_prop->cp_l;
  cs_real_t hv0 = air_prop->hv0;
  cs_real_t rho_l = air_prop->rho_l;
  cs_real_t visc = fp->viscl0;
  cs_real_t  p0 = fp->p0;
  cs_real_t lambda_h = air_prop->lambda_h;
  cs_real_t droplet_diam  = air_prop->droplet_diam;
  cs_real_t sigma = air_prop->sigma;

  if (evap_model != CS_CTWR_NONE) {

    /* =========================
       Phase change source terms
       ========================= */

    cs_math_3_normalize(cs_glob_physical_constants->gravity, vertical);

    vertical[0] *= -1;
    vertical[1] *= -1;
    vertical[2] *= -1;

    horizontal[0] = vertical[0] -1.;
    horizontal[1] = vertical[1] -1.;
    horizontal[2] = vertical[2] -1.;

    /* Phase change in the packing zones
       Between the liquid film and the humid air
       ========================================= */

    for (int ict = 0; ict < _n_ct_zones; ict++) {

      cs_ctwr_zone_t *ct = _ct_zone[ict];

      /* We skip this if we are in injection zone */
      if (ct->type == CS_CTWR_INJECTION)
        continue;

      /* Packing zone characteristics */
      cs_real_t a_0 = ct->xap;
      cs_real_t xnp = ct->xnp;
      int zone_type = ct->type;

      const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;

      for (cs_lnum_t j = 0; j < ct->n_cells; j++) {

        cs_lnum_t cell_id = ze_cell_ids[j];
        /* Identify packing cells ids */
        if (ct->type != CS_CTWR_INJECTION)
          packing_cell[cell_id] = ict;

        /* For correlations, T_h cannot be greater than T_l */
        cs_real_t temp_h = CS_MIN(t_h[cell_id], t_l[cell_id]);

        /* Saturation humidity at humid air temperature */
        cs_real_t x_s_th = cs_air_x_sat(temp_h, p0);

        /* Saturation humidity at injected liquid temperature */
        cs_real_t x_s_tl = cs_air_x_sat(t_l[cell_id], p0);

        /*--------------------------------------------*
         * Counter or cross flow packing zone         *
         *--------------------------------------------*/

        if (zone_type == CS_CTWR_COUNTER_CURRENT) {
          /* Counter flow packing */
          v_air = CS_ABS(cs_math_3_dot_product(vel_h[cell_id], vertical));
        }
        else if (zone_type == CS_CTWR_CROSS_CURRENT) {
          /* Cross flow packing */
          v_air = CS_ABS(cs_math_3_dot_product(vel_h[cell_id], horizontal));
        }

        /* Dry air flux */
        mass_flux_h = rho_h[cell_id] * v_air * (1. - y_w[cell_id]);

        /* Liquid mass flux */
        cs_real_t mass_flux_l = rho_h[cell_id] * y_l[cell_id] * vel_l[cell_id];

        /* Evaporation coefficient 'Beta_x' (kg/m2/s) times exchange surface
         * per unit of volume 'ai' (m2/m3)*/
        cs_real_t beta_x_ai = 0.;
        /* There is evaporation only if we have an injected liquid flow */
        if (mass_flux_l > 0.){
          beta_x_ai = a_0*mass_flux_l*pow((mass_flux_h/mass_flux_l), xnp);
        }

        /* Source terms for the different equations */

        /* Humid air mass source term */
        cs_real_t mass_source = 0.0;
        if (x[cell_id] <= x_s_th) {
          mass_source = beta_x_ai*(x_s_tl - x[cell_id]);
        }
        else {
          mass_source = beta_x_ai*(x_s_tl - x_s_th);
        }
        mass_source = CS_MAX(mass_source, 0.);

        cs_real_t vol_mass_source = mass_source * cell_f_vol[cell_id];
        cs_real_t vol_beta_x_ai = beta_x_ai * cell_f_vol[cell_id];

        /* Humid air thermal source term */
        cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

        /* Global mass source term for continuity (pressure) equation
         * Note that rain is already considered in the bulk, so inner
         * mass transfer between liquid and vapor disappears */
        if (f_id == (CS_F_(p)->id)) {
          /* Warning: not multiplied by Cell volume! no addition neither */
          exp_st[cell_id] = mass_source;

          /* Saving evaporation rate for post-processing */
          evap_rate_pack[cell_id] = mass_source;
        }

        /* Water mass fraction equation except rain */
        else if (f_id == (CS_F_(ym_w)->id)) {
          //TODO add mass_from_rain
          exp_st[cell_id] += vol_mass_source * (1. - f_var[cell_id]);
          imp_st[cell_id] += vol_mass_source;
        }

        /* Injected liquid mass equation (solve in drift model form) */
        else if (f_id == (CS_F_(y_l_pack)->id)) {
          exp_st[cell_id] -= vol_mass_source * y_l[cell_id];
          imp_st[cell_id] += vol_mass_source;
        }

        /* Humid air temperature equation */
        else if (f_id == (CS_F_(t)->id)) {
          // FIXME source term for theta_l instead...
          /* Because the writing is in a non-conservative form */
          cs_real_t l_imp_st = vol_mass_source * cp_h;
          cs_real_t l_exp_st = 0.;
          cs_real_t xlew = _lewis_factor(evap_model, molmassrat,
                                         x[cell_id], x_s_tl);
          if (x[cell_id] <= x_s_th) {
            /* Implicit term */
            l_imp_st += vol_beta_x_ai * (xlew * cp_h
                                          + (x_s_tl - x[cell_id]) * cp_v
                                            / (1. + x[cell_id]));
            l_exp_st = l_imp_st * (t_l[cell_id] - f_var[cell_id]);
          }
          else {
            cs_real_t coeft = xlew * cp_h;
            /* Implicit term */
            l_imp_st += vol_beta_x_ai * (coeft + (x_s_tl - x_s_th) * cp_l
                                                 / (1. + x[cell_id]));
            /* Explicit term */
            //FIXME : why does a hv0 term appears while we work on temperature ?
            l_exp_st = vol_beta_x_ai * (coeft * t_l[cell_id]
                                        + (x_s_tl - x_s_th)
                                          * (cp_v * t_l[cell_id] + hv0)
                                          / (1. + x[cell_id]))
                       - l_imp_st * f_var[cell_id];
          }
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
          exp_st[cell_id] += l_exp_st;
        }

        /* Injected liquid enthalpy equation (solve in drift model form)
         * NB: it is in fact "y_l x h_l" */
        else if (f_id == (CS_F_(h_l)->id)) {
          /* Liquid temperature in Kelvin */
          cs_real_t t_l_k = t_l[cell_id]
                            + cs_physical_constants_celsius_to_kelvin;
          cs_real_t l_exp_st = 0.;
          /* Implicit term */
          cs_real_t l_imp_st = vol_mass_source;
          l_exp_st -= l_imp_st * f_var[cell_id];
          cs_real_t xlew = _lewis_factor(evap_model,
                                         molmassrat,
                                         x[cell_id],
                                         x_s_tl);
          /* Under saturated */
          if (x[cell_id] <= x_s_th) {
            /* Explicit term */
            l_exp_st -= vol_beta_x_ai * ((x_s_tl - x[cell_id])
                                         * (cp_v * t_l_k + hv0)
                                         + xlew * cp_h
                                           * (t_l[cell_id] - t_h[cell_id]));
          }
          /* Over saturated */
          else {
            cs_real_t coefh = xlew * cp_h;
            /* Explicit term */
            l_exp_st += vol_beta_x_ai * (coefh * (t_h[cell_id] - t_l[cell_id])
                                         + (x_s_tl - x_s_th) / (1. + x[cell_id])
                                           * (cp_l * t_h[cell_id]
                                              - (cp_v * t_l[cell_id] + hv0)));
          }
          /* Because we deal with an increment */
          exp_st[cell_id] += l_exp_st;
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);

          /* Saving thermal power for post-processing */
          thermal_power_pack[cell_id] = -(l_exp_st + l_imp_st * f_var[cell_id])
                                         / cell_f_vol[cell_id];
        }
      } /* end loop over the cells of a packing zone */

    } /* end loop over all the packing zones */

    /* Phase change in the whole domain
       Between the rain drops and the humid air
       ======================================== */

    if (cfld_yp != NULL) {
      cs_real_t *y_rain = (cs_real_t *)cfld_yp->val;
      cs_real_t *temp_rain = (cs_real_t *)cs_field_by_name("t_rain")->val;

      for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

        if (y_rain[cell_id] > 0.) {

          /* For correlations, T_h cannot be greater than T_p */
          cs_real_t temp_h = CS_MIN(t_h[cell_id], temp_rain[cell_id]);

          /* Saturation humidity at the temperature of the humid air */
          cs_real_t x_s_th = cs_air_x_sat(temp_h, p0);

          /* Saturation humidity at the temperature of the rain drop  */
          cs_real_t x_s_tl = cs_air_x_sat(temp_rain[cell_id], p0);

          cs_real_3_t *drift_vel_rain
            = (cs_real_3_t *restrict)(cfld_drift_vel->val);
          cs_real_t drift_vel_mag = cs_math_3_norm(drift_vel_rain[cell_id]);
          cs_real_t xlew = _lewis_factor(evap_model,
                                         molmassrat,
                                         x[cell_id],
                                         x_s_tl);
          cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

          /* Rain droplets Reynolds number */
          cs_real_t rey = rho_h[cell_id] * drift_vel_mag * droplet_diam / visc;

          /* Prandtl number */
          cs_real_t pr = cp_h * visc / lambda_h;

          /* Nusselt number correlations */
          /* Ranz-Marshall or Hughmark when rey <= 776.06 && pr <= 250. */
          cs_real_t nusselt = 2.+0.6*sqrt(rey)*pow(pr,(1./3.));
          /* Hughmark when rey > 776.06 && pr <= 250. */
          if (rey > 776.06 && pr <= 250.) {
            nusselt = 2. + 0.27*pow(rey, 0.62)*pow(pr,(1./3.));
          }

          /* Convective exchange coefficient 'a_c' */
          cs_real_t a_c = (nusselt * lambda_h) / droplet_diam;

          /* beta_x coefficient */
          cs_real_t beta_x = a_c / (xlew * cp_h);

          /* Exchange surface area per unit volume based on the total droplets
           * surface in the cell
           * NOTE: Use rho_h to compute the number of particles per unit volume
           * since conservation equation for Y_p based on rho_h
           *   --> Should really be rho_mixture!?
           * Use the symmetric relationship:
           *   a_i = 6*alpha_p*(1.-alpha_p)/droplet_diam
           * where alpha_p is the droplets volume fraction
           * - this kills transfer when there is only one phase (pure humid air
           *   or pure rain) */
          cs_real_t vol_frac_rain = y_rain[cell_id] * rho_h[cell_id] / rho_l;
          if (vol_frac_rain >= 1.0)
            vol_frac_rain = 1.0;
          cs_real_t a_i =  6.0 * vol_frac_rain * (1.0 - vol_frac_rain)
                           / droplet_diam;

          /* Evaporation coefficient 'Beta_x' times exchange surface 'ai' */
          cs_real_t beta_x_ai = beta_x * a_i;

          /* Source terms for the different equations */

          /* Humid air mass source term */
          cs_real_t mass_source = 0.0;
          if (x[cell_id] <= x_s_th) {
            mass_source = beta_x_ai * (x_s_tl - x[cell_id]);
          }
          else {
            mass_source = beta_x_ai * (x_s_tl - x_s_th);
          }
          mass_source = CS_MAX(mass_source, 0.);

          cs_real_t vol_mass_source = mass_source * cell_f_vol[cell_id];
          cs_real_t vol_beta_x_ai = beta_x_ai * cell_f_vol[cell_id];

          /* Humid air thermal source term */
          cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);
          cs_real_t le_f = _lewis_factor(evap_model, molmassrat,
                                         x[cell_id], x_s_tl);
          cs_real_t st_th = 0.;
          if (x[cell_id] <= x_s_th){
            st_th = beta_x_ai * ((x_s_tl - x[cell_id]) / (1 + x[cell_id])
                                 * (cp_v * temp_rain[cell_id] + hv0)
                           + le_f * cp_h * (temp_rain[cell_id] - t_h[cell_id]));
          }
          else{
            st_th = beta_x_ai * ((x_s_tl - x_s_th) / (1 + x[cell_id])
                                 * (cp_v * temp_rain[cell_id] + hv0)
                            + le_f * cp_h * (temp_rain[cell_id] - t_h[cell_id]))
                                 + le_f * (x[cell_id] - x_s_th)
                                   * ((cp_l - cp_v) * temp_rain[cell_id] + hv0);
          }

          /* Global mass source term for continuity (pressure) equation
           * Note that if rain were already considered in the bulk, then inner
           * mass transfer between liquid and vapor would disappear */

          //FIXME: Consider putting rain in the bulk
          //       --> implication on 'c' variables different from 'h'
          if (f_id == (CS_F_(p)->id)) {
            /* Warning: not multiplied by Cell volume! no addition neither */
            // FIXME: Addition needed to avoid deleting packing mass
            //        source term ?
            exp_st[cell_id] += mass_source;

            /* Saving evaporation rate for post-processing */
            evap_rate_rain[cell_id] = mass_source;
          }

          /* Water (vapor + condensate) in gas mass fraction equation
             except rain */
          else if (f_id == (CS_F_(ym_w)->id)) {
            exp_st[cell_id] += vol_mass_source*(1. - f_var[cell_id]);
            imp_st[cell_id] += vol_mass_source;
          }

          /* Rain drop mass equation (solve in drift model form) */
          else if (f_id == cfld_yp->id) {
            exp_st[cell_id] -= vol_mass_source * y_rain[cell_id];
            imp_st[cell_id] += vol_mass_source;
          }

          /* Humid air temperature equation */
          else if (f_id == (CS_F_(t)->id)) {
            /* Because the writing is in a non-conservative form */
            cs_real_t l_imp_st = vol_mass_source * cp_h;
            cs_real_t l_exp_st = 0.;
            if (x[cell_id] <= x_s_th) {
              /* Implicit term */
              l_imp_st += vol_beta_x_ai * (xlew * cp_h
                                           + (x_s_tl - x[cell_id]) * cp_v
                                             / (1. + x[cell_id]));
              /* Explicit term */
              l_exp_st = l_imp_st * (temp_rain[cell_id] - f_var[cell_id]);
            }
            else {
              cs_real_t coeft = xlew * cp_h;
              /* Implicit term */
              l_imp_st += vol_beta_x_ai * (coeft + (x_s_tl - x_s_th)
                                                   * cp_l / (1. + x[cell_id]));
              /* Explicit term */
              l_exp_st = vol_beta_x_ai * (coeft * temp_rain[cell_id]
                         + (x_s_tl - x_s_th) * (cp_v * temp_rain[cell_id] + hv0)
                           / (1. + x[cell_id])) - l_imp_st * f_var[cell_id];
            }
            imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
            exp_st[cell_id] += l_exp_st;
          }

          /* Rain temperature equation (solve in drift model form)
           * NB: it should be in fact "y_rain x T_rain" */  //FIX ME
          else if (f_id == cfld_yt_rain->id) {
            /* Implicit term */
            //          cs_real_t l_imp_st = vol_mass_source * cp_l;
            cs_real_t l_imp_st = vol_mass_source;
            if (x[cell_id] <= x_s_th) {
              cs_real_t coefh = vol_beta_x_ai * (xlew * cp_h
                                                 + (x_s_tl - x[cell_id]) * cp_v
                                                 / (1. + x[cell_id]));
              /* Note: Rain temperature is currently not treated as a
               * temperature field -> division by cp_l needed for implicit and
               * explicit source terms */
              coefh /= cp_l;
              exp_st[cell_id] += coefh * (t_h[cell_id] - temp_rain[cell_id]);
            }
            else {
              cs_real_t coefh = xlew * cp_h;
              /* Note: Rain temperature is currently not treated as a
               * temperature field -> division by cp_l needed for implicit and
               * explicit source terms */
              coefh /= cp_l;
              //FIXME: why does a hv0 term appear
              //       while we work on temperature ?
              exp_st[cell_id] += vol_beta_x_ai * (coefh * (t_h[cell_id]
                                                  - temp_rain[cell_id])
                                         + (x_s_tl - x_s_th) / (1. + x[cell_id])
                                           * (cp_l * t_h[cell_id]
                                  - (cp_v * temp_rain[cell_id] + hv0))) / cp_l;
            }
            /* Because we deal with an increment */
            exp_st[cell_id] -= l_imp_st * f_var[cell_id];
            imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
            /* Note: y_p_t_l not treated as temperature field
             * so no multiplication by cp_l
             * */

            /* Saving thermal power for post-processing */
            if (temp_rain[cell_id] > 0.){
              thermal_power_rain[cell_id] = st_th;
            }
          }

        } /* End if (y_rain > 0) */

      } /* End loop over all the cells of the domain */
    }
  } /* End evaporation model active */

  if (ct_opt->has_rain) {

    /* Generate rain from packing zones which are leaking
       ================================================== */

    cs_real_t *liq_mass_frac = CS_F_(y_l_pack)->val; /* Liquid mass fraction
                                                        in packing */
    /* Inner mass flux of liquids (in the packing) */
    cs_real_t *liq_mass_flow
      = cs_field_by_name("inner_mass_flux_y_l_packing")->val;

    for (int ict = 0; ict < _n_ct_zones; ict++) {

      cs_ctwr_zone_t *ct = _ct_zone[ict];

      if (ct->xleak_fac > 0.0 && ct->type != CS_CTWR_INJECTION) {

        /* Rain generation source terms
           ============================ */

        for (cs_lnum_t i = 0; i < ct->n_outlet_faces; i++) {

          /* Leak face_id */
          cs_lnum_t face_id = ct->outlet_faces_ids[i];
          cs_lnum_t cell_id_leak, cell_id_rain;

          /* Convention: outlet is positive mass flux
           * Then upwind cell for liquid is i_face_cells[][0] */
          int sign = 1;
          if (liq_mass_flow[face_id] < 0) {
            sign = -1;
            cell_id_leak = i_face_cells[face_id][1];
            cell_id_rain = i_face_cells[face_id][0];
          }
          else {
            cell_id_leak = i_face_cells[face_id][0];
            cell_id_rain = i_face_cells[face_id][1];
          }

          /* Note: vol_mass_source must not be multiplied by
           * cell_f_vol[cell_id_rain]
           * because mass source computed from liq_mass_flow is
           * already in kg/s associated to the facing rain cell */
          cs_real_t vol_mass_source = ct->xleak_fac
            * liq_mass_frac[cell_id_leak] * sign * liq_mass_flow[face_id];

          /* Global bulk mass - continuity */
          //FIXME - Ignore for now because drops are not in the bulk
          //        if (f_id == (CS_F_(p)->id)) {
          //          /* Warning: not multiplied by Cell volume! */
          //          exp_st[cell_id] = mass_source;
          //        }

          /* Injected liquid mass equation for rain zones
             (solve in drift model form) */
          if (f_id == cfld_yp->id) {
            /* Because we deal with an increment */
            exp_st[cell_id_rain] += vol_mass_source
                                    * (1. - f_var[cell_id_rain]);
            imp_st[cell_id_rain] += vol_mass_source;
          }
          /* Rain temperature */
          else if (f_id == cfld_yt_rain->id) {
            // FIXME: There should be a y_p factor in there so that
            // mass and enthalpy are compatible
            /* The transported variable is y_rain * T_rain */
            /* Since it is treated as a scalar, no multiplication by cp_l is
             * required */
            /* For temperature equation of the rain */
            exp_st[cell_id_rain] += vol_mass_source
                                    * (t_l[cell_id_leak] - f_var[cell_id_rain]);
            imp_st[cell_id_rain] += vol_mass_source;
          }

        } /* End of loop through outlet cells of the packing zone */

      } /* End of leaking zone test */

      /* Testing if we are in an rain injection zone */
      else if (ct->xleak_fac > 0.0 && ct->type == CS_CTWR_INJECTION) {
        const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
        cs_real_t inj_vol = ct->vol_f;

        for (cs_lnum_t j = 0; j < ct->n_cells; j++) {
          cs_lnum_t cell_id = ze_cell_ids[j];

          cs_real_t vol_mass_source = cell_f_vol[cell_id] / inj_vol
                                      * ct->q_l_bc * ct->xleak_fac;
          cs_real_t t_inj = ct->t_l_bc;

          /* Injected liquid mass equation for rain zones
             (solve in drift model form) */
          if (f_id == cfld_yp->id) {
            /* Because we deal with an increment */
            exp_st[cell_id] += vol_mass_source * (1. - f_var[cell_id]);
            imp_st[cell_id] += vol_mass_source;
          }
          /* Rain temperature */
          else if (f_id == cfld_yt_rain->id) {
            // FIXME: There should be a y_p factor in there so that
            // mass and enthalpy are compatible
            /* The transported variable is y_rain * T_rain */
            /* Since it is treated as a scalar, no multiplication by cp_l is
             * required */
            /* For temperature equation of the rain */
            exp_st[cell_id] += vol_mass_source * (t_inj - f_var[cell_id]);
            imp_st[cell_id] += vol_mass_source;
          }
        }
      }
    } /* End of loop through the packing zones */

    cs_real_t *y_rain = (cs_real_t *)cfld_yp->val;
    cs_real_t *temp_rain = (cs_real_t *)cs_field_by_name("t_rain")->val;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      cs_lnum_t cell_id_0 = i_face_cells[face_id][0];
      cs_lnum_t cell_id_1 = i_face_cells[face_id][1];

      /* one of neigh. cells is in packing */
      if (packing_cell[cell_id_0] != -1 || packing_cell[cell_id_1] != -1) {

        int ct_id = CS_MAX(packing_cell[cell_id_0], packing_cell[cell_id_1]);
        cs_ctwr_zone_t *ct = _ct_zone[ct_id];

        /* Rain sink term in packing zones */
        //TODO : Add rain leak portion inside packing
        if (f_id == cfld_yp->id) {
          if (packing_cell[cell_id_0] != -1) {
            imp_st[cell_id_0] += CS_MAX(imasfl_r[face_id], 0.);
            exp_st[cell_id_0] -= CS_MAX(imasfl_r[face_id],0) * f_var[cell_id_0];
          }
          if (packing_cell[cell_id_1] != -1) {
            imp_st[cell_id_1] += CS_MAX(-imasfl_r[face_id], 0.);
            exp_st[cell_id_1] -= CS_MAX(-imasfl_r[face_id],0) * f_var[cell_id_1];
          }
        }

        if (f_id == cfld_yt_rain->id) {
          if (packing_cell[cell_id_0] != -1) {
            imp_st[cell_id_0] += CS_MAX(imasfl_r[face_id], 0.);
            exp_st[cell_id_0] -= CS_MAX(imasfl_r[face_id],0) * f_var[cell_id_0];
          }
          if (packing_cell[cell_id_1] != -1) {
            imp_st[cell_id_1] += CS_MAX(-imasfl_r[face_id], 0.);
            exp_st[cell_id_1] -= CS_MAX(-imasfl_r[face_id],0) * f_var[cell_id_1];
          }
        }

        /* Liquid source term in packing zones from rain */
        if (f_id == CS_F_(y_l_pack)->id) {
          if (packing_cell[cell_id_0] != -1)
            exp_st[cell_id_0] += CS_MAX(imasfl_r[face_id],0) * cfld_yp->val[cell_id_0];
          if (packing_cell[cell_id_1] != -1)
            exp_st[cell_id_1] += CS_MAX(-imasfl_r[face_id],0) * cfld_yp->val[cell_id_1];
        }

        if (f_id == CS_F_(h_l)->id) {
          if (packing_cell[cell_id_0] != -1) {
            exp_st[cell_id_0] += CS_MAX(imasfl_r[face_id],0)
                                 * cfld_yp->val[cell_id_0]
                                 * cs_liq_t_to_h(temp_rain[cell_id_0]);
          }

          if (packing_cell[cell_id_1] != -1) {
            exp_st[cell_id_1] += CS_MAX(imasfl_r[face_id],0)
                                 * cfld_yp->val[cell_id_1]
                                 * cs_liq_t_to_h(temp_rain[cell_id_1]);
          }
        }
      }
    }

  } /* End of test on whether to generate rain */

  /* Source terms for rain drops velocity
   * ==================================== */
  if (ct_opt->solve_rain_velocity) {
    int class_id = 1;
    char vg_lim_name[80];
    sprintf(vg_lim_name, "vg_lim_p_%02d", class_id);

    /* Drops terminal velocity fields */
    cs_real_3_t *vg_lim_p = (cs_real_3_t *)cs_field_by_name(vg_lim_name)->val;
    cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_p");
    cs_real_t *cpro_taup = NULL;
    if (cfld_taup != NULL)
      cpro_taup = cfld_taup->val;

    /* Continuous phase drift velocity */
    cs_real_3_t *vd_c = (cs_real_3_t *)cs_field_by_name("vd_c")->val;

    /* Rain drift velocity variables */
    char f_name[80];
    sprintf(f_name, "v_p_%02d", class_id);
    cs_field_t *f_vp = cs_field_by_name(f_name);

    if (f_id == f_vp->id) {
      cs_real_33_t *_imp_st = (cs_real_33_t *)imp_st;
      cs_real_3_t *_exp_st = (cs_real_3_t *)exp_st;
      cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
      cs_real_3_t *vp = (cs_real_3_t *)f_vp->val;

      for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          /* Explicit source term */
          _exp_st[cell_id][i] += rho_h[cell_id] * cell_f_vol[cell_id]
                              * 1. / cpro_taup[cell_id]
                              * (vel[cell_id][i] + vd_c[cell_id][i]
                                  + vg_lim_p[cell_id][i] - vp[cell_id][i]);
          /* Implicit source term: only diagonal terms */
          _imp_st[cell_id][i][i] += rho_h[cell_id] * cell_f_vol[cell_id]
                                       / cpro_taup[cell_id];
        }
      }
    }
    /* Interfacial pressure drop due to air / rain friction */
    if (f_id == (CS_F_(vel)->id)) {
      /* Casting implicit source term on a 3x3 symmetric matrix */
      cs_real_33_t *_imp_st = (cs_real_33_t *)imp_st;
      cs_real_3_t *_exp_st = (cs_real_3_t *)exp_st;

      /* Rain mass fraction field */
      cs_real_t *y_rain = (cs_real_t *)cfld_yp->val;
      /* Rain drift and velocity fields */
      sprintf(f_name, "vd_p_%02d", class_id);
      cs_real_3_t *cfld_drift
        = (cs_real_3_t *)cs_field_by_name("drift_vel_y_p")->val;
      cs_real_3_t *vp = (cs_real_3_t *)f_vp->val;

      /* Gravity norm */
      cs_real_t g = cs_math_3_norm(cs_glob_physical_constants->gravity);
      for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
        if (y_rain[cell_id] > 0.){
          /* Droplet drift and absolute velocity */
          cs_real_t drift = cs_math_3_norm(cfld_drift[cell_id]);
          cs_real_t v_drop = cs_math_3_norm(vp[cell_id]);

          /* Droplet Reynolds and Eotvos number */
          cs_real_t re_d = rho_h[cell_id] * drift * droplet_diam / visc;
          cs_real_t e_o = g * droplet_diam * (rho_l - rho_h[cell_id]) / sigma;
          /* Sphere drag coefficient */
          cs_real_t cd = 0.;
          if (re_d > 0.){
            cd = (24. / re_d) * (1. + 0.15 * pow(re_d, 0.685));
          }
          /* Droplet terminal velocity */
          cs_real_t v_term = pow((4. * rho_l * droplet_diam * g)
              / (3. * cd * rho_h[cell_id]), 0.5);
          /* Droplet deformation / elongation */
          cs_real_t e_tau = 1. / (1. + 0.148 * pow(e_o, 0.85));
          //FIXME : check positivity of E
          cs_real_t E =   1. - cs_math_pow2(CS_MIN(v_drop / v_term, 1.))
                        * (1. - e_tau);

          /* Total drag coefficient for deformed drop */
          cs_real_t cd_tot = cd * (1. - 0.17185 * (1. - E)
              + 6.692 * cs_math_pow2(1. - E)
              - 6.605 * cs_math_pow3(1. - E));

          /* Air / droplets interfacial area density calculation */
          cs_real_t vol_frac_rain = y_rain[cell_id] * rho_h[cell_id] / rho_l;
          if (vol_frac_rain >= 1.0)
            vol_frac_rain = 1.0;
          cs_real_t a_i =  6.0 * vol_frac_rain * (1.0 - vol_frac_rain)
            / droplet_diam;

          /* Droplet relaxation time */
          cs_real_t tau_d = rho_l * cs_math_pow2(droplet_diam) / (18. * visc);
          /* Final head loss coefficient */
          cs_real_t k_drop = rho_l * (cd_tot * re_d / 24.) * droplet_diam * a_i
            / (6. * tau_d);
          for (int k = 0; k < 3; k++){
            _imp_st[cell_id][k][k] += -cell_f_vol[cell_id] * k_drop;
            _exp_st[cell_id][k] += cell_f_vol[cell_id] * k_drop * vp[cell_id][k];
          }
        }
      }
    }

  } /* End of solve_rain variable check */

  BFT_FREE(packing_cell);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change mass source term from the evaporating liquid to the
 *        bulk, humid air.
 *
 * Careful, this is different from an injection source term, which would
 * normally be handled with a 'cs_equation_add_volume_mass_injection_' function.
 *
 * \param[out]  mass_source     Mass source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_bulk_mass_source_term(cs_real_t         mass_source[])
{
  cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  /* Compute the mass exchange term */
  cs_real_t *imp_st;

  BFT_MALLOC(imp_st, n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    mass_source[cell_id] = 0.0;
    imp_st[cell_id] = 0.0;
  }

  /* Bulk mass source term is stored for pressure */

  cs_ctwr_source_term(CS_F_(p)->id,
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

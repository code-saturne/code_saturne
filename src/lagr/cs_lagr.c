/*============================================================================
 * Methods for particle parameters
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

/*============================================================================
 * Functions dealing with lagrangian initialization
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_math.h"
#include "cs_mesh_location.h"

#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_volume_zone.h"

#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_physical_model.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_default.h"
#include "cs_prototypes.h"

#include "cs_gui_particles.h"
#include "cs_gui_util.h"

#include "cs_lagr.h"

#include "cs_lagr_lec.h"
#include "cs_lagr_geom.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_clogging.h"
#include "cs_lagr_injection.h"
#include "cs_lagr_gradients.h"
#include "cs_lagr_car.h"
#include "cs_lagr_coupling.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_resuspension.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_print.h"
#include "cs_lagr_poisson.h"
#include "cs_lagr_post.h"
#include "cs_lagr_sde.h"
#include "cs_lagr_sde_model.h"
#include "cs_lagr_prototypes.h"

#include "cs_random.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Fixed-size parameters */

cs_lagr_const_dim_t _lagr_const_dim
  = {.nusbrd = 10,
     .ndlaim = 10,
     .ncharm2 = 5,
     .nlayer = 5};

const cs_lagr_const_dim_t *cs_glob_lagr_const_dim
  = (const cs_lagr_const_dim_t *)(&_lagr_const_dim);

/* General dimensions */

cs_lagr_dim_t _lagr_dim = {.ntersl = 0,
                           .n_boundary_stats = 0};

cs_lagr_dim_t *cs_glob_lagr_dim = &_lagr_dim;

/* Time and Lagrangian-Eulerian coupling scheme */

static cs_lagr_time_scheme_t _lagr_time_scheme
  = {.iilagr = 0,
     .isttio = 0,
     .isuila = 0,
     .t_order = 0,
     .modcpl = 0,
     .idirla = 0,
     .idistu = 0,
     .idiffl = 0,
     .ilapoi = 0,
     .iadded_mass = 0,
     .added_mass_const = 0};

/* Main Lagragian physical model parameters */

static cs_lagr_model_t  _lagr_model
  = {.physical_model = 0,
     .n_temperature_layers = 1,
     .deposition = 0,
     .dlvo = 0,
     .roughness = 0,
     .resuspension = 0,
     .clogging = 0,
     .consolidation = 0,
     .precipitation = 0,
     .fouling = 0,
     .n_stat_classes = 0,
     .n_user_variables = 0};

/* particle counter structure and associated pointer */

static cs_lagr_particle_counter_t _lagr_particle_counter
  = {.n_g_cumulative_total = 0,
     .n_g_cumulative_failed = 0,
     .n_g_total = 0,
     .n_g_new = 0,
     .n_g_exit = 0,
     .n_g_deposited = 0,
     .n_g_fouling = 0,
     .n_g_resuspended = 0,
     .n_g_failed = 0,
     .w_total = 0.,
     .w_new = 0.,
     .w_exit = 0.,
     .w_deposited = 0.,
     .w_fouling = 0.,
     .w_resuspended = 0.};

/* lagr specific physics structure and associated pointer */
static cs_lagr_specific_physics_t _cs_glob_lagr_specific_physics
  = {0, 0, 0, 0, 0};
cs_lagr_specific_physics_t *cs_glob_lagr_specific_physics
  = &_cs_glob_lagr_specific_physics;

/* lagr reentrained model structure and associated pointer */
static cs_lagr_reentrained_model_t _cs_glob_lagr_reentrained_model
  = {0, 0, 0, 0, 0, 0, 0};
cs_lagr_reentrained_model_t *cs_glob_lagr_reentrained_model
  = &_cs_glob_lagr_reentrained_model;

/* lagr precipitation model structure and associated pointer */
static cs_lagr_precipitation_model_t _cs_glob_lagr_precipitation_model
  = {0, 0, 0, NULL, NULL, NULL};
cs_lagr_precipitation_model_t *cs_glob_lagr_precipitation_model
  = &_cs_glob_lagr_precipitation_model;

/* lagr clogging model structure and associated pointer */
static cs_lagr_clogging_model_t _cs_glob_lagr_clogging_model = {0, 0, 0, 0};
cs_lagr_clogging_model_t *cs_glob_lagr_clogging_model
   = &_cs_glob_lagr_clogging_model;

/* lagr clogging model structure and associated pointer */
static cs_lagr_consolidation_model_t _cs_glob_lagr_consolidation_model
   = {0, 0, 0, 0};
cs_lagr_consolidation_model_t *cs_glob_lagr_consolidation_model
   = &_cs_glob_lagr_consolidation_model;

/*! current time step status */

static cs_lagr_time_step_t _cs_glob_lagr_time_step
  = {.nor = 0,
     .dtp = 0.,
     .ttclag = 0.};

/* lagr source terms structure and associated pointer */
static cs_lagr_source_terms_t _cs_glob_lagr_source_terms
  = {.ltsdyn = 0,
     .ltsmas = 0,
     .ltsthe = 0,
     .itsli = 0,
     .itske = 0,
     .itste = 0,
     .itsti = 0,
     .itsmas = 0,
     .itsmv1 = NULL,
     .itsmv2 = NULL,
     .itsco = 0,
     .itsfp4 = 0,
     .nstits = 0,
     .npts = 0,
     .ntxerr = 0,
     .vmax = 0,
     .tmamax = 0,
     .st_val = NULL};
cs_lagr_source_terms_t *cs_glob_lagr_source_terms
= &_cs_glob_lagr_source_terms;

/* lagr encrustation structure and associated pointer */
static cs_lagr_encrustation_t _cs_glob_lagr_encrustation
  = {0, 0, NULL, NULL, NULL, NULL, 0} ;
cs_lagr_encrustation_t *cs_glob_lagr_encrustation
  = &_cs_glob_lagr_encrustation;

/* lagr physico chemical structure and associated pointer */
static cs_lagr_physico_chemical_t _cs_glob_lagr_physico_chemical
= {0, 0, 0, 0, 0, 0, 0};
cs_lagr_physico_chemical_t *cs_glob_lagr_physico_chemical
  = &_cs_glob_lagr_physico_chemical;

/* lagr brownian structure and associated pointer */
static cs_lagr_brownian_t _cs_glob_lagr_brownian = {0};
cs_lagr_brownian_t *cs_glob_lagr_brownian = &_cs_glob_lagr_brownian;

/* lagr boundary interactions structure and associated pointer */
static cs_lagr_boundary_interactions_t _cs_glob_lagr_boundary_interactions
  = {.nusbor = 0,
     .npstf = 0,
     .npstft = 0,
     .has_part_impact_nbr = 0,
     .iflmbd = 0,
     .iangbd= 0,
     .ivitbd = 0,
     .iclgst = 0,
     .iencnbbd = 0,
     .iencmabd = 0,
     .iencdibd = 0,
     .iencckbd = 0,
     .inbr = -1,
     .iflm = -1,
     .iang = -1,
     .ivit = -1,
     .ires = -1,
     .iflres = -1,
     .iencnb = -1,
     .iencma = -1,
     .iencdi = -1,
     .iencck = -1,
     .iusb = NULL,
     .imoybr = NULL,
     .inclg = -1,
     .inclgt = -1,
     .iclogt = -1,
     .iscovc = -1,
     .ihdepm = -1,
     .ihdiam = -1,
     .ihsum = -1,
     .tstatp = 0.,
     .nombrd = NULL};

cs_lagr_boundary_interactions_t *cs_glob_lagr_boundary_interactions
  = &_cs_glob_lagr_boundary_interactions;

/* lagr extra modules and associated pointer */

static cs_lagr_extra_module_t _lagr_extra_module
  = {.iturb = 0,
     .itytur = 0,
     .ncharb = 0,
     .ncharm = 0,
     .radiative_model = 0,
     .icp = -1,
     .diftl0 = 0,
     .cmu = 0,
     .visls0 = 0,
     .ustar = NULL,
     .cromf = NULL,
     .pressure = NULL,
     .scal_t = NULL,
     .temperature = NULL,
     .t_gaz = NULL,
     .vel = NULL,
     .viscl = NULL,
     .cpro_viscls = NULL,
     .cpro_cp = NULL,
     .luminance = NULL,
     .x_oxyd = NULL,
     .x_eau = NULL,
     .x_m = NULL,
     .cvar_k = NULL,
     .cvar_ep = NULL,
     .cvar_omg = NULL,
     .cvar_r11 = NULL,
     .cvar_r22 = NULL,
     .cvar_r33 = NULL,
     .cvar_rij = NULL};

cs_lagr_extra_module_t *cs_glob_lagr_extra_module = &_lagr_extra_module;

/* lagr coal combustion structure and associated pointer */

static cs_lagr_coal_comb_t _lagr_coal_comb
  = {0,    0,    0,
     0,    0.0 , 0.0,
     0,    NULL,
     0,    NULL, NULL,
     0,    NULL, NULL, NULL, NULL,
     NULL, NULL, NULL, NULL, NULL,
     NULL, NULL, NULL, NULL, NULL};

/* boundary and volume condition data */

static cs_lagr_zone_data_t  *_boundary_conditions = NULL;
static cs_lagr_zone_data_t  *_volume_conditions = NULL;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointers to global structures; some may be const-qualified in the future
   to set to read-only ouside accessor functions, but this is not done
   at this stage, as the API is not yet stabilized */

cs_lagr_time_scheme_t       *cs_glob_lagr_time_scheme = &_lagr_time_scheme;
cs_lagr_model_t             *cs_glob_lagr_model = &_lagr_model;

const cs_lagr_particle_counter_t  *cs_glob_lagr_particle_counter
                                     = &_lagr_particle_counter;

int cs_glob_lagr_log_frequency_n = 1; /* log every frequency_n time steps */

cs_lagr_time_step_t *cs_glob_lagr_time_step = &_cs_glob_lagr_time_step;

cs_lagr_coal_comb_t *cs_glob_lagr_coal_comb = &_lagr_coal_comb;

cs_real_t *bound_stat = NULL;

const cs_lagr_zone_data_t     *cs_glob_lagr_boundary_conditions = NULL;
const cs_lagr_zone_data_t     *cs_glob_lagr_volume_conditions = NULL;

cs_lagr_internal_condition_t  *cs_glob_lagr_internal_conditions = NULL;

/* Geometry helper arrays */
/*------------------------*/

/*! Unit normals and offsets of boundary faces */
cs_real_4_t  *cs_glob_lagr_b_u_normal = NULL;

/*! Projection matrices for global to local coordinates on boundary faces */

cs_real_33_t  *cs_glob_lagr_b_face_proj = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_lagr_params_pointers(cs_int_t  **p_iilagr,
                          cs_int_t  **p_idepst,
                          cs_int_t  **p_ipreci);

void
cs_f_lagr_dim_pointers(cs_int_t     **p_ntersl);

void
cs_f_lagr_clogging_model_pointers(cs_real_t **jamlim,
                                  cs_real_t **mporos,
                                  cs_real_t **csthpp);

void
cs_f_lagr_consolidation_model_pointers(cs_lnum_t **iconsol,
                                       cs_real_t **rate_consol,
                                       cs_real_t **slope_consol,
                                       cs_real_t **force_consol);

void
cs_f_lagr_source_terms_pointers(cs_int_t **p_ltsdyn,
                                cs_int_t **p_ltsmas,
                                cs_int_t **p_ltsthe,
                                cs_int_t **p_itsli,
                                cs_int_t **p_itske,
                                cs_int_t **p_itste,
                                cs_int_t **p_itsti,
                                cs_int_t **p_itsmas,
                                cs_int_t **p_itsco,
                                cs_int_t **p_itsmv1,
                                cs_int_t **p_itsmv2,
                                cs_int_t  *dim_itsmv1,
                                cs_int_t  *dim_itsmv2);

void
cs_f_lagr_specific_physics(int        *iirayo,
                           int        *ncharb,
                           int        *ncharm,
                           cs_real_t  *diftl0);

void
cs_f_lagr_coal_comb(cs_int_t   *ih2o,
                    cs_int_t   *io2,
                    cs_int_t   *ico,
                    cs_int_t   *iatc,
                    cs_real_t  *prefth,
                    cs_real_t  *trefth,
                    cs_int_t   *natom,
                    cs_real_t  *wmolat,
                    cs_int_t   *ngazem,
                    cs_real_t  *wmole,
                    cs_int_t   *iym1,
                    cs_int_t   *ncharm,
                    cs_real_t  *a1ch,
                    cs_real_t  *h02ch,
                    cs_real_t  *e1ch,
                    cs_real_t  *a2ch,
                    cs_real_t  *e2ch,
                    cs_real_t  *y1ch,
                    cs_real_t  *y2ch,
                    cs_real_t  *cp2ch,
                    cs_real_t  *ahetch,
                    cs_real_t  *ehetch,
                    cs_real_t  *rho0ch,
                    cs_real_t  *xwatch,
                    cs_real_t  *xashch,
                    cs_real_t  *thcdch);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Map some options
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   p_iilagr --> lagrangian model type
 *   p_idepo  --> deposition option flag
 *   p_ipreci --> precipitation option flag
 *----------------------------------------------------------------------------*/

void
cs_f_lagr_params_pointers(cs_int_t  **p_iilagr,
                          cs_int_t  **p_idepst,
                          cs_int_t  **p_ipreci)
{
  *p_iilagr = &_lagr_time_scheme.iilagr;
  *p_idepst = &_lagr_model.deposition;
  *p_ipreci  = &_lagr_model.precipitation;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global lagr dim structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ntersl  --> pointer to cs_glob_lagr_dim->ntersl
 *----------------------------------------------------------------------------*/

void
cs_f_lagr_dim_pointers(cs_int_t   **p_ntersl)
{
  *p_ntersl = &(_lagr_dim.ntersl);
}

void
cs_f_lagr_clogging_model_pointers(cs_real_t **jamlim,
                                  cs_real_t **mporos,
                                  cs_real_t **csthpp)
{
  *jamlim = &cs_glob_lagr_clogging_model->jamlim;
  *mporos = &cs_glob_lagr_clogging_model->mporos;
  *csthpp = &cs_glob_lagr_clogging_model->csthpp;
}

void
cs_f_lagr_consolidation_model_pointers(cs_lnum_t **iconsol,
                                       cs_real_t **rate_consol,
                                       cs_real_t **slope_consol,
                                       cs_real_t **force_consol)
{
  *iconsol      = &cs_glob_lagr_consolidation_model->iconsol;
  *rate_consol  = &cs_glob_lagr_consolidation_model->rate_consol;
  *slope_consol = &cs_glob_lagr_consolidation_model->slope_consol;
  *force_consol = &cs_glob_lagr_consolidation_model->force_consol;
}

void
cs_f_lagr_source_terms_pointers(cs_int_t **p_ltsdyn,
                                cs_int_t **p_ltsmas,
                                cs_int_t **p_ltsthe,
                                cs_int_t **p_itsli,
                                cs_int_t **p_itske,
                                cs_int_t **p_itste,
                                cs_int_t **p_itsti,
                                cs_int_t **p_itsmas,
                                cs_int_t **p_itsco,
                                cs_int_t **p_itsmv1,
                                cs_int_t **p_itsmv2,
                                cs_int_t  *dim_itsmv1,
                                cs_int_t  *dim_itsmv2)
{
  *p_ltsdyn = &cs_glob_lagr_source_terms->ltsdyn;
  *p_ltsmas = &cs_glob_lagr_source_terms->ltsmas;
  *p_ltsthe = &cs_glob_lagr_source_terms->ltsthe;
  *p_itsli  = &cs_glob_lagr_source_terms->itsli;
  *p_itske  = &cs_glob_lagr_source_terms->itske;
  *p_itste  = &cs_glob_lagr_source_terms->itste;
  *p_itsti  = &cs_glob_lagr_source_terms->itsti;
  *p_itsmas = &cs_glob_lagr_source_terms->itsmas;
  *p_itsco  = &cs_glob_lagr_source_terms->itsco;

  if (cs_glob_lagr_source_terms->itsmv1 == NULL)
    BFT_MALLOC(cs_glob_lagr_source_terms->itsmv1,
               cs_glob_lagr_const_dim->ncharm2, int);
  *p_itsmv1 = cs_glob_lagr_source_terms->itsmv1;
  *dim_itsmv1 = cs_glob_lagr_const_dim->ncharm2;

  if (cs_glob_lagr_source_terms->itsmv2 == NULL)
    BFT_MALLOC(cs_glob_lagr_source_terms->itsmv2,
               cs_glob_lagr_const_dim->ncharm2, int);
  *p_itsmv2 = cs_glob_lagr_source_terms->itsmv2;
  *dim_itsmv2 = cs_glob_lagr_const_dim->ncharm2;
}

void
cs_f_lagr_specific_physics(int        *iirayo,
                           int        *ncharb,
                           int        *ncharm,
                           cs_real_t  *diftl0)
{
  _lagr_extra_module.iturb  = cs_glob_turb_model->iturb;
  _lagr_extra_module.itytur = cs_glob_turb_model->itytur;
  _lagr_extra_module.ncharb = *ncharb;
  _lagr_extra_module.ncharm = *ncharm;
  _lagr_extra_module.icp    = cs_glob_fluid_properties->icp;

  _lagr_extra_module.radiative_model = *iirayo;
  _lagr_extra_module.diftl0 = *diftl0;
  _lagr_extra_module.cmu    = cs_turb_cmu;
}

void
cs_f_lagr_coal_comb(cs_int_t   *ih2o,
                    cs_int_t   *io2,
                    cs_int_t   *ico,
                    cs_int_t   *iatc,
                    cs_real_t  *prefth,
                    cs_real_t  *trefth,
                    cs_int_t   *natom,
                    cs_real_t  *wmolat,
                    cs_int_t   *ngazem,
                    cs_real_t  *wmole,
                    cs_int_t   *iym1,
                    cs_int_t   *ncharm,
                    cs_real_t  *a1ch,
                    cs_real_t  *h02ch,
                    cs_real_t  *e1ch,
                    cs_real_t  *a2ch,
                    cs_real_t  *e2ch,
                    cs_real_t  *y1ch,
                    cs_real_t  *y2ch,
                    cs_real_t  *cp2ch,
                    cs_real_t  *ahetch,
                    cs_real_t  *ehetch,
                    cs_real_t  *rho0ch,
                    cs_real_t  *xwatch,
                    cs_real_t  *xashch,
                    cs_real_t  *thcdch)
{
  cs_glob_lagr_coal_comb->ih2o   = *ih2o;
  cs_glob_lagr_coal_comb->io2    = *io2;
  cs_glob_lagr_coal_comb->ico    = *ico;

  cs_glob_lagr_coal_comb->iatc   = *iatc;
  cs_glob_lagr_coal_comb->prefth = *prefth;
  cs_glob_lagr_coal_comb->trefth = *trefth;

  cs_glob_lagr_coal_comb->natom = *natom;
  cs_glob_lagr_coal_comb->wmolat = wmolat;

  cs_glob_lagr_coal_comb->ngazem = *ngazem;
  cs_glob_lagr_coal_comb->wmole  = wmole;
  cs_glob_lagr_coal_comb->iym1   = iym1;

  cs_glob_lagr_coal_comb->ncharm = *ncharm;
  cs_glob_lagr_coal_comb->a1ch   = a1ch;
  cs_glob_lagr_coal_comb->h02ch  = h02ch;
  cs_glob_lagr_coal_comb->e1ch   = e1ch;
  cs_glob_lagr_coal_comb->a2ch   = a2ch;
  cs_glob_lagr_coal_comb->e2ch   = e2ch;
  cs_glob_lagr_coal_comb->y1ch   = y1ch;
  cs_glob_lagr_coal_comb->y2ch   = y2ch;
  cs_glob_lagr_coal_comb->cp2ch  = cp2ch;
  cs_glob_lagr_coal_comb->ahetch = ahetch;
  cs_glob_lagr_coal_comb->ehetch = ehetch;
  cs_glob_lagr_coal_comb->rho0ch = rho0ch;
  cs_glob_lagr_coal_comb->xwatch = xwatch;
  cs_glob_lagr_coal_comb->xashch = xashch;
  cs_glob_lagr_coal_comb->thcdch = thcdch;

}

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static void
_lagr_map_fields_default(void)
{
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0) {
    _lagr_extra_module.cromf       = cs_field_by_name_try("rho_gas");
  }
  else {
    _lagr_extra_module.cromf       = cs_field_by_name_try("density");
  }

  _lagr_extra_module.pressure    = cs_field_by_name_try("pressure");

  _lagr_extra_module.luminance   = cs_field_by_name_try("luminance");
  if (cs_field_by_name_try("velocity_1") != NULL) {
    /* we are probably using NEPTUNE_CFD */
    _lagr_extra_module.vel         = cs_field_by_name_try("velocity_1");

    _lagr_extra_module.cvar_k      = cs_field_by_name_try("TurbKineEner_k_1");
    _lagr_extra_module.cvar_ep     = cs_field_by_name_try("epsilon_1");
    _lagr_extra_module.cvar_omg    = NULL;
    _lagr_extra_module.cvar_r11    = cs_field_by_name_try("lagr_r11");
    _lagr_extra_module.cvar_r22    = cs_field_by_name_try("lagr_r22");
    _lagr_extra_module.cvar_r33    = cs_field_by_name_try("lagr_r33");
    _lagr_extra_module.cvar_rij    = cs_field_by_name_try("lagr_rij");
    _lagr_extra_module.viscl       = cs_field_by_name_try
                                       ("molecular_viscosity_1");
    _lagr_extra_module.scal_t      = cs_field_by_name_try("enthalpy_1");
    _lagr_extra_module.cpro_viscls = cs_field_by_name_try
                                       ("thermal_conductivity_1");
    _lagr_extra_module.cpro_cp     = cs_field_by_name_try("specific_heat_1");
    _lagr_extra_module.temperature = cs_field_by_name_try("lagr_temperature");
    _lagr_extra_module.t_gaz       = NULL;
    _lagr_extra_module.x_oxyd      = NULL;
    _lagr_extra_module.x_eau       = NULL;
    _lagr_extra_module.x_m         = NULL;
    _lagr_extra_module.cromf       = cs_field_by_name_try("density_1");
    /* TODO FIXME */
    _lagr_extra_module.visls0      = 0.;

    _lagr_extra_module.ustar   = cs_field_by_name_try("wall_friction_velocity");
  }
  else {
    /* we use Code_Saturne */
    _lagr_extra_module.vel         = cs_field_by_name_try("velocity");
    _lagr_extra_module.cvar_k      = cs_field_by_name_try("k");
    _lagr_extra_module.cvar_ep     = cs_field_by_name_try("epsilon");
    _lagr_extra_module.cvar_omg    = cs_field_by_name_try("omega");
    _lagr_extra_module.cvar_r11    = cs_field_by_name_try("r11");
    _lagr_extra_module.cvar_r22    = cs_field_by_name_try("r22");
    _lagr_extra_module.cvar_r33    = cs_field_by_name_try("r33");
    _lagr_extra_module.cvar_rij    = cs_field_by_name_try("rij");
    _lagr_extra_module.viscl       = cs_field_by_name_try("molecular_viscosity");
    _lagr_extra_module.cpro_viscls = NULL;

    if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE)
        _lagr_extra_module.scal_t    = cs_field_by_name_try("temperature");
    else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY)
        _lagr_extra_module.scal_t    = cs_field_by_name_try("enthalpy");
    else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TOTAL_ENERGY)
        _lagr_extra_module.scal_t    = cs_field_by_name_try("total_energy");
    else
        _lagr_extra_module.scal_t    = NULL;

    if (_lagr_extra_module.scal_t != NULL) {
        _lagr_extra_module.visls0
            = cs_field_get_key_double(_lagr_extra_module.scal_t,
                    cs_field_key_id("diffusivity_ref"));

        int l_id = cs_field_get_key_int(_lagr_extra_module.scal_t,
                cs_field_key_id("diffusivity_id"));
        if (l_id >= 0)
            _lagr_extra_module.cpro_viscls = cs_field_by_id(l_id);
    }

    _lagr_extra_module.cpro_cp     = cs_field_by_name_try("specific_heat");
    _lagr_extra_module.temperature = cs_field_by_name_try("temperature");
    _lagr_extra_module.t_gaz       = cs_field_by_name_try("t_gas");
    _lagr_extra_module.x_oxyd      = cs_field_by_name_try("ym_o2");
    _lagr_extra_module.x_eau       = cs_field_by_name_try("ym_h2o");
    _lagr_extra_module.x_m         = cs_field_by_name_try("xm");

    _lagr_extra_module.ustar  = cs_field_by_name_try("ustar");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reallocate a zone injection data structure for a given zone and set.
 *
 * The data is reallocated (expanded) if needed, and section relative
 * to the given class id is set to default values.
 *
 * Note that the structure is allocated in one block, large enough to contain
 * the base structure and model-dependent arrays it may point to.
 *
 * \param[in]        location_id       id of associated location
 * \param[in]        zone_id           id of associated zone
 * \param[in]        set_id            id of requested class
 * \param[in, out]   n_injection_sets  number of injection sets for this zone
 * \param[in, out]   injection_sets    zone injection data for this zone
 */
/*----------------------------------------------------------------------------*/

static void
_zone_injection_set_init(int                        location_id,
                         int                        zone_id,
                         int                        set_id,
                         int                       *n_injection_sets,
                         cs_lagr_injection_set_t  **injection_sets)
{
  /* Simply reset if allocation is large enough */

  if (set_id < *n_injection_sets) {
    cs_lagr_injection_set_t *zis
      = (cs_lagr_injection_set_t *)((*injection_sets) + set_id);
    cs_lagr_injection_set_default(zis);
  }

  /* Reallocate memory so as to maintain
     sub-arrays at the end of the structure */

  else {

    cs_lagr_injection_set_t *_zis  = *injection_sets;
    BFT_REALLOC(_zis, set_id+1, cs_lagr_injection_set_t);

    for (int i = *n_injection_sets; i <= set_id; i++) {

      cs_lagr_injection_set_t *zis = (cs_lagr_injection_set_t *)(_zis + i);

      memset(zis, 0, sizeof(cs_lagr_injection_set_t));

      zis->location_id = location_id;
      zis->zone_id = zone_id;
      zis->set_id = set_id;

      cs_lagr_injection_set_default(zis);

    }

    *n_injection_sets = set_id+1;
    *injection_sets = _zis;

  }
}

/*----------------------------------------------------------------------------
 * Initialize or update a zone data structure
 *
 * parameters:
 *   zone_data   <-> pointer to zone data structure pointer
 *   location_id <-- mesh location id
 *   n_zones     <-- number of zones
 *
 * returns:
 *   a new defined cs_lagr_injection_sets_t structure
 *----------------------------------------------------------------------------*/

static void
_update_zone_data_struct(cs_lagr_zone_data_t  **zone_data,
                         int                    location_id,
                         int                    n_zones)
{
  cs_lagr_zone_data_t  *zd = *zone_data;

  if (*zone_data == NULL) {
    BFT_MALLOC(zd, 1, cs_lagr_zone_data_t);
    zd->location_id = location_id;
    zd->n_zones = 0;
    zd->zone_type = NULL;
    zd->n_injection_sets = NULL;
    zd->injection_set = NULL;
    zd->elt_type = NULL;
    zd->particle_flow_rate = NULL;
    *zone_data = zd;
  }

  if (zd->n_zones < n_zones) {
    int n_stats = cs_glob_lagr_model->n_stat_classes + 1;
    BFT_REALLOC(zd->zone_type, n_zones, int);
    BFT_REALLOC(zd->n_injection_sets, n_zones, int);
    BFT_REALLOC(zd->injection_set, n_zones, cs_lagr_injection_set_t *);
    BFT_REALLOC(zd->particle_flow_rate, n_zones*n_stats, cs_real_t);
    for (int i = zd->n_zones; i < n_zones; i++) {
      zd->zone_type[i] = -1;
      zd->n_injection_sets[i] = 0;
      zd->injection_set[i] = NULL;
    }
    for (int i = zd->n_zones*n_stats; i < n_zones*n_stats; i++)
      zd->particle_flow_rate[i] = 0;

    zd->n_zones = n_zones;
  }
}

/*----------------------------------------------------------------------------
 * Initialize a cs_lagr_internal_condition_t structure.
 *
 * returns:
 *   a new defined cs_lagr_internal_condition_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_internal_condition_t *
_create_internal_cond_struct(void)
{
  cs_lagr_internal_condition_t *internal_cond = NULL;
  cs_mesh_t *mesh = cs_glob_mesh;

  BFT_MALLOC(internal_cond, 1, cs_lagr_internal_condition_t);

  BFT_MALLOC(internal_cond->i_face_zone_id, mesh->n_i_faces, int);

  for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++)
    internal_cond->i_face_zone_id[i] = -1;

  return internal_cond;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update (or build) boundary face types.
 */
/*----------------------------------------------------------------------------*/

static void
_update_boundary_face_type(void)
{
  cs_lagr_zone_data_t *bcs = cs_lagr_get_boundary_conditions();

  const cs_mesh_t *mesh = cs_glob_mesh;

  BFT_REALLOC(bcs->elt_type, mesh->n_b_faces, char);

  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++)
    bcs->elt_type[i] = 0;

  for (int z_id = 0; z_id < bcs->n_zones; z_id++) {

    if (bcs->zone_type[z_id] < 0) /* ignore undefined zones */
      continue;

    char z_type = bcs->zone_type[z_id];

    const cs_zone_t  *z = cs_boundary_zone_by_id(z_id);
    for (cs_lnum_t i = 0; i < z->n_elts; i++)
      bcs->elt_type[z->elt_ids[i]] = z_type;

  }

  /* Check definitions */

  {
    int *bc_flag;
    BFT_MALLOC(bc_flag, mesh->n_b_faces, int);

    for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++)
      bc_flag[i] = bcs->elt_type[i];

    cs_boundary_conditions_error(bc_flag, _("Lagrangian boundary type"));

    BFT_FREE(bc_flag);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointers to lagrangian arrays
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   dim_bound_stat   --> dimensions for bound_stat pointer
 *   p_bound_stat     --> bound_stat pointer
 *----------------------------------------------------------------------------*/

void
cs_lagr_init_c_arrays(int          dim_cs_glob_lagr_source_terms[2],
                      cs_real_t  **p_cs_glob_lagr_source_terms)
{
  cs_lnum_t  n_b_faces = cs_glob_mesh->n_b_faces;
  int   n_boundary_stats = cs_glob_lagr_dim->n_boundary_stats;

  assert(bound_stat == NULL);

  if (n_boundary_stats > 0)
    BFT_MALLOC(bound_stat, n_b_faces * n_boundary_stats, cs_real_t);

  BFT_MALLOC(cs_glob_lagr_source_terms->st_val,
             cs_glob_lagr_dim->ntersl * cs_glob_mesh->n_cells_with_ghosts,
             cs_real_t);

  *p_cs_glob_lagr_source_terms     = cs_glob_lagr_source_terms->st_val;
  dim_cs_glob_lagr_source_terms[0] = cs_glob_mesh->n_cells_with_ghosts;
  dim_cs_glob_lagr_source_terms[1] = cs_glob_lagr_dim->ntersl;
}

/*----------------------------------------------------------------------------
 * Free lagrangian arrays
 *
 * This function is intended for use by Fortran wrappers.
 *----------------------------------------------------------------------------*/

void
cs_lagr_finalize(void)
{
  int  n_boundary_stats = cs_glob_lagr_dim->n_boundary_stats;

  if (n_boundary_stats > 0) {
    assert(bound_stat != NULL);
    BFT_FREE(bound_stat);
  }

  BFT_FREE(cs_glob_lagr_precipitation_model->nbprec);
  BFT_FREE(cs_glob_lagr_precipitation_model->solub);

  BFT_FREE(cs_glob_lagr_precipitation_model->mp_diss);

  BFT_FREE(cs_glob_lagr_source_terms->st_val);

  /* geometry */

  BFT_FREE(cs_glob_lagr_b_u_normal);
  BFT_FREE(cs_glob_lagr_b_face_proj);

  /* encrustation pointers */

  BFT_FREE(cs_glob_lagr_encrustation->enc1);
  BFT_FREE(cs_glob_lagr_encrustation->enc2);
  BFT_FREE(cs_glob_lagr_encrustation->tprenc);
  BFT_FREE(cs_glob_lagr_encrustation->visref);

  /* boundary interaction pointers */

  BFT_FREE(cs_glob_lagr_boundary_interactions->iusb);
  BFT_FREE(cs_glob_lagr_boundary_interactions->imoybr);

  for (int i = 0; i < cs_glob_lagr_dim->n_boundary_stats; i++)
    BFT_FREE(cs_glob_lagr_boundary_interactions->nombrd[i]);
  BFT_FREE(cs_glob_lagr_boundary_interactions->nombrd);

  /* Source terms */

  BFT_FREE(cs_glob_lagr_source_terms->itsmv1);
  BFT_FREE(cs_glob_lagr_source_terms->itsmv2);

  /* Statistics */

  cs_lagr_stat_finalize();

  /* Also close log file (TODO move this) */

  cs_lagr_print_finalize();

  /* Close tracking structures */

  cs_lagr_tracking_finalize();

  cs_lagr_finalize_zone_conditions();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to injection set structure.
 *
 * This access method ensures the strucure is initialized for the given
 * zone and injection set.
 *
 * \param[in]  zone_data  pointer to boundary or volume conditions structure
 * \param[in]  zone_id    zone id
 * \param[in]  set_id     injection set id
 *
 * \return pointer to injection set data structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_injection_set_t *
cs_lagr_get_injection_set(cs_lagr_zone_data_t  *zone_data,
                          int                   zone_id,
                          int                   set_id)
{
  assert(zone_data != NULL);
  assert(zone_id >= 0 && set_id >= 0);
  assert(zone_id < zone_data->n_zones);

  if (set_id >= zone_data->n_injection_sets[zone_id])
    _zone_injection_set_init(zone_data->location_id,
                             zone_id,
                             set_id,
                             &(zone_data->n_injection_sets[zone_id]),
                             &(zone_data->injection_set[zone_id]));

  return &(zone_data->injection_set[zone_id][set_id]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize injection set data structure fields to defaults.
 *
 * \param[in, out]   zis  pointer to structure to initialize
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_injection_set_default(cs_lagr_injection_set_t  *zis)
{
  assert(zis != NULL);

  zis->n_inject             =  0;
  zis->injection_frequency  =  0;

  zis->injection_profile_func = NULL;
  zis->injection_profile_input = NULL;

  zis->velocity_profile     = -2;
  zis->temperature_profile  = -2;

  if (cs_glob_lagr_model->physical_model == 2)
    zis->coal_number        = -2;

  zis->cluster              =  0;

  zis->velocity_magnitude   = - cs_math_big_r;

  for (int  i = 0; i < 3; i++)
    zis->velocity[i]        = - cs_math_big_r;

  zis->stat_weight          = - cs_math_big_r;
  zis->diameter             = - cs_math_big_r;
  zis->diameter_variance    = - cs_math_big_r;
  zis->density              = - cs_math_big_r;

  zis->temperature      = - cs_math_big_r;
  zis->cp               = - cs_math_big_r;
  zis->emissivity       = - cs_math_big_r;

  zis->flow_rate     = 0.0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get read/write pointer to global particle counter
 *
 * \return
 *   pointer to lagrangian particle counter structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_counter_t *
cs_lagr_get_particle_counter(void)
{
  return &_lagr_particle_counter;
}

/*----------------------------------------------------------------------------*/
/*!
  \brief Update global particle counter
 *
 * All fields handled in the local particle set are updated relative
 * to that data (using global sums).
 *
 * \return  pointer to lagrangian particle counter structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_particle_counter_t *
cs_lagr_update_particle_counter(void)
{
  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;
  cs_lagr_particle_counter_t *pc = &_lagr_particle_counter;

  cs_gnum_t gcount[7] = {p_set->n_particles,
                         p_set->n_part_new,
                         p_set->n_part_out,
                         p_set->n_part_dep,
                         p_set->n_part_fou,
                         p_set->n_part_resusp,
                         p_set->n_failed_part};

  cs_real_t wsum[6] = {p_set->weight,
                       p_set->weight_new,
                       p_set->weight_out,
                       p_set->weight_dep,
                       p_set->weight_fou,
                       p_set->weight_resusp};

  cs_parall_counter(gcount, 7);
  cs_parall_sum(6, CS_REAL_TYPE, wsum);

  pc->n_g_total = gcount[0];
  pc->n_g_new = gcount[1];
  pc->n_g_exit = gcount[2];
  pc->n_g_deposited = gcount[3];
  pc->n_g_fouling = gcount[4];
  pc->n_g_resuspended = gcount[5];
  pc->n_g_failed = gcount[6];

  pc->w_total = wsum[0];
  pc->w_new = wsum[1];
  pc->w_exit = wsum[2];
  pc->w_deposited = wsum[3];
  pc->w_fouling = wsum[4];
  pc->w_resuspended = wsum[5];

  return pc;
}

/*----------------------------------------------------------------------------*/
/*! \brief Provide access to cs_lagr_specific_physics_t
 *
 * \return
 *   pointer to lagrangian specific physics options
 *----------------------------------------------------------------------------*/

cs_lagr_specific_physics_t *
cs_get_lagr_specific_physics(void)
{
  return &_cs_glob_lagr_specific_physics;
}

/*----------------------------------------------------------------------------*/
/*! \brief Provide access to cs_lagr_reentrained_model_t
 *
 * \return
 *   pointer to lagrangian reentrained model options
 *----------------------------------------------------------------------------*/

cs_lagr_reentrained_model_t *
cs_get_lagr_reentrained_model(void)
{
  return &_cs_glob_lagr_reentrained_model;
}

/*----------------------------------------------------------------------------*/
/*! \brief Provide access to cs_lagr_precipitation_model_t
 *
 * \return
 *   pointer to lagrangian precipitation model options
 *----------------------------------------------------------------------------*/

cs_lagr_precipitation_model_t *
cs_get_lagr_precipitation_model(void)
{
  return &_cs_glob_lagr_precipitation_model;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_clogging_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_clogging_model_t *
cs_get_lagr_clogging_model(void)
{
  return &_cs_glob_lagr_clogging_model;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_consolidation_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_consolidation_model_t *
cs_get_lagr_consolidation_model(void)
{
  return &_cs_glob_lagr_consolidation_model;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_time_step_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_time_step_t *
cs_get_lagr_time_step(void)
{
  return &_cs_glob_lagr_time_step;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_source_terms_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_source_terms_t *
cs_get_lagr_source_terms(void)
{
  return &_cs_glob_lagr_source_terms;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_encrustation_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_encrustation_t *
cs_get_lagr_encrustation(void)
{
  return &_cs_glob_lagr_encrustation;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_physico_chemical_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_physico_chemical_t *
cs_get_lagr_physico_chemical(void)
{
  return &_cs_glob_lagr_physico_chemical;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_brownian_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_brownian_t *
cs_get_lagr_brownian(void)
{
  return &_cs_glob_lagr_brownian;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main internal conditions structure.
 *
 * \return
 *   pointer to current internal_contditions or NULL
 */
/*----------------------------------------------------------------------------*/

cs_lagr_internal_condition_t  *
cs_lagr_get_internal_conditions(void)
{
  /* Define a structure with default parameters if not done yet */

  if (cs_glob_lagr_internal_conditions == NULL)
    cs_glob_lagr_internal_conditions = _create_internal_cond_struct();

  if (cs_glob_lagr_internal_conditions->i_face_zone_id == NULL) {
    BFT_MALLOC(cs_glob_lagr_internal_conditions->i_face_zone_id,
               cs_glob_mesh->n_i_faces,
               int);

    for (cs_lnum_t i = 0; i < cs_glob_mesh->n_i_faces; i++)
      cs_glob_lagr_internal_conditions->i_face_zone_id[i] = -1;

  }

  return cs_glob_lagr_internal_conditions;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main boundary conditions structure.
 *
 * \return  pointer to current boundary zone data structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_zone_data_t  *
cs_lagr_get_boundary_conditions(void)
{
  /* Define a structure with default parameters if not done yet */

  cs_lnum_t n_zones = cs_boundary_zone_n_zones();

  _update_zone_data_struct(&_boundary_conditions,
                           CS_MESH_LOCATION_BOUNDARY_FACES,
                           n_zones);

  cs_glob_lagr_boundary_conditions = _boundary_conditions;

  return _boundary_conditions;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main volume conditions structure.
 *
 * \return pointer to current volume zone data structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_zone_data_t  *
cs_lagr_get_volume_conditions(void)
{
  /* Define a structure with default parameters if not done yet */

  cs_lnum_t n_zones = cs_volume_zone_n_zones();

  _update_zone_data_struct(&_volume_conditions,
                           CS_MESH_LOCATION_CELLS,
                           n_zones);

  cs_glob_lagr_volume_conditions = _volume_conditions;

  return _volume_conditions;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main internal conditions structure.
 *
 * \return pointer to current internal conditions structure
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_finalize_zone_conditions(void)
{
  cs_lagr_zone_data_t  *zda[2] = {_boundary_conditions,
                                  _volume_conditions};

  for (int i = 0; i < 2; i++) {

    cs_lagr_zone_data_t  *zd = zda[i];

    if (zd != NULL) {

      BFT_FREE(zd->zone_type);
      for (int j = 0; j < zd->n_zones; j++)
        BFT_FREE(zd->injection_set[j]);
      BFT_FREE(zd->injection_set);
      BFT_FREE(zd->n_injection_sets);

      BFT_FREE(zd->elt_type);
      BFT_FREE(zd->particle_flow_rate);

      BFT_FREE(zda[i]);

    }

  }
}

/*----------------------------------------------------------------------------
 * Destroy finalize the global cs_lagr_internal_condition_t structure.
 *----------------------------------------------------------------------------*/

void cs_lagr_finalize_internal_cond(void)
{
  cs_lagr_internal_condition_t  *internal_cond
    = cs_glob_lagr_internal_conditions;
  if (internal_cond != NULL) {
    BFT_FREE(internal_cond->i_face_zone_id);
    BFT_FREE(internal_cond);
  }
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_boundary_interactions_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_boundary_interactions_t *
cs_get_lagr_boundary_interactions(void)
{
  return &_cs_glob_lagr_boundary_interactions;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_extra_module_t
 *
 *----------------------------------------------------------------------------*/

cs_lagr_extra_module_t *
cs_get_lagr_extra_module(void)
{
  return &_lagr_extra_module;
}

/*----------------------------------------------------------------------------
 * Prepare for execution of the Lagrangian model.
 *
 * This should be called before the fist call to cs_lagr_solve_time_step.
 *
 *  parameters:
 *    dt     <-- time step (per cell)
 *----------------------------------------------------------------------------*/

void
cs_lagr_solve_initialize(const cs_real_t  *dt)
{
  CS_UNUSED(dt);

  /* For frozen field:
     values at previous time step = values at current time step */

  if (cs_glob_lagr_time_scheme->iilagr == 3) {

    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++){

      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE)
        cs_field_current_to_previous(f);

    }

  }

  /* First call (initializations)
     ---------------------------- */

  _lagr_map_fields_default();

  cs_lagr_tracking_initialize();

  /* first initializations */

  cs_lagr_post_init();

  /* Read particle restart data */

  if (cs_glob_lagr_time_scheme->iilagr > 0)
    cs_lagr_restart_read_p();
}

/*----------------------------------------------------------------------------
 * Execute one time step of the Lagrangian model.
 *
 * This is the main function for that model.
 *
 *  parameters:
 *    itypfb <-- boundary face types
 *    dt     <-- time step (per cell)
 *----------------------------------------------------------------------------*/

void
cs_lagr_solve_time_step(const int         itypfb[],
                        const cs_real_t  *dt)
{
  static int ipass = 0;
  cs_time_step_t *ts = cs_get_glob_time_step();

  int  mode;

  cs_lagr_boundary_interactions_t *lag_bdi = cs_glob_lagr_boundary_interactions;
  cs_lagr_model_t *lagr_model = cs_glob_lagr_model;
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_particle_counter_t *part_c = cs_lagr_get_particle_counter();

  cs_lnum_t ncelet = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_lnum_t *ifabor = cs_glob_mesh->b_face_cells;
  cs_real_t *surfbo = cs_glob_mesh_quantities->b_face_surf;
  cs_real_t *surfbn = cs_glob_mesh_quantities->b_face_normal;

  /* Allocate temporary arrays */
  cs_real_3_t *gradpr;
  BFT_MALLOC(gradpr, ncelet, cs_real_3_t);
  cs_real_t *w1, *w2;
  BFT_MALLOC(w1, ncelet, cs_real_t);
  BFT_MALLOC(w2, ncelet, cs_real_t);

  /* Allocate other arrays depending on user options    */
  cs_real_33_t *gradvf = NULL;
  if (cs_glob_lagr_time_scheme->modcpl > 0)
    BFT_MALLOC(gradvf, ncelet, cs_real_33_t);

  cs_real_t *tempp = NULL;
  if (   lagr_model->clogging == 1
      || lagr_model->roughness == 1
      || lagr_model->dlvo == 1)
    BFT_MALLOC(tempp, cs_glob_mesh->n_cells, cs_real_t);

  ipass ++;

  static cs_real_t *vislen = NULL;
  if ((lagr_model->deposition == 1) && (ipass == 1)) {
    BFT_MALLOC(vislen, n_b_faces, cs_real_t);
    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
      vislen[ifac] = -cs_math_big_r;
  }

  /* Initialization
     -------------- */

  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;

  /* First call (initializations)
     ---------------------------- */

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1) {

    /* If the deposition model is activated */

    if (lagr_model->deposition >= 1) {

      cs_real_t ustarmoy = 0.0;
      cs_real_t surftot  = 0.0;
      cs_real_3_t  dtmp;
      dtmp[0]  = 0.0;

      /* Boundary faces data  */

      cs_lagr_geom();

      /* Average friction velocity calculation    */

      for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

        if (itypfb[ifac] == CS_SMOOTHWALL || itypfb[ifac] == CS_ROUGHWALL) {

          cs_lnum_t iel   = ifabor[ifac];

          /* the density pointer according to the flow location */
          cs_real_t romf   = extra->cromf->val[iel];
          cs_real_t visccf = extra->viscl->val[iel] / romf;

          cs_real_t _ustar = CS_MAX(extra->ustar->val[ifac], 1e-15);

          cs_real_t surfb = surfbo[ifac];
          ustarmoy     = ustarmoy + surfb * _ustar;
          surftot      = surftot + surfb;
          vislen[ifac] = visccf / _ustar; // FIXME to be coherent with wall fn: y/y+
          dtmp[0]      = dtmp[0] + 1.0;

        }

      }

      if (cs_glob_rank_id >= 0) {
        dtmp[1] = ustarmoy;
        dtmp[2] = surftot;

        cs_parall_sum(3, CS_DOUBLE, dtmp);

        ustarmoy = dtmp[1];
        surftot  = dtmp[2];

      }

      if (dtmp[0] > 0.0)
        ustarmoy = ustarmoy / surftot;

      /*  Average friction velocity display  */

      if (cs_glob_rank_id <= 0)
        bft_printf(_("   ** LAGRANGIAN MODULE:\n"
                     "   ** deposition submodel\n"
                     "---------------------------------------------\n"
                     "** Mean friction velocity  (ustar) = %7.3f\n"
                     "---------------------------------------------\n"),
                   ustarmoy);

    }

  }

  /* Update boundary condition types;
     in most cases, this should be useful only at the first iteration,
     but we prefer to be safe in case of advanced uses with time-dependent
     zones or zone types */

  /* User initialization by set and boundary */

  if (cs_gui_file_is_loaded() && ts->nt_cur == ts->nt_prev +1)
    cs_gui_particles_bcs();

  /* User initialization by set and zone */

  cs_user_lagr_boundary_conditions(itypfb);
  cs_user_lagr_volume_conditions();

  _update_boundary_face_type();

  /* Initialize counter */

  part_c->n_g_total = 0;
  part_c->n_g_new = 0;
  part_c->n_g_exit = 0;
  part_c->n_g_deposited = 0;
  part_c->n_g_fouling = 0;
  part_c->n_g_resuspended = 0;
  part_c->n_g_failed = 0;
  part_c->w_total = 0;
  part_c->w_new = 0;
  part_c->w_exit = 0;
  part_c->w_deposited = 0;
  part_c->w_fouling = 0;
  part_c->w_resuspended = 0;

  /* Initialization for the dlvo, roughness and clogging  model
     ---------------------------------------------------------- */

  if (   lagr_model->dlvo == 1
      || lagr_model->roughness == 1
      || lagr_model->clogging == 1) {

    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {

      if (extra->scal_t != NULL) {

        if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
            && cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
          tempp[iel] = extra->scal_t->val[iel]
            + cs_physical_constants_celsius_to_kelvin;

        else if (   cs_glob_thermal_model->itherm ==
                            CS_THERMAL_MODEL_TEMPERATURE
                 && cs_glob_thermal_model->itpscl ==
                            CS_TEMPERATURE_SCALE_KELVIN)
          tempp[iel] = extra->scal_t->val[iel];

        else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

          mode = 1;
          CS_PROCF(usthht,USTHHT)(&mode, &extra->scal_t->val[iel], &tempp[iel]);

        }

      }

      else
        tempp[iel] = cs_glob_fluid_properties->t0;

    }

  }

  /* Initialization for the dlvo model
     --------------------------------- */

  if (lagr_model->dlvo == 1)
    cs_lagr_dlvo_init(cs_glob_lagr_physico_chemical->epseau,
                      cs_glob_lagr_physico_chemical->fion,
                      tempp,
                      cs_glob_lagr_physico_chemical->valen,
                      cs_glob_lagr_physico_chemical->phi_p,
                      cs_glob_lagr_physico_chemical->phi_s,
                      cs_glob_lagr_physico_chemical->cstham,
                      cs_glob_lagr_clogging_model->csthpp,
                      cs_glob_lagr_physico_chemical->lambda_vdw);

  /* Initialization for the roughness surface model
     ---------------------------------------------- */

  if (lagr_model->roughness == 1)
    roughness_init(&cs_glob_lagr_physico_chemical->epseau,
                   &cs_glob_lagr_physico_chemical->fion,
                   tempp,
                   &cs_glob_lagr_physico_chemical->valen,
                   &cs_glob_lagr_physico_chemical->phi_p,
                   &cs_glob_lagr_physico_chemical->phi_s,
                   &cs_glob_lagr_physico_chemical->cstham,
                   &cs_glob_lagr_physico_chemical->lambda_vdw,
                   &cs_glob_lagr_reentrained_model->espasg,
                   &cs_glob_lagr_reentrained_model->denasp,
                   &cs_glob_lagr_reentrained_model->rayasp,
                   &cs_glob_lagr_reentrained_model->rayasg);

  /* Initialization for the clogging model
     ------------------------------------- */

  if (lagr_model->clogging == 1)
    cloginit(&cs_glob_lagr_physico_chemical->epseau,
             &cs_glob_lagr_physico_chemical->fion,
             &cs_glob_lagr_clogging_model->jamlim,
             &cs_glob_lagr_clogging_model->mporos,
             &cs_glob_lagr_clogging_model->diam_mean,
             tempp,
             &cs_glob_lagr_physico_chemical->valen,
             &cs_glob_lagr_physico_chemical->phi_p,
             &cs_glob_lagr_physico_chemical->phi_s,
             &cs_glob_lagr_physico_chemical->cstham,
             &cs_glob_lagr_clogging_model->csthpp,
             &cs_glob_lagr_physico_chemical->lambda_vdw);

  /* Update for new particles which entered the domain
     ------------------------------------------------- */

  /* At the first time step we initialize particles to
     values at current time step and not at previous time step, because
     values at previous time step = initialization */

  int  iprev;
  if (ts->nt_cur == 1) /* Use fields at current time step */
    iprev = 0;

  else                 /* Use fields at previous time step */
    iprev = 1;

  cs_lagr_injection(iprev, itypfb, vislen);

  /* Management of advancing time
     ---------------------------- */

  cs_glob_lagr_time_step->dtp = dt[0];
  cs_glob_lagr_time_step->ttclag += cs_glob_lagr_time_step->dtp;

  part_c = cs_lagr_update_particle_counter();

  if (part_c->n_g_total > 0) {

  /* Record particle's starting cell and rank, and update rebound time id */

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      cs_lagr_particles_set_lnum_n
        (p_set, ip, 1, CS_LAGR_CELL_NUM,
         cs_lagr_particles_get_lnum(p_set, ip, CS_LAGR_CELL_NUM));

      cs_lagr_particles_set_lnum_n(p_set, ip, 1, CS_LAGR_RANK_ID,
                                   cs_glob_rank_id);

      cs_lnum_t rebound_id
        = cs_lagr_particles_get_lnum(p_set, ip, CS_LAGR_REBOUND_ID);
      if (rebound_id >= 0)
        cs_lagr_particles_set_lnum(p_set, ip, CS_LAGR_REBOUND_ID,
                                   rebound_id + 1);

    }

    /* Pressure and fluid velocity gradients
       ------------------------------------- */

    /* At the first time step we initialize particles to
       values at current time step and not at previous time step, because
       values at previous time step = initialization (null gradients) */

    mode = 1;
    if (ts->nt_cur == 1)
      mode = 0;

    cs_lagr_gradients(mode, gradpr, gradvf);

    /* Particles progression
       --------------------- */

    bool go_on = true;
    cs_lnum_t n_particles_prev = p_set->n_particles - p_set->n_part_new;
    while (go_on) {

      cs_glob_lagr_time_step->nor
        = cs_glob_lagr_time_step->nor % cs_glob_lagr_time_scheme->t_order;
      cs_glob_lagr_time_step->nor++;

      /* Allocations     */

      cs_lnum_t nresnew = 0;

      cs_real_t   *taup;
      cs_real_33_t *bx;
      cs_real_3_t *tlag, *piil;
      BFT_MALLOC(taup, p_set->n_particles, cs_real_t);
      BFT_MALLOC(tlag, p_set->n_particles, cs_real_3_t);
      BFT_MALLOC(piil, p_set->n_particles, cs_real_3_t);
      BFT_MALLOC(bx, p_set->n_particles, cs_real_33_t);

      cs_real_t *tsfext = NULL;
      if (cs_glob_lagr_time_scheme->iilagr == 2)
        BFT_MALLOC(tsfext, p_set->n_particles, cs_real_t);

      cs_real_t *cpgd1 = NULL, *cpgd2 = NULL, *cpght = NULL;
      if (   cs_glob_lagr_time_scheme->iilagr == 2
          && lagr_model->physical_model == 2
          && cs_glob_lagr_source_terms->ltsthe == 1) {

        BFT_MALLOC(cpgd1, p_set->n_particles, cs_real_t);
        BFT_MALLOC(cpgd2, p_set->n_particles, cs_real_t);
        BFT_MALLOC(cpght, p_set->n_particles, cs_real_t);

      }

      cs_real_t *tempct = NULL;
      if (   (   lagr_model->physical_model == 1
              && cs_glob_lagr_specific_physics->itpvar == 1)
          || lagr_model->physical_model == 2)
        BFT_MALLOC(tempct, p_set->n_particles * 2, cs_real_t);

      cs_real_t *terbru = NULL;
      if (cs_glob_lagr_brownian->lamvbr == 1)
        BFT_MALLOC(terbru, p_set->n_particles, cs_real_t);

      /* Copy results from previous step */

      if (cs_glob_lagr_time_step->nor == 1) {
        /* Current to previous but not on new particles at the first time
         * because the user may have changed their position */
        for (cs_lnum_t ip = 0; ip < n_particles_prev; ip++)
          cs_lagr_particles_current_to_previous(p_set, ip);

        n_particles_prev = p_set->n_particles;
      }

      /* Computation of the fluid's pressure and velocity gradient
         at n+1 (with values at current time step) */
      if (   cs_glob_lagr_time_step->nor == 2
          && cs_glob_lagr_time_scheme->iilagr != 3)
        cs_lagr_gradients(0, gradpr, gradvf);

      /* use fields at previous or current time step */
      if (cs_glob_lagr_time_step->nor == 1)
        /* Use fields at previous time step    */
        iprev = 1;

      else
        iprev = 0;

      /* Retrieve bx values associated with particles from previous pass */

      if (   cs_glob_lagr_time_scheme->t_order == 2
          && cs_glob_lagr_time_step->nor == 2) {

        for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

          cs_real_t *jbx1 = cs_lagr_particles_attr(p_set, ip, CS_LAGR_TURB_STATE_1);

          for (cs_lnum_t ii = 0; ii < 3; ii++) {

            bx[ip][ii][0] = jbx1[ii];

          }

        }

      }

      cs_lagr_car(iprev,
                  dt,
                  taup,
                  tlag,
                  piil,
                  bx,
                  tempct,
                  gradpr,
                  gradvf,
                  w1,
                  w2);

      /* Integration of SDEs: position, fluid and particle velocity */

      cs_lagr_sde(cs_glob_lagr_time_step->dtp,
                  (const cs_real_t *)taup,
                  (const cs_real_3_t *)tlag,
                  (const cs_real_3_t *)piil,
                  (const cs_real_33_t *)bx,
                  tsfext,
                  (const cs_real_3_t *)gradpr,
                  (const cs_real_33_t *)gradvf,
                  terbru,
                  vislen,
                  &nresnew);

      /* Save bx values associated with particles for next pass */

      if (   cs_glob_lagr_time_scheme->t_order == 2
          && cs_glob_lagr_time_step->nor == 1) {

        for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

          cs_real_t *jbx1 = cs_lagr_particles_attr(p_set, ip, CS_LAGR_TURB_STATE_1);

          for (int  ii = 0; ii < 3; ii++)
            jbx1[ii] = bx[ip][ii][0];

        }

      }

      /* Integration of SDE related to physical models */

      if (lagr_model->physical_model == 1 || lagr_model->physical_model == 2) {

        if (cs_glob_lagr_time_step->nor == 1)
          /* Use fields at previous time step    */
          iprev   = 1;

        else
          /* Use fields at current time step     */
          iprev   = 0;

        cs_lagr_sde_model(tempct, cpgd1, cpgd2, cpght);

      }

      /* Integration of additional user variables
         -----------------------------------------*/

      if (cs_glob_lagr_model->n_user_variables > 0)
        cs_user_lagr_sde(dt, taup, tlag, tempct);

      /* Reverse coupling: compute source terms
         -------------------------------------- */

      if (   cs_glob_lagr_time_scheme->iilagr == 2
          && cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order)
        cs_lagr_coupling(taup, tempct, tsfext, cpgd1, cpgd2, cpght, w1, w2);

      /* Deallocate arrays whose size is based on p_set->n_particles
         (which may change next) */

      BFT_FREE(tlag);
      BFT_FREE(taup);
      BFT_FREE(piil);
      BFT_FREE(bx);

      if (cs_glob_lagr_time_scheme->iilagr == 2)
        BFT_FREE(tsfext);

      if (   cs_glob_lagr_time_scheme->iilagr == 2
          && lagr_model->physical_model == 2
          && cs_glob_lagr_source_terms->ltsthe == 1) {
        BFT_FREE(cpgd1);
        BFT_FREE(cpgd2);
        BFT_FREE(cpght);
      }

      if (   (   lagr_model->physical_model == 1
              && cs_glob_lagr_specific_physics->itpvar == 1)
          || lagr_model->physical_model == 2)
        BFT_FREE(tempct);

      if (cs_glob_lagr_brownian->lamvbr == 1)
        BFT_FREE(terbru);

      p_set->n_particles += nresnew;

      /* Location of particles - boundary conditions for particle positions
         ------------------------------------------------------------------ */

      if (cs_glob_lagr_time_step->nor == 1) {

        /* In unsteady case, reset boundary statistics */

        if (   cs_glob_lagr_time_scheme->isttio == 0
            || (   cs_glob_lagr_time_scheme->isttio == 1
                &&    cs_glob_time_step->nt_cur
                   <= cs_glob_lagr_stat_options->nstist)) {

          lag_bdi->tstatp = 0.0;
          lag_bdi->npstf  = 0;

          for (int  ii = 0; ii < cs_glob_lagr_dim->n_boundary_stats; ii++) {

            for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
              bound_stat[ii * n_b_faces + ifac] = 0.0;

          }

        }

        lag_bdi->tstatp += cs_glob_lagr_time_step->dtp;
        lag_bdi->npstf++;
        lag_bdi->npstft++;

        cs_lagr_tracking_particle_movement(vislen);

      }

      /* Update residence time */

      if (cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order) {

        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

          if (cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_CELL_NUM) - 1 >= 0) {
            cs_real_t res_time
              = cs_lagr_particles_get_real(p_set, npt,
                                           CS_LAGR_RESIDENCE_TIME)
                + cs_glob_lagr_time_step->dtp;
            cs_lagr_particles_set_real
              (p_set, npt, CS_LAGR_RESIDENCE_TIME, res_time);
          }

        }

      }

      /* Compute adhesion for reentrainement model
         ----------------------------------------- */

      if (lagr_model->resuspension > 0)
        cs_lagr_resuspension();

      /* Compute statistics
         ------------------ */

      /* Calculation of consolidation:
       * linear increase of the consolidation height with deposit time */

      if (cs_glob_lagr_consolidation_model->iconsol > 0) {

        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

          if (cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_DEPOSITION_FLAG) > 0) {

            cs_real_t p_depo_time
              = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_DEPO_TIME);
            cs_real_t p_consol_height
              = CS_MIN(cs_lagr_particles_get_real(p_set, npt, CS_LAGR_HEIGHT),
                       cs_glob_lagr_consolidation_model->rate_consol * p_depo_time );
            cs_lagr_particles_set_real(p_set, npt,
                                       CS_LAGR_CONSOL_HEIGHT, p_consol_height);
            cs_lagr_particles_set_real(p_set, npt, CS_LAGR_DEPO_TIME,
                                       p_depo_time + cs_glob_lagr_time_step->dtp);
          }
          else {
            cs_lagr_particles_set_real(p_set, npt, CS_LAGR_CONSOL_HEIGHT, 0.0);
          }

        }

      }

      if (   cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order
          && cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt)
        cs_lagr_stat_update();

      /* Statistics for clogging */

      if (   cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order
          && lagr_model->clogging == 1
          && cs_glob_lagr_consolidation_model->iconsol == 1) {

        /* Height and time of deposit     */

        for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

          bound_stat[lag_bdi->iclogt * n_b_faces + ifac] = 0.0;
          bound_stat[lag_bdi->iclogh * n_b_faces + ifac] = 0.0;
          bound_stat[lag_bdi->ihdiam * n_b_faces + ifac] = 0.0;

        }

        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

          if (cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_DEPOSITION_FLAG) == 1) {

            cs_lnum_t face_id = cs_lagr_particles_get_lnum
                                  (p_set, npt, CS_LAGR_NEIGHBOR_FACE_ID);

            cs_real_t p_diam = cs_lagr_particles_get_real
                                 (p_set, npt, CS_LAGR_DIAMETER);

            bound_stat[lag_bdi->iclogt * n_b_faces + face_id]
              += cs_lagr_particles_get_real(p_set, npt, CS_LAGR_DEPO_TIME);
            bound_stat[lag_bdi->iclogh * n_b_faces + face_id]
              +=  cs_lagr_particles_get_real(p_set, npt, CS_LAGR_CONSOL_HEIGHT)
                  * cs_math_pi * cs_math_sq(p_diam) * 0.25 / surfbn[face_id];

          }

        }

        for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

          if (bound_stat[lag_bdi->inclg * n_b_faces + ifac] > 0) {

            bound_stat[lag_bdi->iclogt * n_b_faces + ifac]
              /=  bound_stat[lag_bdi->inclg * n_b_faces + ifac];
            bound_stat[lag_bdi->ihdiam * n_b_faces + ifac]
              =   bound_stat[lag_bdi->ihsum  * n_b_faces + ifac]
                / bound_stat[lag_bdi->inclgt * n_b_faces + ifac];

          }
          else if (bound_stat[lag_bdi->inclg * n_b_faces + ifac] <= 0) {

            bound_stat[lag_bdi->iclogt * n_b_faces + ifac] = 0.0;
            bound_stat[lag_bdi->ihdiam * n_b_faces + ifac] = 0.0;

          }
          else {

            /* FIXME */

            bft_printf("   ** LAGRANGIAN MODULE:\n"
                       "   ** Error in cs_lagr.c: inclg < 0 ! \n"
                       "---------------------------------------------\n\n\n"
                       "** Ifac = %d  and inclg = %g\n"
                       "-------------------------------------------------\n",
                       ifac, bound_stat[lag_bdi->inclg * n_b_faces + ifac]);
          }

        }

      }

      /* Poisson equation
         ---------------- */

      if (   cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order
          && cs_glob_lagr_time_scheme->ilapoi == 1)
        cs_lagr_poisson(itypfb);

      /* Loop again ?
         ---------- */

      if (   cs_glob_lagr_time_scheme->t_order != 2
          || cs_glob_lagr_time_step->nor != 1)
        go_on = false; //exit the while loop

    }

  }  /* end if number of particles > 0 */

  /* Optional user modification of Lagrangian variables at end of iteration */

  cs_user_lagr_extra_operations(dt);

  /* Update particle counter */
  /*-------------------------*/

  part_c = cs_lagr_update_particle_counter();
  part_c->n_g_cumulative_total += part_c->n_g_new;
  part_c->n_g_cumulative_failed += part_c->n_g_failed;

  /* Logging
     ------- */

  int  modntl;
  if (ipass == 1)
    modntl = 0;

  else if (cs_glob_lagr_log_frequency_n > 0)
    modntl = ts->nt_cur % cs_glob_lagr_log_frequency_n;

  else if (cs_glob_lagr_log_frequency_n == -1 && ts->nt_cur == ts->nt_max)
    modntl = 0;

  else
    modntl = 1;

  if (modntl == 0)
    cs_lagr_print(ts->t_cur);

  /* Free memory */

  BFT_FREE(gradpr);
  BFT_FREE(w1);
  BFT_FREE(w2);
  if (cs_glob_lagr_time_scheme->modcpl > 0)
    BFT_FREE(gradvf);

  if (   lagr_model->deposition == 1
      && ts->nt_cur == ts->nt_max)
    BFT_FREE(vislen);

  if (   lagr_model->clogging == 1
      || lagr_model->roughness == 1
      || lagr_model->dlvo == 1)
    BFT_FREE(tempp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

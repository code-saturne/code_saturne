/*============================================================================
 * Methods for particle parameters
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

#include "cs_array.h"
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
#include "cs_lagr_aux_mean_fluid_quantities.h"
#include "cs_lagr_car.h"
#include "cs_lagr_coupling.h"
#include "cs_lagr_new.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_resuspension.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_print.h"
#include "cs_lagr_poisson.h"
#include "cs_lagr_post.h"
#include "cs_lagr_sde.h"
#include "cs_lagr_sde_model.h"
#include "cs_lagr_orientation.h"
#include "cs_lagr_prototypes.h"
#include "cs_lagr_agglo.h"
#include "cs_lagr_fragmentation.h"

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
  = {.iilagr = CS_LAGR_OFF,
     .isttio = 0,
     .isuila = 1,
     .t_order = 0,
     .extended_t_scheme = 0,
     .interpol_field = 1,
     .ilapoi = 0,
     .iadded_mass = 0,
     .added_mass_const = 0};

/* Main Lagrangian physical model parameters */

static cs_lagr_model_t  _lagr_model
  = {.physical_model = CS_LAGR_PHYS_OFF,
     .n_temperature_layers = 1,
     .modcpl = 1,
     .idistu = -1,
     .idiffl = -1,
     .deposition = 0,
     .dlvo = 0,
     .roughness = 0,
     .resuspension = 0,
     .clogging = 0,
     .shape = 0,
     .consolidation = 0,
     .precipitation = 0,
     .fouling = 0,
     .n_stat_classes = 0,
     .agglomeration = 0,
     .fragmentation = 0,
     .n_user_variables = 0,
     .viscous_terms = false};

/* particle counter structure and associated pointer */

static cs_lagr_particle_counter_t _lagr_particle_counter
  = {.n_g_cumulative_total = 0,
     .n_g_cumulative_failed = 0,
     .n_g_total = 0,
     .n_g_new = 0,
     .n_g_exit = 0,
     .n_g_merged = 0,
     .n_g_deposited = 0,
     .n_g_fouling = 0,
     .n_g_resuspended = 0,
     .n_g_failed = 0,
     .w_total = 0.,
     .w_new = 0.,
     .w_exit = 0.,
     .w_merged = 0.,
     .w_deposited = 0.,
     .w_fouling = 0.,
     .w_resuspended = 0.};

/* lagr specific physics structure and associated pointer */
static cs_lagr_specific_physics_t _cs_glob_lagr_specific_physics
  = {0, 0, 0, 0, 0};
cs_lagr_specific_physics_t *cs_glob_lagr_specific_physics
  = &_cs_glob_lagr_specific_physics;

/* lagr reentrained model structure and associated pointer */
static cs_lagr_reentrained_model_t _lagr_reentrained_model
  = {.ireent = 0,
     .iflow = 0,
     .espasg = 0.,
     .denasp = 0.,
     .modyeq = 0.,
     .rayasp = 0.,
     .rayasg = 0.};
cs_lagr_reentrained_model_t *cs_glob_lagr_reentrained_model
  = &_lagr_reentrained_model;

/* lagr precipitation model structure and associated pointer */
static cs_lagr_precipitation_model_t _cs_glob_lagr_precipitation_model
  = {0, 0, 0, NULL, NULL, NULL};
cs_lagr_precipitation_model_t *cs_glob_lagr_precipitation_model
  = &_cs_glob_lagr_precipitation_model;

/* lagr clogging model structure and associated pointer */
static cs_lagr_clogging_model_t _cs_glob_lagr_clogging_model = {0, 0, 0, 0};
cs_lagr_clogging_model_t *cs_glob_lagr_clogging_model
  = &_cs_glob_lagr_clogging_model;

/* lagr non-spherical model structure and associated pointer */
static cs_lagr_shape_model_t _cs_glob_lagr_shape_model = {0};
cs_lagr_shape_model_t *cs_glob_lagr_shape_model
  = &_cs_glob_lagr_shape_model;

/* lagr agglomeration model structure and associated pointer */
static cs_lagr_agglomeration_model_t _cs_glob_lagr_agglomeration_model
  = {
    .n_max_classes = 10000,
    .min_stat_weight = 0.,
    .max_stat_weight = 0.,
    .scalar_kernel = 0.,
    .base_diameter = 0.};

cs_lagr_agglomeration_model_t *cs_glob_lagr_agglomeration_model
  = &_cs_glob_lagr_agglomeration_model;

/* lagr fragmentation model structure and associated pointer */
static cs_lagr_fragmentation_model_t _cs_glob_lagr_fragmentation_model
  = {
    .scalar_kernel = 0.,
    .base_diameter = 0.,
    NULL};

cs_lagr_fragmentation_model_t *cs_glob_lagr_fragmentation_model
  = &_cs_glob_lagr_fragmentation_model;

/* lagr consolidation model structure and associated pointer */
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
  = {.npstf = 0,
     .npstft = 0,
     .has_part_impact_nbr = 0,
     .iclgst = 0,
     .inbr = -1,
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
     .cmu = 0,
     .visls0 = 0,
     .ustar = NULL,
     .cromf = NULL,
     .pressure = NULL,
     .scal_t = NULL,
     .temperature = NULL,
     .vel = NULL,
     .viscl = NULL,
     .cpro_viscls = NULL,
     .cpro_cp = NULL,
     .rad_energy = NULL,
     .x_oxyd = NULL,
     .x_eau = NULL,
     .x_m = NULL,
     .cvar_k = NULL,
     .cvar_ep = NULL,
     .cvar_omg = NULL,
     .cvar_rij = NULL,
     .grad_pr = NULL,
     .grad_vel = NULL,
     .lagr_time = NULL,
     .grad_lagr_time = NULL};

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

/*! Projection matrices for global to local coordinates on boundary faces */

cs_real_33_t  *cs_glob_lagr_b_face_proj = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_lagr_params_pointers(int  **p_iilagr,
                          int  **p_idepst,
                          int  **p_iflow,
                          int  **p_ipreci);

void
cs_f_lagr_dim_pointers(int  **p_ntersl);

void
cs_f_lagr_clogging_model_pointers(cs_real_t  **jamlim,
                                  cs_real_t  **mporos,
                                  cs_real_t  **csthpp);

void
cs_f_lagr_shape_model_pointers(cs_real_t **param_chmb);

void
cs_f_lagr_agglomeration_model_pointers( cs_lnum_t  **n_max_classes,
                                        cs_real_t  **min_stat_weight,
                                        cs_real_t  **max_stat_weight,
                                        cs_real_t  **scalar_kernel,
                                        cs_real_t  **base_diameter );

void
cs_f_lagr_consolidation_model_pointers(cs_lnum_t  **iconsol,
                                       cs_real_t  **rate_consol,
                                       cs_real_t  **slope_consol,
                                       cs_real_t  **force_consol);

void
cs_f_lagr_source_terms_pointers(int  **p_ltsdyn,
                                int  **p_ltsmas,
                                int  **p_ltsthe,
                                int  **p_itsli,
                                int  **p_itske,
                                int  **p_itste,
                                int  **p_itsti,
                                int  **p_itsmas);

void
cs_f_lagr_specific_physics(int        *iirayo,
                           int        *ncharb,
                           int        *ncharm);

void
cs_f_lagr_coal_comb(int        *ih2o,
                    int        *io2,
                    int        *ico,
                    int        *iatc,
                    cs_real_t  *prefth,
                    cs_real_t  *trefth,
                    int        *natom,
                    cs_real_t  *wmolat,
                    int        *ngazem,
                    cs_real_t  *wmole,
                    int        *iym1,
                    int        *ncharm,
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
 *   p_iflow  -->
 *   p_ipreci --> precipitation option flag
 *----------------------------------------------------------------------------*/

void
cs_f_lagr_params_pointers(int  **p_iilagr,
                          int  **p_idepst,
                          int  **p_iflow,
                          int  **p_ipreci)
{
  *p_iilagr = &_lagr_time_scheme.iilagr;
  *p_idepst = &_lagr_model.deposition;
  *p_iflow= &_lagr_reentrained_model.iflow;
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
cs_f_lagr_dim_pointers(int  **p_ntersl)
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
cs_f_lagr_shape_model_pointers(cs_real_t **param_chmb)
{
  *param_chmb = &cs_glob_lagr_shape_model->param_chmb;
}

void
cs_f_lagr_agglomeration_model_pointers(cs_lnum_t **n_max_classes,
                                       cs_real_t **min_stat_weight,
                                       cs_real_t **max_stat_weight,
                                       cs_real_t **scalar_kernel,
                                       cs_real_t **base_diameter )
{
  *n_max_classes    = &cs_glob_lagr_agglomeration_model->n_max_classes;
  *min_stat_weight = &cs_glob_lagr_agglomeration_model->min_stat_weight;
  *max_stat_weight = &cs_glob_lagr_agglomeration_model->max_stat_weight;
  *scalar_kernel   = &cs_glob_lagr_agglomeration_model->scalar_kernel;
  *base_diameter   = &cs_glob_lagr_agglomeration_model->base_diameter;
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
cs_f_lagr_source_terms_pointers(int  **p_ltsdyn,
                                int  **p_ltsmas,
                                int  **p_ltsthe,
                                int  **p_itsli,
                                int  **p_itske,
                                int  **p_itste,
                                int  **p_itsti,
                                int  **p_itsmas)
{
  *p_ltsdyn = &cs_glob_lagr_source_terms->ltsdyn;
  *p_ltsmas = &cs_glob_lagr_source_terms->ltsmas;
  *p_ltsthe = &cs_glob_lagr_source_terms->ltsthe;
  *p_itsli  = &cs_glob_lagr_source_terms->itsli;
  *p_itske  = &cs_glob_lagr_source_terms->itske;
  *p_itste  = &cs_glob_lagr_source_terms->itste;
  *p_itsti  = &cs_glob_lagr_source_terms->itsti;
  *p_itsmas = &cs_glob_lagr_source_terms->itsmas;
}

void
cs_f_lagr_specific_physics(int        *iirayo,
                           int        *ncharb,
                           int        *ncharm)
{
  cs_turb_model_t  *turb_model = cs_get_glob_turb_model();

  if (turb_model == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Turbulence modelling is not set.", __func__);

  _lagr_extra_module.iturb  = turb_model->iturb;
  _lagr_extra_module.itytur = turb_model->itytur;
  if (ncharb != NULL)
    _lagr_extra_module.ncharb = *ncharb;
  if (ncharm != NULL)
    _lagr_extra_module.ncharm = *ncharm;
  _lagr_extra_module.icp    = cs_glob_fluid_properties->icp;

  _lagr_extra_module.radiative_model = *iirayo;
  _lagr_extra_module.cmu    = cs_turb_cmu;
}

void
cs_f_lagr_coal_comb(int        *ih2o,
                    int        *io2,
                    int        *ico,
                    int        *iatc,
                    cs_real_t  *prefth,
                    cs_real_t  *trefth,
                    int        *natom,
                    cs_real_t  *wmolat,
                    int        *ngazem,
                    cs_real_t  *wmole,
                    int        *iym1,
                    int        *ncharm,
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
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0) {
    _lagr_extra_module.cromf   = cs_field_by_name_try("rho_gas");
  }
  else {
    _lagr_extra_module.cromf   = cs_field_by_name_try("density");
  }

  _lagr_extra_module.pressure  = cs_field_by_name_try("pressure");

  _lagr_extra_module.rad_energy = cs_field_by_name_try("rad_energy");

  _lagr_extra_module.lagr_time = cs_field_by_name_try("lagr_time");

  _lagr_extra_module.cvar_k = cs_field_by_name_try("k");
  /* using LES */
  if (_lagr_extra_module.cvar_k == NULL)
    _lagr_extra_module.cvar_k = cs_field_by_name_try("k_sgs");
  if (_lagr_extra_module.cvar_k == NULL)
    _lagr_extra_module.cvar_k = cs_field_by_name_try("lagr_k");

  _lagr_extra_module.cvar_ep = cs_field_by_name_try("epsilon");
  /* using LES */
  if (_lagr_extra_module.cvar_ep == NULL)
    _lagr_extra_module.cvar_ep = cs_field_by_name_try("epsilon_sgs");
  if (_lagr_extra_module.cvar_ep == NULL)
    _lagr_extra_module.cvar_ep = cs_field_by_name_try("lagr_epsilon");

  if (cs_field_by_name_try("velocity_1") != NULL) {
    /* we are probably using NEPTUNE_CFD */
    _lagr_extra_module.vel         = cs_field_by_name_try("lagr_velocity");

    _lagr_extra_module.cvar_omg    = NULL;
    _lagr_extra_module.cvar_rij    = cs_field_by_name_try("lagr_rij");
    _lagr_extra_module.viscl       = cs_field_by_name_try
                                       ("lagr_molecular_viscosity");
    _lagr_extra_module.scal_t      = cs_field_by_name_try("lagr_enthalpy");
    _lagr_extra_module.cpro_viscls = cs_field_by_name_try
                                       ("lagr_thermal_conductivity");
    _lagr_extra_module.cpro_cp     = cs_field_by_name_try("lagr_specific_heat");
    _lagr_extra_module.temperature = cs_field_by_name_try("lagr_temperature");
    _lagr_extra_module.x_oxyd      = NULL;
    _lagr_extra_module.x_eau       = NULL;
    _lagr_extra_module.x_m         = NULL;
    _lagr_extra_module.cromf       = cs_field_by_name_try("lagr_density");
    /* TODO FIXME */
    _lagr_extra_module.visls0      = 0.;

    _lagr_extra_module.ustar
      = cs_field_by_name_try("lagr_wall_friction_velocity");
  }
  else {
    /* we use code_saturne */
    _lagr_extra_module.vel         = cs_field_by_name_try("velocity");
    _lagr_extra_module.cvar_omg    = cs_field_by_name_try("omega");
    _lagr_extra_module.cvar_rij    = cs_field_by_name_try("rij");
    _lagr_extra_module.viscl       = cs_field_by_name_try("molecular_viscosity");
    _lagr_extra_module.cpro_viscls = NULL;

    _lagr_extra_module.scal_t = cs_thermal_model_field();

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

    _lagr_extra_module.x_oxyd      = cs_field_by_name_try("ym_o2");
    _lagr_extra_module.x_eau       = cs_field_by_name_try("ym_h2o");
    _lagr_extra_module.x_m         = cs_field_by_name_try("xm");

    _lagr_extra_module.ustar  = cs_field_by_name_try("boundary_ustar");
    if (_lagr_extra_module.ustar == NULL)
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Obtain the number of mesh cells occupied by at least one particle.
 *
 * \param[in]   p_set   pointer to particle data structure
 * \param[in]   start   start position in p_set
 * \param[in]   end     end position in p_set
 *
 * \return
 *   integer that gives the number of cells occupied by at least one particle
 *   in a particle sub-set
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_get_n_occupied_cells(const cs_lagr_particle_set_t  *p_set,
                      cs_lnum_t                      start,
                      cs_lnum_t                      end)
{
  /*Initialization:
    counter for number of cells that contain particles */
  cs_lnum_t counter_particle_cells = 0;

  /* Main loop to update the counter */
  if (end - start > 1) {
    counter_particle_cells = 1;
    cs_lnum_t prev_cell_id
      = cs_lagr_particles_get_lnum(p_set, start, CS_LAGR_CELL_ID);
    cs_lnum_t curr_cell_id = prev_cell_id;
    for (cs_lnum_t p = start + 1; p < end; ++p) {
      curr_cell_id = cs_lagr_particles_get_lnum(p_set, p, CS_LAGR_CELL_ID);
      if (prev_cell_id != curr_cell_id)
        counter_particle_cells++;
      prev_cell_id = curr_cell_id;
    }
  }
  else if (p_set->n_particles == 1) {
    counter_particle_cells = 1;
  }

  return counter_particle_cells;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Obtain the index of cells that are occupied by at least one
 *         particle and the list of indices of the particle set that contain
 *         the first element in a cell
 *
 * The two arrays that contain this information need to be pre-allocated to
 * size n_occupied_cells and n_occupied_cells+1 respectively.
 *
 * The value of particle_gaps at index i contains the first particle in
 * in cell[i+1]. The last element of particle_gaps contains the size of p_set.
 *
 * \param[in]   p_set              id of associated location
 * \param[in]   start              id of associated zone
 * \param[in]   end                id of requested class
 * \param[in]   n_occupied_cells   number of cells that are occupied by particles
 * \param[out]  occupied_cell_ids  occupied_cell ids
 * \param[out]  particle_gaps      preallocated list of size n_occupied_cells+1
 *                                 that will contain the starting indices of
 *                                 particles that are in the current cell
 */
/*----------------------------------------------------------------------------*/

static void
_occupied_cells(cs_lagr_particle_set_t  *p_set,
                cs_lnum_t                start,
                cs_lnum_t                end,
                cs_lnum_t                n_occupied_cells,
                cs_lnum_t                occupied_cell_ids[],
                cs_lnum_t                particle_gaps[])
{
  if (end - start >= 1) {
    /* Initialization */
    cs_lnum_t prev_cell_id
      = cs_lagr_particles_get_lnum(p_set, start, CS_LAGR_CELL_ID);
    cs_lnum_t curr_cell_id = prev_cell_id;
    cs_lnum_t counter = 0;
    occupied_cell_ids[0] = curr_cell_id;
    particle_gaps[0] = 0;

    counter = 1;

    /* Update lists */
    for (cs_lnum_t part = start + 1; part < end; ++part) {
      curr_cell_id
        = cs_lagr_particles_get_lnum(p_set, part, CS_LAGR_CELL_ID);
      if (prev_cell_id != curr_cell_id) {
        occupied_cell_ids[counter] = curr_cell_id;
        particle_gaps[counter] = part;
        counter++;
      }
      prev_cell_id = curr_cell_id;
    }
    particle_gaps[n_occupied_cells] = p_set->n_particles;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Obtain the number of particles to be deleted
 *
 * \param[in]        p_set             pointer to particle data structure
 * \param[in]        start             start position in p_set
 * \param[in]        end               end position in p_se
 *
 * \return:
 *   integer that gives the number of particles to be deleted in a particle sub-set
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_get_n_deleted(cs_lagr_particle_set_t  *p_set,
               cs_lnum_t                start,
               cs_lnum_t                end)
{
  cs_lnum_t res = 0;

  for (cs_lnum_t idx = start; idx < end; ++idx) {
    if (cs_lagr_particles_get_flag(p_set, idx, CS_LAGR_PART_TO_DELETE))
      res++;
  }
  return res;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize lagrangian arrays
 *----------------------------------------------------------------------------*/

void
cs_lagr_init_arrays(void)
{
  cs_lnum_t  n_b_faces = cs_glob_mesh->n_b_faces;
  int   n_boundary_stats = cs_glob_lagr_dim->n_boundary_stats;

  assert(bound_stat == NULL);

  if (n_boundary_stats > 0)
    BFT_MALLOC(bound_stat, n_b_faces * n_boundary_stats, cs_real_t);

  BFT_MALLOC(cs_glob_lagr_source_terms->st_val,
             cs_glob_lagr_dim->ntersl * cs_glob_mesh->n_cells_with_ghosts,
             cs_real_t);
  for (cs_lnum_t i = 0; i < cs_glob_lagr_dim->ntersl; i++) {
    cs_real_t *st =   cs_glob_lagr_source_terms->st_val
                   + i*cs_glob_mesh->n_cells_with_ghosts;
    cs_array_real_fill_zero(cs_glob_mesh->n_cells_with_ghosts, st);
  }
}

/*----------------------------------------------------------------------------
 * Free lagrangian arrays
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

  BFT_FREE(cs_glob_lagr_b_face_proj);

  /* encrustation pointers */

  BFT_FREE(cs_glob_lagr_encrustation->enc1);
  BFT_FREE(cs_glob_lagr_encrustation->enc2);
  BFT_FREE(cs_glob_lagr_encrustation->tprenc);
  BFT_FREE(cs_glob_lagr_encrustation->visref);

  /* boundary interaction pointers */

  for (int i = 0; i < cs_glob_lagr_dim->n_boundary_stats; i++)
    BFT_FREE(cs_glob_lagr_boundary_interactions->nombrd[i]);
  BFT_FREE(cs_glob_lagr_boundary_interactions->nombrd);

  /* Statistics */

  cs_lagr_stat_finalize();

  /* Also close log file (TODO move this) */

  cs_lagr_print_finalize();

  /* Close tracking structures */

  cs_lagr_tracking_finalize();

  cs_lagr_finalize_zone_conditions();

  /* Fluid gradients */
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  BFT_FREE(extra->grad_pr);
  if (extra->grad_vel != NULL)
    BFT_FREE(extra->grad_vel);

  if (extra->grad_lagr_time != NULL)
    BFT_FREE(extra->grad_lagr_time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create additional fields needed by the Lagrangien model
 *
 * Most additional fields can be defined directly in
 * \ref cs_lagr_options_definition, but some fields may be mapped to
 * different fields based on the calling module (i.e. code_saturne or
 * neptune_cfd), and possibly defined after that call.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_add_fields(void)
{
  if (cs_glob_lagr_time_scheme->iilagr <= CS_LAGR_OFF)
    return;

  /* Add Lagrangian integral time for Lagrangian computation */

  const int k_vis = cs_field_key_id("post_vis");
  const int k_log = cs_field_key_id("log");

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_field_t *f = NULL;

  /* Add Lagrangian integral time */

  f = cs_field_create("lagr_time",
                      CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                      CS_MESH_LOCATION_CELLS,
                      1,
                      false);

  cs_field_set_key_int(f, k_log, 1);
  cs_field_set_key_int(f, k_vis, 1);

  /* Add TKE for DRSM model */

  if (extra->itytur == 3 || extra->itytur == 4) {
    f = cs_field_by_name_try("k");
    if (f == NULL)
      f = cs_field_by_name_try("k_sgs");
    if (f == NULL)
      f = cs_field_by_name_try("lagr_k");
    if (f == NULL) {
      f = cs_field_find_or_create("k_sgs",
                                  CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                  CS_MESH_LOCATION_CELLS,
                                  1,
                                  true);
      cs_field_set_key_int(f, k_log, 1);
    }
  }

  /* Add Dissipation for LES or k-omega */

  f = cs_field_by_name_try("epsilon");
  if (f == NULL)
    f = cs_field_by_name_try("epsilon_sgs");
  if (f == NULL)
    f = cs_field_by_name_try("lagr_epsilon");
  if (f == NULL) {
    f = cs_field_find_or_create("epsilon_sgs",
                                CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                CS_MESH_LOCATION_CELLS,
                                1,
                                true);
    cs_field_set_key_int(f, k_log, 1);
  }
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

  /* Fluid velocity by default */
  zis->velocity_profile     = -1;
  /* Fluid temperature by default  */
  zis->temperature_profile  = 0;

  if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL)
    zis->coal_number        = -2;

  zis->cluster              =  0;

  /* Agglomeration/fragmentation default*/
  zis->aggregat_class_id = 1;
  zis->aggregat_fractal_dim = 3;

  zis->velocity_magnitude   = - cs_math_big_r;

  for (int  i = 0; i < 3; i++)
    zis->velocity[i]        = - cs_math_big_r;

  /* For spheroids without inertia  */
  /* Default shape: sphere */
  zis->shape = CS_LAGR_SHAPE_SPHERE_MODEL;

  /* Angular velocity */
  for (int i = 0; i < 3; i++)
    zis->angular_vel[i] = 0.;

  /* Spheroids radii a b c */
  for (int i = 0; i < 3; i++)
    zis->radii[i] = 0.;

  /* Shape parameters */
  for (int i = 0; i < 4; i++)
    zis->shape_param[i] = 0.;

  /* First three Euler parameters */
  for (int i = 0; i < 3; i++)
    zis->euler[i] = 0.;

  zis->euler[3] = 1.;

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

  cs_gnum_t gcount[] = {p_set->n_particles,
                        p_set->n_part_new,
                        p_set->n_part_merged,
                        p_set->n_part_out,
                        p_set->n_part_dep,
                        p_set->n_part_fou,
                        p_set->n_part_resusp,
                        p_set->n_failed_part};

  cs_real_t wsum[] = {p_set->weight,
                      p_set->weight_new,
                      p_set->weight_merged,
                      p_set->weight_out,
                      p_set->weight_dep,
                      p_set->weight_fou,
                      p_set->weight_resusp};

  cs_lnum_t size_count = sizeof(gcount) / sizeof(gcount[0]);
  cs_lnum_t size_sum = sizeof(wsum) / sizeof(wsum[0]);

  cs_parall_counter(gcount, size_count);
  cs_parall_sum(size_sum, CS_REAL_TYPE, wsum);

  pc->n_g_total = gcount[0];
  pc->n_g_new = gcount[1];
  pc->n_g_merged = gcount[2];
  pc->n_g_exit = gcount[3];
  pc->n_g_deposited = gcount[4];
  pc->n_g_fouling = gcount[5];
  pc->n_g_resuspended = gcount[6];
  pc->n_g_failed = gcount[7];

  pc->w_total = wsum[0];
  pc->w_new = wsum[1];
  pc->w_merged = wsum[2];
  pc->w_exit = wsum[3];
  pc->w_deposited = wsum[4];
  pc->w_fouling = wsum[5];
  pc->w_resuspended = wsum[6];

  return pc;
}

/*----------------------------------------------------------------------------*/
/*! \brief Provide access to cs_lagr_specific_physics_t
 *
 * \return
 *   pointer to lagrangian specific physics options
 */
/*----------------------------------------------------------------------------*/

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
 */
/*----------------------------------------------------------------------------*/

cs_lagr_reentrained_model_t *
cs_get_lagr_reentrained_model(void)
{
  return &_lagr_reentrained_model;
}

/*----------------------------------------------------------------------------*/
/*! \brief Provide access to cs_lagr_precipitation_model_t
 *
 * \return  pointer to lagrangian precipitation model options
 */
/*----------------------------------------------------------------------------*/

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
 * Provide access to cs_lagr_shape_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_shape_model_t *
cs_get_lagr_shape_model(void)
{
  return &_cs_glob_lagr_shape_model;
}

/*----------------------------------------------------------------------------
 * Provide access to cs_lagr_agglomeration_model_t
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_lagr_agglomeration_model_t *
cs_get_lagr_agglomeration_model(void)
{
  return &_cs_glob_lagr_agglomeration_model;
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
 * The structure is allocated on demand, when this function is first called.
 *
 * \return pointer to current internal_conditions structure
 */
/*----------------------------------------------------------------------------*/

cs_lagr_internal_condition_t  *
cs_lagr_get_internal_conditions(void)
{
  if (cs_glob_mesh->time_dep >= CS_MESH_TRANSIENT_CONNECT) {
    bft_error(__FILE__, __LINE__, 0,
              "%s: Lagrangian internal conditions are not currenty\n"
              "compatible with a transient (time-varying) mesh.", __func__);
  }

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

  /* Allocate pressure and velocity gradients */
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lnum_t ncelet = cs_glob_mesh->n_cells_with_ghosts;

  if (   cs_glob_lagr_time_scheme->extended_t_scheme !=0
      && cs_glob_lagr_model->idistu == 1)
      BFT_MALLOC(extra->grad_lagr_time, ncelet, cs_real_3_t);

  if (   cs_glob_lagr_model->modcpl > 0
      || cs_glob_lagr_model->shape > 0
      || cs_glob_lagr_time_scheme->interpol_field != 0)
      BFT_MALLOC(extra->grad_vel, ncelet, cs_real_33_t);

  /* For frozen field:
     values at previous time step = values at current time step */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_FROZEN_CONTINUOUS_PHASE) {

    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {

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

  if (cs_glob_lagr_time_scheme->iilagr != CS_LAGR_OFF)
    cs_lagr_restart_read_p();

  /* Read statistics restart data */

  cs_lagr_stat_restart_read();
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
  const cs_time_step_t *ts = cs_glob_time_step;
  const cs_mesh_t *mesh = cs_glob_mesh;

  cs_lagr_boundary_interactions_t *lag_bdi = cs_glob_lagr_boundary_interactions;
  cs_lagr_model_t *lagr_model = cs_glob_lagr_model;
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_particle_counter_t *part_c = cs_lagr_get_particle_counter();

  cs_lnum_t n_b_faces = mesh->n_b_faces;

  cs_lnum_t *ifabor = mesh->b_face_cells;
  cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

  /* Allocate arrays depending on user options */

  cs_real_t *tempp = NULL;
  if (   lagr_model->clogging == 1
      || lagr_model->roughness == 1
      || lagr_model->dlvo == 1)
    BFT_MALLOC(tempp, mesh->n_cells, cs_real_t);

  ipass ++;

  static cs_real_t *vislen = NULL;
  if ((lagr_model->deposition == 1) && (ipass == 1)) {
    BFT_MALLOC(vislen, n_b_faces, cs_real_t);
    for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++)
      vislen[ifac] = -cs_math_big_r;
  }

  /* mask for for deposition-related states */

  int deposition_mask =   CS_LAGR_PART_DEPOSITED | CS_LAGR_PART_ROLLING
                        | CS_LAGR_PART_IMPOSED_MOTION;

  /* Initialization
     -------------- */

  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;

  /* First call (initializations)
     ---------------------------- */

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1) {

    /* If the deposition model is activated */

    if (lagr_model->deposition > 0) {

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

          cs_real_t surfb = b_face_surf[ifac];
          ustarmoy     = ustarmoy + surfb * _ustar;
          surftot      = surftot + surfb;
          vislen[ifac] = visccf / _ustar; // nu /u*
          //FIXME to be coherent with wall fn: y/y+
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

  /* Prepare statistics for this time step */

  cs_lagr_stat_prepare();

  /* Update boundary condition types;
     in most cases, this should be useful only at the first iteration,
     but we prefer to be safe in case of advanced uses with time-dependent
     zones or zone types */

  /* User initialization by set and boundary */

  if (ts->nt_cur == ts->nt_prev +1)
    cs_gui_particles_bcs();

  /* User initialization by set and zone */

  cs_user_lagr_boundary_conditions(itypfb);
  cs_user_lagr_volume_conditions();

  _update_boundary_face_type();

  /* Initialize counter */

  part_c->n_g_total = 0;
  part_c->n_g_new = 0;
  part_c->n_g_exit = 0;
  part_c->n_g_merged = 0;
  part_c->n_g_deposited = 0;
  part_c->n_g_fouling = 0;
  part_c->n_g_resuspended = 0;
  part_c->n_g_failed = 0;
  part_c->w_total = 0;
  part_c->w_new = 0;
  part_c->w_exit = 0;
  part_c->w_merged = 0;
  part_c->w_deposited = 0;
  part_c->w_fouling = 0;
  part_c->w_resuspended = 0;

  /* particles->n_part_new: handled in injection step */
  p_set->weight = 0.0;
  p_set->n_part_out = 0;
  p_set->n_part_dep = 0;
  p_set->n_part_fou = 0;
  p_set->weight_out = 0.0;
  p_set->weight_dep = 0.0;
  p_set->weight_fou = 0.0;
  p_set->n_failed_part = 0;
  p_set->weight_failed = 0.0;

  /* Initialization for the dlvo, roughness and clogging  model
     ---------------------------------------------------------- */

  if (   lagr_model->dlvo == 1
      || lagr_model->roughness == 1
      || lagr_model->clogging == 1) {

    if (extra->temperature != NULL) {

      if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS) {
        for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++) {
          tempp[iel] =    extra->temperature->val[iel]
                        + cs_physical_constants_celsius_to_kelvin;
        }
      }
      else {
        for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++)
          tempp[iel] =    extra->temperature->val[iel];
      }

    }

    else {
      if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS) {
        for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++) {
          tempp[iel] =    cs_glob_fluid_properties->t0
                        + cs_physical_constants_celsius_to_kelvin;
        }
      }
      else {
        for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++)
          tempp[iel] = cs_glob_fluid_properties->t0;
      }
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

  /* Initialization for the nonsphere model
     ------------------------------------- */

  if (lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_STOC_MODEL)
    cs_glob_lagr_shape_model->param_chmb = 1.0;

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

  /* Fields at current time step if using NEPTUNE_CFD */
  if (cs_field_by_name_try("velocity_1") != NULL)
    iprev = 0;

  cs_lagr_injection(iprev, itypfb, vislen);

  /* Initialization for the agglomeration/fragmentation models
     --------------------------------------------------------- */

  /* The algorithms are based on discrete values of particle radii
     (CS_LAGR_DIAMETER).
     It uses a minimum diameter, corresponding to monomers
     (unbreakable particles).
     Aggregates are formed of multiples of these monomers. */

  /* Evaluation of the minimum diameter */
  cs_real_t minimum_particle_diam = 0.;

  if (cs_glob_lagr_model->agglomeration == 1)
    minimum_particle_diam = cs_glob_lagr_agglomeration_model->base_diameter;

  if (cs_glob_lagr_model->fragmentation == 1)
    minimum_particle_diam = cs_glob_lagr_fragmentation_model->base_diameter;

  /* Management of advancing time
     ---------------------------- */

  cs_glob_lagr_time_step->dtp = dt[0];
  cs_glob_lagr_time_step->ttclag += cs_glob_lagr_time_step->dtp;

  part_c = cs_lagr_update_particle_counter();

  /* Global number of particles > 0 */
  if (part_c->n_g_total > 0) {

  /* Record particle's starting cell and rank, and update rebound time id */

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      cs_lagr_particles_set_lnum_n
        (p_set, ip, 1, CS_LAGR_CELL_ID,
         cs_lagr_particles_get_lnum(p_set, ip, CS_LAGR_CELL_ID));

      cs_lagr_particles_set_lnum_n(p_set, ip, 1, CS_LAGR_RANK_ID,
                                   cs_glob_rank_id);

      cs_lnum_t rebound_id
        = cs_lagr_particles_get_lnum(p_set, ip, CS_LAGR_REBOUND_ID);
      if (rebound_id >= 0)
        cs_lagr_particles_set_lnum(p_set, ip, CS_LAGR_REBOUND_ID,
                                   rebound_id + 1);

    }

    /* Compute the Lagrangian time */
    /* Pressure fluid velocity and Lagrangian time  gradients
       ----------------------------------------------------- */

    /* At the first time step we initialize particles to
       values at current time step and not at previous time step, because
       values at previous time step = initialization (zero gradients) */

    if (extra->itytur == 3 && extra->cvar_k->n_time_vals > 1) {
      /* save previous value dor the kinetic energy */
      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++)
        extra->cvar_k->vals[1][cell_id]
          = 0.5 * (  extra->cvar_rij->vals[1][6*cell_id]
                   + extra->cvar_rij->vals[1][6*cell_id + 1]
                   + extra->cvar_rij->vals[1][6*cell_id + 2]);
    }

    /* First pass allocate and compute it */
    if (extra->grad_pr == NULL) {
      BFT_MALLOC(extra->grad_pr, cs_glob_mesh->n_cells_with_ghosts, cs_real_3_t);

      cs_lagr_aux_mean_fluid_quantities(extra->lagr_time,
                                        extra->grad_pr,
                                        extra->grad_vel,
                                        extra->grad_lagr_time);
    }
    else if (   cs_glob_lagr_time_scheme->iilagr
             != CS_LAGR_FROZEN_CONTINUOUS_PHASE) {

      if (mesh->time_dep >= CS_MESH_TRANSIENT_CONNECT) {
        cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
        BFT_REALLOC(extra->grad_pr, n_cells_ext, cs_real_3_t);
        if (extra->grad_vel != NULL)
          BFT_REALLOC(extra->grad_vel, n_cells_ext, cs_real_33_t);
      }

      cs_lagr_aux_mean_fluid_quantities(extra->lagr_time,
                                        extra->grad_pr,
                                        extra->grad_vel,
                                        extra->grad_lagr_time);
    }

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

      cs_real_t *taup;
      cs_real_33_t *bx;
      cs_real_3_t *tlag, *piil;
      BFT_MALLOC(taup, p_set->n_particles, cs_real_t);
      BFT_MALLOC(tlag, p_set->n_particles, cs_real_3_t);
      BFT_MALLOC(piil, cs_glob_mesh->n_cells, cs_real_3_t);
      BFT_MALLOC(bx, p_set->n_particles, cs_real_33_t);

      cs_real_t *tsfext = NULL;
      if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING)
        BFT_MALLOC(tsfext, p_set->n_particles, cs_real_t);

      cs_real_3_t *beta = NULL;
      if (cs_glob_lagr_time_scheme->extended_t_scheme != 0)
        BFT_MALLOC(beta, p_set->n_particles, cs_real_3_t);

      cs_real_t *cpgd1 = NULL, *cpgd2 = NULL, *cpght = NULL;
      if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
          && lagr_model->physical_model == CS_LAGR_PHYS_COAL
          && cs_glob_lagr_source_terms->ltsthe == 1) {

        BFT_MALLOC(cpgd1, p_set->n_particles, cs_real_t);
        BFT_MALLOC(cpgd2, p_set->n_particles, cs_real_t);
        BFT_MALLOC(cpght, p_set->n_particles, cs_real_t);

      }

      cs_real_t *tempct = NULL;
      if (   (   lagr_model->physical_model == CS_LAGR_PHYS_HEAT
              && cs_glob_lagr_specific_physics->itpvar == 1)
          || lagr_model->physical_model == CS_LAGR_PHYS_COAL)
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
          && cs_glob_lagr_time_scheme->iilagr != CS_LAGR_FROZEN_CONTINUOUS_PHASE)
        cs_lagr_aux_mean_fluid_quantities(extra->lagr_time,
                                          extra->grad_pr,
                                          extra->grad_vel,
                                          extra->grad_lagr_time);

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

          cs_real_t *jbx1 = cs_lagr_particles_attr(p_set, ip,
                                                   CS_LAGR_TURB_STATE_1);

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
                  beta,
                  extra->grad_pr,
                  extra->grad_vel,
                  extra->grad_lagr_time);

      /* Integration of SDEs: position, fluid and particle velocity */

      cs_lagr_sde(cs_glob_lagr_time_step->dtp,
                  (const cs_real_t *)taup,
                  (const cs_real_3_t *)tlag,
                  (const cs_real_3_t *)piil,
                  (const cs_real_33_t *)bx,
                  tsfext,
                  (const cs_real_3_t *)extra->grad_pr,
                  (const cs_real_33_t *)extra->grad_vel,
                  terbru,
                  (const cs_real_t *)vislen,
                  beta,
                  &nresnew);

      /* Integration of SDEs for orientation of spheroids without inertia */
      if (lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_STOC_MODEL) {
        cs_lagr_orientation_dyn_spheroids(iprev,
                                          cs_glob_lagr_time_step->dtp,
                                          (const cs_real_33_t *)extra->grad_vel);
      }
      /* Integration of Jeffrey equations for ellispoids */
      else if (lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_JEFFERY_MODEL) {
        cs_lagr_orientation_dyn_jeffery(cs_glob_lagr_time_step->dtp,
                                        (const cs_real_33_t *)extra->grad_vel);
      }


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

      if (lagr_model->physical_model == CS_LAGR_PHYS_HEAT
          || lagr_model->physical_model == CS_LAGR_PHYS_COAL) {

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

      /* Integration of agglomeration and fragmentation
         ----------------------------------------------*/

      /* Agglomeration and fragmentation preparation */

      /* Preparation: find cells occupied by particles (number)
                      generate lists of these cells
                      generate list particles indexes (sublists within a cell) */

      cs_lnum_t n_occupied_cells;

      cs_lnum_t *occupied_cell_ids = NULL;
      cs_lnum_t *particle_list = NULL;

      if (   cs_glob_lagr_model->agglomeration == 1
          || cs_glob_lagr_model->fragmentation == 1 ) {

        n_occupied_cells
          = _get_n_occupied_cells(p_set, 0, p_set->n_particles);

        BFT_MALLOC(occupied_cell_ids, n_occupied_cells, cs_lnum_t);
        BFT_MALLOC(particle_list, n_occupied_cells+1, cs_lnum_t);

        _occupied_cells(p_set, 0, p_set->n_particles,
                        n_occupied_cells,
                        occupied_cell_ids,
                        particle_list);

      }

      /* Compute agglomeration and fragmentation
         (avoid second pass if second order scheme is used) */
      if (   cs_glob_lagr_time_step->nor == 1
          && ((cs_glob_lagr_model->agglomeration == 1) ||
              (cs_glob_lagr_model->fragmentation == 1))) {

        /* Initialize lists (ids of cells and particles) */
        cs_lnum_t *cell_particle_idx;

        BFT_MALLOC(cell_particle_idx, n_occupied_cells+1, cs_lnum_t);
        cell_particle_idx[0] = 0;

        cs_lnum_t enter_parts = p_set->n_particles;

        /* Loop on all cells that contain at least one particle */
        for (cs_lnum_t icell = 0; icell < n_occupied_cells; ++icell) {

          cs_lnum_t cell_id = occupied_cell_ids[icell];

          /* Particle indices: between start_part and end_part (list) */
          cs_lnum_t start_part = particle_list[icell];
          cs_lnum_t end_part = particle_list[icell+1];

          cs_lnum_t init_particles = p_set->n_particles;

          /* Treat agglomeration */
          if (cs_glob_lagr_model->agglomeration == 1) {

            cs_lagr_agglomeration(cell_id,
                                  dt[0],
                                  minimum_particle_diam,
                                  start_part,
                                  end_part);
          }

          /* Save number of created particles */

          cs_lnum_t inserted_parts_agglo = p_set->n_particles - init_particles;

          /* Create local buffer (deleted particles at the end) */

          cs_lnum_t local_size = end_part - start_part;
          cs_lnum_t deleted_parts = _get_n_deleted(p_set, start_part, end_part);
          size_t swap_buffer_size =   p_set->p_am->extents
                                    * (local_size - deleted_parts);
          size_t swap_buffer_deleted = p_set->p_am->extents * deleted_parts;

          /* Create buffers for deleted particles */
          unsigned char * swap_buffer, *deleted_buffer;
          BFT_MALLOC(swap_buffer, swap_buffer_size, unsigned char);
          BFT_MALLOC(deleted_buffer, swap_buffer_deleted, unsigned char);

          /* Update buffer for existing particles */
          cs_lnum_t count_del = 0, count_swap = 0;
          for (cs_lnum_t i = start_part; i < end_part; ++i) {
            if (cs_lagr_particles_get_flag(p_set, i,
                                           CS_LAGR_PART_TO_DELETE)) {
              memcpy(deleted_buffer + p_set->p_am->extents * count_del,
                     p_set->p_buffer + p_set->p_am->extents * i,
                     p_set->p_am->extents);
              count_del++;
            }
            else {
              memcpy(swap_buffer + p_set->p_am->extents * count_swap,
                     p_set->p_buffer + p_set->p_am->extents * i,
                     p_set->p_am->extents);
              count_swap++;
            }
          }

          memcpy(p_set->p_buffer + p_set->p_am->extents * start_part,
                 swap_buffer, swap_buffer_size);
          memcpy(  p_set->p_buffer
                 + p_set->p_am->extents * (local_size-deleted_parts+start_part),
                 deleted_buffer, swap_buffer_deleted);

          BFT_FREE(deleted_buffer);
          BFT_FREE(swap_buffer);

          /* Treat fragmentation */
          init_particles = p_set->n_particles;

          if (cs_glob_lagr_model->fragmentation == 1) {
            cs_lagr_fragmentation(dt[0],
                                  minimum_particle_diam,
                                  start_part,
                                  end_part - deleted_parts,
                                  init_particles,
                                  p_set->n_particles);
          }
          cs_lnum_t inserted_parts_frag = p_set->n_particles - init_particles;

          cell_particle_idx[icell+1] =   cell_particle_idx[icell]
                                     + inserted_parts_agglo + inserted_parts_frag;
        }

        p_set->n_particles = enter_parts;

        /* Introduce new particles (uniformly in the cell) */
        cs_lagr_new_v(p_set,
                      n_occupied_cells,
                      occupied_cell_ids,
                      cell_particle_idx);
        p_set->n_particles += cell_particle_idx[n_occupied_cells];

        BFT_FREE(cell_particle_idx);
      }

      BFT_FREE(occupied_cell_ids);
      BFT_FREE(particle_list);

      /* Reverse coupling: compute source terms
         -------------------------------------- */

      if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
          && cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order)
        cs_lagr_coupling(taup, tempct, tsfext, cpgd1, cpgd2, cpght);

      /* Deallocate arrays whose size is based on p_set->n_particles
         (which may change next) */

      BFT_FREE(tlag);
      BFT_FREE(taup);
      BFT_FREE(piil);
      BFT_FREE(bx);

      if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING)
        BFT_FREE(tsfext);

      if (beta != NULL)
        BFT_FREE(beta);

      if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
          && lagr_model->physical_model == CS_LAGR_PHYS_COAL
          && cs_glob_lagr_source_terms->ltsthe == 1) {
        BFT_FREE(cpgd1);
        BFT_FREE(cpgd2);
        BFT_FREE(cpght);
      }

      if (   (   lagr_model->physical_model == CS_LAGR_PHYS_HEAT
              && cs_glob_lagr_specific_physics->itpvar == 1)
          || lagr_model->physical_model == CS_LAGR_PHYS_COAL)
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

          if (cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_CELL_ID) >= 0) {
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

          if (cs_lagr_particles_get_flag(p_set, npt, deposition_mask)) {

            cs_real_t p_depo_time
              = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_DEPO_TIME);
            cs_real_t p_consol_height
              = CS_MIN(cs_lagr_particles_get_real(p_set, npt, CS_LAGR_HEIGHT),
                       cs_glob_lagr_consolidation_model->rate_consol * p_depo_time);
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

          if (cs_lagr_particles_get_flag(p_set, npt, deposition_mask)) {

            cs_lnum_t face_id = cs_lagr_particles_get_lnum
                                  (p_set, npt, CS_LAGR_NEIGHBOR_FACE_ID);

            cs_real_t p_diam = cs_lagr_particles_get_real
                                 (p_set, npt, CS_LAGR_DIAMETER);

            bound_stat[lag_bdi->iclogt * n_b_faces + face_id]
              += cs_lagr_particles_get_real(p_set, npt, CS_LAGR_DEPO_TIME);
            bound_stat[lag_bdi->iclogh * n_b_faces + face_id]
              +=  cs_lagr_particles_get_real(p_set, npt, CS_LAGR_CONSOL_HEIGHT)
                  * cs_math_pi * cs_math_sq(p_diam) * 0.25 / b_face_surf[face_id];

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
                       "** Ifac = %ld  and inclg = %g\n"
                       "-------------------------------------------------\n",
                       (long)ifac, bound_stat[lag_bdi->inclg * n_b_faces + ifac]);
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

  else if (cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt)
    cs_lagr_stat_update();

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

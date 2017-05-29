/*============================================================================
 * Methods for particle parameters
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
#include "cs_prototypes.h"

#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_physical_model.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_default.h"

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
     .nflagm = 100,
     .ndlaim = 10,
     .ncharm2 = 5,
     .nlayer = 5};

const cs_lagr_const_dim_t *cs_glob_lagr_const_dim
  = (const cs_lagr_const_dim_t *)(&_lagr_const_dim);

/* General dimensions */

cs_lagr_dim_t _lagr_dim = {.ntersl = 0,
                           .nvisbr = 0};

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
static cs_lagr_consolidation_model_t _cs_glob_lagr_consolidation_model = {0, 0, 0, 0};
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
     .itsvx = 0,
     .itsvy = 0,
     .itsvz = 0,
     .itsli = 0,
     .itske = 0,
     .itsr11 = 0,
     .itsr12 = 0,
     .itsr13 = 0,
     .itsr22 = 0,
     .itsr23 = 0,
     .itsr33 = 0,
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
     .inbrbd = 0,
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
     .iirayo = 0,
     .icp = -1,
     .diftl0 = 0,
     .cmu = 0,
     .visls0 = 0,
     .uetbor = NULL,
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
     .cvar_r33 = NULL};

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

cs_lagr_zone_class_data_t *_lagr_zone_class_data = NULL;

cs_lagr_bdy_condition_t  *cs_glob_lagr_bdy_conditions = NULL;

cs_lagr_internal_condition_t  *cs_glob_lagr_internal_conditions = NULL;

/* Array dimensions for lagrangien user data per zone per class */

int  cs_glob_lagr_nzone_max  = 0;
int  cs_glob_lagr_nclass_max = 0;

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
cs_f_lagr_reentrained_model_pointers(cs_int_t         **ireent,
                                     cs_int_t         **iflow,
                                     cs_real_t        **espasg,
                                     cs_real_t        **denasp,
                                     cs_real_t        **modyeq,
                                     cs_real_t        **rayasp,
                                     cs_real_t        **rayasg);

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
                                cs_int_t **p_itsvx,
                                cs_int_t **p_itsvy,
                                cs_int_t **p_itsvz,
                                cs_int_t **p_itsli,
                                cs_int_t **p_itske,
                                cs_int_t **p_itsr11,
                                cs_int_t **p_itsr12,
                                cs_int_t **p_itsr13,
                                cs_int_t **p_itsr22,
                                cs_int_t **p_itsr23,
                                cs_int_t **p_itsr33,
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
cs_f_lagr_reentrained_model_pointers(cs_int_t         **ireent,
                                     cs_int_t         **iflow,
                                     cs_real_t        **espasg,
                                     cs_real_t        **denasp,
                                     cs_real_t        **modyeq,
                                     cs_real_t        **rayasp,
                                     cs_real_t        **rayasg)
{
  *ireent = &cs_glob_lagr_reentrained_model->ireent;
  *iflow  = &cs_glob_lagr_reentrained_model->iflow;
  *espasg = &cs_glob_lagr_reentrained_model->espasg;
  *denasp = &cs_glob_lagr_reentrained_model->denasp;
  *modyeq = &cs_glob_lagr_reentrained_model->modyeq;
  *rayasp = &cs_glob_lagr_reentrained_model->rayasp;
  *rayasg = &cs_glob_lagr_reentrained_model->rayasg;
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
                                cs_int_t **p_itsvx,
                                cs_int_t **p_itsvy,
                                cs_int_t **p_itsvz,
                                cs_int_t **p_itsli,
                                cs_int_t **p_itske,
                                cs_int_t **p_itsr11,
                                cs_int_t **p_itsr12,
                                cs_int_t **p_itsr13,
                                cs_int_t **p_itsr22,
                                cs_int_t **p_itsr23,
                                cs_int_t **p_itsr33,
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
  *p_itsvx  = &cs_glob_lagr_source_terms->itsvx;
  *p_itsvy  = &cs_glob_lagr_source_terms->itsvy;
  *p_itsvz  = &cs_glob_lagr_source_terms->itsvz;
  *p_itsli  = &cs_glob_lagr_source_terms->itsli;
  *p_itske  = &cs_glob_lagr_source_terms->itske;
  *p_itsr11 = &cs_glob_lagr_source_terms->itsr11;
  *p_itsr12 = &cs_glob_lagr_source_terms->itsr12;
  *p_itsr13 = &cs_glob_lagr_source_terms->itsr13;
  *p_itsr22 = &cs_glob_lagr_source_terms->itsr22;
  *p_itsr23 = &cs_glob_lagr_source_terms->itsr23;
  *p_itsr33 = &cs_glob_lagr_source_terms->itsr33;
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

  _lagr_extra_module.iirayo = *iirayo;
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

  if (cs_field_by_name_try("pressure") != NULL) {
    /* we use Code_Saturne */
    _lagr_extra_module.pressure    = cs_field_by_name_try("pressure");
    _lagr_extra_module.vel         = cs_field_by_name_try("velocity");
    _lagr_extra_module.cvar_k      = cs_field_by_name_try("k");
    _lagr_extra_module.cvar_ep     = cs_field_by_name_try("epsilon");
    _lagr_extra_module.cvar_omg    = cs_field_by_name_try("omega");
    _lagr_extra_module.cvar_r11    = cs_field_by_name_try("r11");
    _lagr_extra_module.cvar_r22    = cs_field_by_name_try("r22");
    _lagr_extra_module.cvar_r33    = cs_field_by_name_try("r33");
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
                    cs_field_key_id("scalar_diffusivity_ref"));

        int l_id = cs_field_get_key_int(_lagr_extra_module.scal_t,
                cs_field_key_id("scalar_diffusivity_id"));
        if (l_id >= 0)
            _lagr_extra_module.cpro_viscls = cs_field_by_id(l_id);
    }

    _lagr_extra_module.cpro_cp     = cs_field_by_name_try("specific_heat");
    _lagr_extra_module.temperature = cs_field_by_name_try("temperature");
    _lagr_extra_module.t_gaz       = cs_field_by_name_try("t_gas");
    _lagr_extra_module.luminance   = cs_field_by_name_try("luminance");
    _lagr_extra_module.x_oxyd      = cs_field_by_name_try("ym_o2");
    _lagr_extra_module.x_eau       = cs_field_by_name_try("ym_h2o");
    _lagr_extra_module.x_m         = cs_field_by_name_try("xm");

    cs_field_t *f = cs_field_by_name_try("ustar");
    if (f != NULL)
        _lagr_extra_module.uetbor = f->val;
    else
        _lagr_extra_module.uetbor = NULL;
  }
  else {
    /* we use NEPTUNE_CFD */
    _lagr_extra_module.pressure    = cs_field_by_name_try("Pressure");
    _lagr_extra_module.vel         = cs_field_by_name_try("lagr_velocity");
    _lagr_extra_module.cvar_k      = cs_field_by_name_try("lagr_k");
    _lagr_extra_module.cvar_ep     = cs_field_by_name_try("lagr_epsilon");
    _lagr_extra_module.cvar_omg    = NULL;
    _lagr_extra_module.cvar_r11    = cs_field_by_name_try("lagr_r11");
    _lagr_extra_module.cvar_r22    = cs_field_by_name_try("lagr_r22");
    _lagr_extra_module.cvar_r33    = cs_field_by_name_try("lagr_r33");
    _lagr_extra_module.viscl       = cs_field_by_name_try("lagr_molecular_viscosity");
    _lagr_extra_module.scal_t      = cs_field_by_name_try("lagr_enthalpy");
    _lagr_extra_module.cpro_viscls = cs_field_by_name_try("lagr_thermal_conductivity");
    _lagr_extra_module.cpro_cp     = cs_field_by_name_try("lagr_specific_heat");
    _lagr_extra_module.temperature = cs_field_by_name_try("lagr_temperature");
    _lagr_extra_module.t_gaz       = NULL;
    _lagr_extra_module.luminance   = cs_field_by_name_try("luminance");
    _lagr_extra_module.x_oxyd      = NULL;
    _lagr_extra_module.x_eau       = NULL;
    _lagr_extra_module.x_m         = NULL;
    _lagr_extra_module.cromf       = cs_field_by_name_try("lagr_density");
    /* TODO FIX ME */
    _lagr_extra_module.visls0      = 0.;

    cs_field_t *f = cs_field_by_name_try("wall_friction_velocity");
    if (f != NULL)
        _lagr_extra_module.uetbor  = f->val;
    else
        _lagr_extra_module.uetbor  = NULL;
    }
}

static void
_cs_lagr_free_zone_class_data(cs_lagr_zone_class_data_t *zone_class_data)
{
  assert(zone_class_data != NULL);

  if (cs_glob_lagr_model->physical_model == 1)
    BFT_FREE(zone_class_data->temperature);

  else if (cs_glob_lagr_model->physical_model == 2) {

    BFT_FREE(zone_class_data->coke_density);
    BFT_FREE(zone_class_data->temperature);
    BFT_FREE(zone_class_data->coal_mass_fraction);
    BFT_FREE(zone_class_data->coke_mass_fraction);

  }
}

static void
_cs_lagr_free_all_zone_class_data(void)
{
  if (_lagr_zone_class_data != NULL) {
    for (int  i = 0; i < cs_glob_lagr_nzone_max * cs_glob_lagr_nclass_max ; i++)
      _cs_lagr_free_zone_class_data(&(_lagr_zone_class_data[i]));
    BFT_FREE(_lagr_zone_class_data);
  }
}

/*----------------------------------------------------------------------------
 * Zone and class data structure allocation for a given class and zone.
 * Reallocation of the main array if necessary.
 *----------------------------------------------------------------------------*/

static cs_lagr_zone_class_data_t *
_cs_lagr_allocate_zone_class_data(int  iclass,
                                  int  izone)
{
  if (izone >= cs_glob_lagr_nzone_max || iclass >= cs_glob_lagr_nclass_max) {

    int  old_lagr_nzone  = cs_glob_lagr_nzone_max;
    int  old_lagr_nclass = cs_glob_lagr_nclass_max;

    if (izone >= cs_glob_lagr_nzone_max)
      cs_glob_lagr_nzone_max  = CS_MAX(izone +1, cs_glob_lagr_nzone_max  + 5);
    if (iclass >= cs_glob_lagr_nclass_max)
      cs_glob_lagr_nclass_max = CS_MAX(iclass +1, cs_glob_lagr_nclass_max + 1);

    BFT_REALLOC(_lagr_zone_class_data,
                cs_glob_lagr_nzone_max * cs_glob_lagr_nclass_max,
                cs_lagr_zone_class_data_t);

    if (cs_glob_lagr_nzone_max != old_lagr_nzone) {
      for (int  ii = old_lagr_nclass-1; ii > 0; ii--) { /* no-op for ii = 0 */
        for (int  jj = old_lagr_nzone-1; jj >= 0; jj--) {
          _lagr_zone_class_data[cs_glob_lagr_nzone_max * ii + jj]
            = _lagr_zone_class_data[old_lagr_nzone * ii + jj];
          memset(_lagr_zone_class_data + (old_lagr_nzone * ii + jj),
                 0,
                 sizeof(cs_lagr_zone_class_data_t));
        }
      }
    }

    for (int ii = old_lagr_nclass; ii < cs_glob_lagr_nclass_max; ii++)
      for (int jj = 0; jj < cs_glob_lagr_nzone_max; jj++) {
        memset(_lagr_zone_class_data + (cs_glob_lagr_nzone_max * ii + jj),
               0,
               sizeof(cs_lagr_zone_class_data_t));
    }
    for (int ii = 0; ii < cs_glob_lagr_nclass_max; ii++)
      for (int jj = old_lagr_nzone; jj < cs_glob_lagr_nzone_max; jj++) {
        memset(_lagr_zone_class_data + (cs_glob_lagr_nzone_max * ii + jj),
               0,
               sizeof(cs_lagr_zone_class_data_t));
    }
  }

  cs_lagr_zone_class_data_t *zone_class_data
    = cs_lagr_get_zone_class_data(iclass, izone);

  assert(zone_class_data != NULL);

  /* On first call for this zone */

  if (zone_class_data->nb_part == 0) {

    if (cs_glob_lagr_model->physical_model == 1) {
      BFT_MALLOC(zone_class_data->temperature,
                 1,
                 cs_real_t);
    }
    else if (cs_glob_lagr_model->physical_model == 2) {

      BFT_MALLOC(zone_class_data->coke_density,
                 cs_glob_lagr_model->n_temperature_layers,
                 cs_real_t);
      BFT_MALLOC(zone_class_data->temperature,
                 cs_glob_lagr_model->n_temperature_layers,
                 cs_real_t);
      BFT_MALLOC(zone_class_data->coal_mass_fraction,
                 cs_glob_lagr_model->n_temperature_layers,
                 cs_real_t);
      BFT_MALLOC(zone_class_data->coke_mass_fraction,
                 cs_glob_lagr_model->n_temperature_layers,
                 cs_real_t);

    }

  }

  return zone_class_data;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a cs_lagr_bdy_condition_t structure.
 *
 * parameters:
 *   n_max_zones     <--  number max. of boundary zones
 *
 * returns:
 *   a new defined cs_lagr_bdy_condition_t structure
 *----------------------------------------------------------------------------*/

static cs_lagr_bdy_condition_t *
_create_bdy_cond_struct(int   n_max_zones)
{
  int   i;

  cs_lagr_bdy_condition_t *bdy_cond = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  BFT_MALLOC(bdy_cond, 1, cs_lagr_bdy_condition_t);

  bdy_cond->n_b_zones = 0;
  bdy_cond->n_b_max_zones = n_max_zones;

  BFT_MALLOC(bdy_cond->particle_flow_rate, n_max_zones, cs_real_t);
  BFT_MALLOC(bdy_cond->b_zone_id, n_max_zones, int);
  BFT_MALLOC(bdy_cond->b_zone_classes, n_max_zones, int);
  BFT_MALLOC(bdy_cond->b_zone_natures, n_max_zones, int);

  for (i = 0; i < n_max_zones; i++) {

    bdy_cond->particle_flow_rate[i] = 0.0;
    bdy_cond->b_zone_id[i] = -1;
    bdy_cond->b_zone_classes[i] = 0;
    bdy_cond->b_zone_natures[i] = -1;

  }

  BFT_MALLOC(bdy_cond->b_face_zone_id, mesh->n_b_faces, int);

  for (i = 0; i < cs_glob_mesh->n_b_faces; i++)
    bdy_cond->b_face_zone_id[i] = -1;

  bdy_cond->steady_bndy_conditions = false;

  return bdy_cond;
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

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

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
  int   nvisbr = cs_glob_lagr_dim->nvisbr;

  assert(bound_stat == NULL);

  if (nvisbr > 0)
    BFT_MALLOC(bound_stat, n_b_faces * nvisbr, cs_real_t);

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
  int  nvisbr = cs_glob_lagr_dim->nvisbr;

  if (nvisbr > 0) {
    assert(bound_stat != NULL);
    BFT_FREE(bound_stat);
  }

  BFT_FREE(cs_glob_lagr_precipitation_model->nbprec);
  BFT_FREE(cs_glob_lagr_precipitation_model->solub);

  BFT_FREE(cs_glob_lagr_precipitation_model->mp_diss);

  BFT_FREE(cs_glob_lagr_source_terms->st_val);

  _cs_lagr_free_all_zone_class_data();

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

  for (int i = 0; i < cs_glob_lagr_dim->nvisbr; i++)
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set temperature parameters for a given class and boundary zone
 *
 * \param[in]   iclass      class number
 * \param[in]   izone       boundary zone number
 * \param[in]   profile     temperature profile
 * \param[in]   temp        pointer to temperature values
 * \param[in]   emissivity  emissivity value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_temperature(int         iclass,
                                   int         izone,
                                   int         profile,
                                   cs_real_t  *temp,
                                   cs_real_t   emissivity)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  assert(cs_glob_lagr_model->physical_model == 1);

  /* temperature, emissivity */
  zonedata->temperature_profile = profile;

  zonedata->temperature[0] = temp[0];
  zonedata->emissivity = emissivity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set temperature parameters for a given class and boundary zone
 *
 * \param[in]   iclass  class number
 * \param[in]   izone   boundary zone number
 * \param[in]   cp      pointer to specific heat value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_cp(int        iclass,
                          int        izone,
                          cs_real_t  cp)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->cp = cp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set coal parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    coal profile
 * \param[in]   number     coal number
 * \param[in]   temp       pointer to temperature array
 * \param[in]   coal_mf    pointer to coal mass fraction
 * \param[in]   coke_mf    pointer to coke mass fraction
 * \param[in]   coke_density  pointer to coke density after pyrolysis
 * \param[in]   water_mf   pointer to water mass fraction
 * \param[in]   shrink_diam  pointer to coke shrinking diameter
 * \param[in]   init_diam  pointer to initial particle diameter
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_coal(int         iclass,
                            int         izone,
                            int         profile,
                            int         number,
                            cs_real_t  *temp,
                            cs_real_t  *coal_mf,
                            cs_real_t  *coke_mf,
                            cs_real_t  *coke_density,
                            cs_real_t   water_mf,
                            cs_real_t   shrink_diam,
                            cs_real_t   init_diam)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->coal_profile = profile;
  zonedata->coal_number  = number;

  for (cs_lnum_t ilayer = 0;
       ilayer < cs_glob_lagr_model->n_temperature_layers;
       ilayer++) {

    if (temp != NULL)
      zonedata->temperature[ilayer] = temp[ilayer];

    if (coke_density != NULL)
      zonedata->coke_density[ilayer] = coke_density[ilayer];

    if (coal_mf != NULL)
      zonedata->coal_mass_fraction[ilayer] = coal_mf[ilayer];

    if (coke_mf !=NULL)
      zonedata->coke_mass_fraction[ilayer] = coke_mf[ilayer];

  }

  zonedata->initial_diameter   = init_diam;
  zonedata->shrinking_diameter = shrink_diam;

  zonedata->water_mass_fraction = water_mf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set coal parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    pointer to flag for flow and stat weight profile
 * \param[in]   weight     pointer to stat weight value
 * \param[in]   flow       pointer to mass flow rate value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_stat(int        iclass,
                            int        izone,
                            int        profile,
                            cs_real_t  weight,
                            cs_real_t  flow)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->distribution_profile = profile;
  zonedata->stat_weight          = weight;
  zonedata->flow_rate            = flow;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set diameter parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    pointer to flag for diameter profile
 * \param[in]   diam       pointer to diameter value
 * \param[in]   diam_dev   pointer to diameter standard deviation value
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_diam (int        iclass,
                             int        izone,
                             int        profile,
                             cs_real_t  diam,
                             cs_real_t  diam_dev)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->diameter_profile  = profile;
  zonedata->diameter          = diam;
  zonedata->diameter_variance = diam_dev;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set injection parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   number     pointer to number of particles to inject
 * \param[in]   freq       pointer to injection frequency
 * \param[in]   stat       pointer to statistical groups id
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_injection(int        iclass,
                                 int        izone,
                                 int        number,
                                 int        freq,
                                 int        stat)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->nb_part = number;
  zonedata->injection_frequency = freq;
  zonedata->cluster = stat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set velocity parameters for a given class and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   profile    pointer to velocity profile
 * \param[in]   velocity   pointer to velocity values array
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_velocity(int        iclass,
                                int        izone,
                                int        profile,
                                cs_real_t  velocity[])
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->velocity_profile = profile;

  if (zonedata->velocity_profile == 0)
    zonedata->velocity_magnitude = velocity[0];

  else if (zonedata->velocity_profile == 1) {

    for (int  i = 0; i < 3; i++)
      zonedata->velocity[i] = velocity[i];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set density for a given class of particls and boundary zone
 *
 * \param[in]   iclass     class number
 * \param[in]   izone      boundary zone number
 * \param[in]   density    pointer to density value
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_density(int        iclass,
                               int        izone,
                               cs_real_t  density)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->density = density;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set fouling index for a given class of particles and boundary zone
 *
 * \param[in]   iclass      class number
 * \param[in]   izone       boundary zone number
 * \param[in]   foul_index  pointer to fouling index value
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_set_zone_class_foul_index(int        iclass,
                                  int        izone,
                                  cs_real_t  foul_index)
{
  cs_lagr_zone_class_data_t *zonedata
    = cs_lagr_get_zone_class_data(iclass, izone);

  zonedata->foul_index = foul_index;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to class/boundary zone parameters structure
 *
 * \param[in]    iclass     particle class number
 * \param[in]    izone      boundary zone number
 *
 * \return
 *   pointer to particle class and boundary zone structure of parameters
 */
/*----------------------------------------------------------------------------*/

cs_lagr_zone_class_data_t *
cs_lagr_get_zone_class_data(int   iclass,
                            int   izone)
{
  assert(_lagr_zone_class_data != NULL);

  return &(_lagr_zone_class_data[iclass*cs_glob_lagr_nzone_max + izone]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a new class/boundary zone parameters structure
 *
 * \param[in]    iclass     particle class number
 * \param[in]    izone      boundary zone number
 *
 * \return
 *   pointer to particle class and boundary zone structure of parameters
 */
/* ----------------------------------------------------------------------------*/

cs_lagr_zone_class_data_t *
cs_lagr_init_zone_class_new(int       iclass,
                            int       izone)
{
  cs_lagr_zone_class_data_t *zonedata
    = _cs_lagr_allocate_zone_class_data(iclass, izone);

  zonedata->nb_part              =  0;
  zonedata->injection_frequency  =  0;
  zonedata->velocity_profile     = -2;
  zonedata->distribution_profile = -2;
  zonedata->temperature_profile  = -2;
  zonedata->diameter_profile     = -2;

  if (cs_glob_lagr_model->physical_model == 2) {

    zonedata->coal_profile       = -2;
    zonedata->coal_number        = -2;

  }

  zonedata->cluster              =  0;

  zonedata->velocity_magnitude   = - cs_math_big_r;

  for (int  i = 0; i < 3; i++)
    zonedata->velocity[i]        = - cs_math_big_r;

  zonedata->stat_weight          = - cs_math_big_r;
  zonedata->diameter             = - cs_math_big_r;
  zonedata->diameter_variance    = - cs_math_big_r;
  zonedata->density              = - cs_math_big_r;

  if (cs_glob_lagr_model->physical_model == 1) {

    if (cs_glob_lagr_specific_physics->itpvar == 1 ) {

      zonedata->temperature[0]   = - cs_math_big_r;
      zonedata->cp               = - cs_math_big_r;
      zonedata->emissivity       = - cs_math_big_r;

    }

  }

  else if (cs_glob_lagr_model->physical_model == 2) {

    zonedata->cp                 = - cs_math_big_r;

    for (cs_lnum_t ilayer = 0;
         ilayer < cs_glob_lagr_model->n_temperature_layers;
         ilayer ++) {

      zonedata->temperature[ilayer]        = - cs_math_big_r;
      zonedata->coal_mass_fraction[ilayer] = - cs_math_big_r;
      zonedata->coke_mass_fraction[ilayer] = - cs_math_big_r;
      zonedata->coke_density[ilayer]       = - cs_math_big_r;

    }

    zonedata->water_mass_fraction = - cs_math_big_r;
    zonedata->shrinking_diameter  = - cs_math_big_r;
    zonedata->initial_diameter    = - cs_math_big_r;

  }

  zonedata->flow_rate     = 0.0;

  return zonedata;
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
    BFT_MALLOC(cs_glob_lagr_internal_conditions->i_face_zone_id, cs_glob_mesh->n_i_faces, int);

    for (cs_lnum_t i = 0; i < cs_glob_mesh->n_i_faces; i++)
      cs_glob_lagr_internal_conditions->i_face_zone_id[i] = -1;

  }

  return cs_glob_lagr_internal_conditions;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to the main boundary conditions structure.
 *
 * \return
 *   pointer to current bdy_contditions or NULL
 */
/*----------------------------------------------------------------------------*/

cs_lagr_bdy_condition_t  *
cs_lagr_get_bdy_conditions(void)
{
  /* Define a structure with default parameters if not done yet */

  if (cs_glob_lagr_bdy_conditions == NULL)
    cs_glob_lagr_bdy_conditions
      = _create_bdy_cond_struct(cs_glob_lagr_const_dim->nflagm);

  return cs_glob_lagr_bdy_conditions;
}

/*----------------------------------------------------------------------------
 * Destroy finalize the global cs_lagr_bdy_condition_t structure.
 *----------------------------------------------------------------------------*/

void
cs_lagr_finalize_bdy_cond(void)
{
  cs_lagr_bdy_condition_t  *bdy_cond = cs_glob_lagr_bdy_conditions;

  if (bdy_cond != NULL) {

    BFT_FREE(bdy_cond->b_zone_id);
    BFT_FREE(bdy_cond->b_zone_natures);
    BFT_FREE(bdy_cond->b_zone_classes);

    BFT_FREE(bdy_cond->b_face_zone_id);

    BFT_FREE(bdy_cond->particle_flow_rate);

    BFT_FREE(cs_glob_lagr_bdy_conditions);

  }
}

/*----------------------------------------------------------------------------
 * Destroy finalize the global cs_lagr_internal_condition_t structure.
 *----------------------------------------------------------------------------*/

void cs_lagr_finalize_internal_cond(void)
{
  cs_lagr_internal_condition_t  *internal_cond = cs_glob_lagr_internal_conditions;
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

  /* ====================================================================   */
  /* 1.  INITIALISATIONS  */
  /* ====================================================================   */

  /* ->Sur Champ fige Lagrangien :
   *   values at previous time step = values at current time step
   *   Rem : cette boucle pourrait etre faite au 1er passage
   *         mais la presence de cs_user_extra_operations incite a la prudence...*/

  if (cs_glob_lagr_time_scheme->iilagr == 3) {

    int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++){

      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE)
        cs_field_current_to_previous(f);

    }

  }

  cs_lagr_particle_set_t *p_set = cs_glob_lagr_particle_set;

  /* First call (initializations)
     ---------------------------- */

  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_prev + 1) {

    /* --> if the deposition model is activated */

    if (lagr_model->deposition >= 1) {

      cs_real_t ustarmoy = 0.0;
      cs_real_t surftot  = 0.0;
      cs_real_3_t  dtmp;
      dtmp[0]  = 0.0;

      /* boundary faces data  */

      cs_lagr_geom();

      /* Average friction velocity calculation    */

      for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

        if (itypfb[ifac] == CS_SMOOTHWALL || itypfb[ifac] == CS_ROUGHWALL) {

          cs_lnum_t iel   = ifabor[ifac];

          /* the density pointer according to the flow location */
          cs_real_t romf   = extra->cromf->val[iel];
          cs_real_t visccf = extra->viscl->val[iel] / romf;

          cs_real_t _ustar = CS_MAX(extra->uetbor[ifac], 1e-15);

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
                     "---------------------------------------------\n\n\n"
                     "** Mean friction velocity  (ustar) =  %7.3F\n"
                     "---------------------------------------------------------------\n"),
                   ustarmoy);

    }

  }

  /* Initialization  */

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

  /* ====================================================================   */
  /* 1.bis  Initialization for the dlvo, roughness and clogging  model */
  /* ====================================================================   */

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

  /* ====================================================================   */
  /*   Initialization for the dlvo model */
  /* ====================================================================   */

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

  /* ====================================================================   */
  /*  Initialization for the roughness surface model    */
  /* ====================================================================   */

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

  /* ====================================================================   */
  /*   Initialization for the clogging model  */
  /* ====================================================================   */

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

  /* ====================================================================   */
  /* 2.  MISE A JOUR DES NOUVELLES PARTICULES ENTREES DANS LE DOMAINE  */
  /* ====================================================================   */
  /* At the first time step we initialize particles to  */
  /* values at current time step and not at previous time step, because     */
  /* values at previous time step = initialization */

  int  iprev;
  if (ts->nt_cur == 1) /* Use fields at current time step     */
    iprev = 0;

  else  /* Use fields at previous time step    */
    iprev = 1;

  cs_lagr_injection(iprev, itypfb, vislen);

  /* ====================================================================   */
  /* 3.  GESTION DU TEMPS QUI PASSE...   */
  /* ====================================================================   */

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

    /* ====================================================================   */
    /* 4.  GRADIENT DE PRESSION ET DE LA VITESSE FLUIDE   */
    /* ====================================================================   */
    /* At the first time step we initialize particles to  */
    /* values at current time step and not at previous time step, because     */
    /* values at previous time step = initialization (null gradients)    */

    mode = 1;
    if (ts->nt_cur == 1)
      mode = 0;

    cs_lagr_gradients(mode, gradpr, gradvf);

    /* ====================================================================   */
    /* 5. PROGRESSION DES PARTICULES  */
    /* ====================================================================   */

    bool go_on = true;
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

      /* --> Recopie des resultats de l'etape precedente :  */

      if (cs_glob_lagr_time_step->nor == 1) {
        for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++)
          cs_lagr_particles_current_to_previous(p_set, ip);
      }

      /* ----> COMPUTATION OF THE FLUID'S PRESSURE AND VELOCITY GRADIENT   */
      /*       AT N+1 (with values at current time step)    */
      if (   cs_glob_lagr_time_step->nor == 2
          && cs_glob_lagr_time_scheme->iilagr != 3)
        cs_lagr_gradients(0, gradpr, gradvf);

      /* use fields at previous or current time step */
      if (cs_glob_lagr_time_step->nor == 1)
        /* Use fields at previous time step    */
        iprev     = 1;

      else
        iprev     = 0;

      /* Retrieve bx values associated with particles from previous pass   */

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

      /* --> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES  */
      /*     POSITION, VITESSE FLUIDE, VITESSE PARTICULE    */

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

      /* --> INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES  */
      /*     LIEES AUX PHYSIQUES PARTICULIERES PARTICULAIRES     */

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

      /* ====================================================================   */
      /* 8.  Reperage des particules - Traitement des conditions aux limites    */
      /*     pour la position des particules */
      /* ====================================================================   */

      if (cs_glob_lagr_time_step->nor == 1) {

        /* -> Si on est en instationnaire, RAZ des statistiques aux frontieres    */

        if (   cs_glob_lagr_time_scheme->isttio == 0
            || (   cs_glob_lagr_time_scheme->isttio == 1
                &&    cs_glob_time_step->nt_cur
                   <= cs_glob_lagr_stat_options->nstist)) {

          lag_bdi->tstatp = 0.0;
          lag_bdi->npstf  = 0;

          for (int  ii = 0; ii < cs_glob_lagr_dim->nvisbr; ii++) {

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

      /* ====================================================================   */
      /* 11.  CALCUL DE L'ADHESION SI MODELE DE REENTRAINEMENT   */
      /* ====================================================================   */

      if (lagr_model->resuspension > 0)
        cs_lagr_resuspension();

      /* ====================================================================   */
      /* 11.  CALCUL STATISTIQUES  */
      /* ====================================================================   */

      /* Calculation of consolidation:
       * linear increase of the consolidation height with deposit time */
      if (cs_glob_lagr_consolidation_model->iconsol > 0) {

        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

          if (cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_DEPOSITION_FLAG) > 0) {

            cs_real_t p_depo_time = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_DEPO_TIME);
            cs_real_t p_consol_height =
              CS_MIN( cs_lagr_particles_get_real(p_set, npt, CS_LAGR_HEIGHT),
                   cs_glob_lagr_consolidation_model->rate_consol * p_depo_time );
            cs_lagr_particles_set_real(p_set, npt, CS_LAGR_CONSOL_HEIGHT, p_consol_height);
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

      /*  STATISTICS FOR CLOGGING  */
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

            bft_printf("   ** LAGRANGIAN MODULE:\n"
                       "   ** Error in cs_lagr.c: inclg < 0 ! \n"
                       "---------------------------------------------\n\n\n"
                       "** Ifac = %d  and inclg = %g\n"
                       "-------------------------------------------------\n",
                       ifac, bound_stat[lag_bdi->inclg * n_b_faces + ifac]);
          }

        }

      }

      /* ====================================================================   */
      /* 12.  Equation de Poisson  */
      /* ====================================================================   */

      if (   cs_glob_lagr_time_step->nor == cs_glob_lagr_time_scheme->t_order
          && cs_glob_lagr_time_scheme->ilapoi == 1)
        cs_lagr_poisson(itypfb);

      /* ====================================================================   */
      /* 14. UN AUTRE TOUR ?  */
      /* ====================================================================   */

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

  /* ====================================================================
   * 18. ECRITURE SUR FICHIERS DES INFORMATIONS SUR LE NOMBRE DE PARTICULES
   *   - nombre de particules dans le domaine
   *   - nombre de particules entrantes
   *   - nombre de particules sorties
   *   - ...
   * ==================================================================== */

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

  /* Free memory     */
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

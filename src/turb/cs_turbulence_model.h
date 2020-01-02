#ifndef __CS_TURBULENCE_MODEL_H__
#define __CS_TURBULENCE_MODEL_H__

/*============================================================================
 * Base turbulence model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * turbulence models
 *----------------------------------------------------------------------------*/

enum {
  CS_TURB_NONE = 0,
  CS_TURB_MIXING_LENGTH = 10,
  CS_TURB_K_EPSILON = 20,
  CS_TURB_K_EPSILON_LIN_PROD = 21,
  CS_TURB_K_EPSILON_LS = 22,
  CS_TURB_K_EPSILON_QUAD = 23,
  CS_TURB_RIJ_EPSILON_LRR = 30,
  CS_TURB_RIJ_EPSILON_SSG = 31,
  CS_TURB_RIJ_EPSILON_EBRSM = 32,
  CS_TURB_LES_SMAGO_CONST = 40,
  CS_TURB_LES_SMAGO_DYN = 41,
  CS_TURB_LES_WALE = 42,
  CS_TURB_V2F_PHI = 50,
  CS_TURB_V2F_BL_V2K = 51,
  CS_TURB_K_OMEGA = 60,
  CS_TURB_SPALART_ALLMARAS = 70
};

/*----------------------------------------------------------------------------
 * turbulence type of model
 *----------------------------------------------------------------------------*/

enum {
/* We also use
  CS_TURB_NONE = 0, */
  CS_TURB_RANS = 1,
  CS_TURB_LES = 2,
  CS_TURB_HYBRID = 3
};

/*----------------------------------------------------------------------------
 * turbulence order of model
 *----------------------------------------------------------------------------*/

enum {
  CS_TURB_ALGEBRAIC = 0,
  CS_TURB_FIRST_ORDER = 1,
  CS_TURB_SECOND_ORDER = 2
};

/*----------------------------------------------------------------------------
 * hybrid models
 *----------------------------------------------------------------------------*/

enum {
  CS_HYBRID_NONE = 0,
  CS_HYBRID_DES  = 1,
  CS_HYBRID_DDES = 2,
  CS_HYBRID_SAS  = 3
};

/* turbulence model general options descriptor */
/*---------------------------------------------*/

typedef struct {

  int           iturb; /* turbulence model
                          CS_TURB_NONE: no turbulence model (laminar flow)
                          CS_TURB_MIXING_LENGTH: mixing length model
                          CS_TURB_K_EPSILON: standard k-epsilon model
                          CS_TURB_K_EPSILON_LIN_PROD: k-epsilon model with
                            Linear Production (LP) correction
                          CS_TURB_K_EPSILON_LS: Launder-Sharma low Re
                            k-epsilon model
                          CS_TURB_K_EPSILON_QUAD: Baglietto et al. low Re
                            k epsilon model
                          CS_TURB_RIJ_EPSILON_LRR: Rij-epsilon (LRR)
                          CS_TURB_RIJ_EPSILON_SSG: Rij-epsilon (SSG)
                          CS_TURB_RIJ_EPSILON_EBRSM: Rij-epsilon (EBRSM)
                          CS_TURB_LES_SMAGO_CONST: LES
                            (constant Smagorinsky model)
                          CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic
                            Smagorisky model)
                          CS_TURB_LES_WALE: LES (WALE)
                          CS_TURB_V2F_PHI: v2f phi-model
                          CS_TURB_V2F_BL_V2K: v2f BL-v2-k
                          CS_TURB_K_OMEGA: k-omega SST
                          CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model */
  int           itytur;       /* class of turbulence model (integer value
                                 iturb/10) */
  int           hybrid_turb;  /* Type of Hybrid Turbulence Model
                                   - CS_HYBRID_NONE : No model
                                   - CS_HYBRID_DES  : Detached Eddy Simulation
                                   - CS_HYBRID_DDES : Delayed Detached Eddy Simulation
                                   - CS_HYBRID_SAM  : Scale Adaptive Model */
  int           type;  /* Type of turbulence modelling:
                          - CS_TURB_NONE: No model
                          - CS_TURB_RANS: RANS modelling
                          - CS_TURB_LES: LES modelling
                          - CS_TURB_HYBRID: RANS -- LES modelling */
  int           order; /* Order of the turbulence model:
                          - CS_TURB_ALGEBRAIC: 0th order algebraik model
                          - CS_TURB_FIRST_ORDER: 1st order Eddy Viscosity type models
                          - CS_TURB_SECOND_ORDER: 2nd order Differential Reynolds Stress type models */
} cs_turb_model_t;


/* Reference values for turbulence structure and associated pointer */
/*------------------------------------------------------------------*/

typedef struct {
  double        almax;        /* characteristic macroscopic length of the
                                 domain */
  double        uref;         /* characteristic flow velocity */
} cs_turb_ref_values_t;


/* rans turbulence model descriptor */
/*----------------------------------*/

typedef struct {

  int           irccor;       /* activation of rotation/curvature correction for
                                 an eddy viscosity turbulence models
                                 - 0: false
                                 - 1: true */
  int           itycor;       /* type of rotation/curvature correction for an
                                 eddy viscosity turbulence models
                                 - 1: Cazalbou correction (default when irccor=1
                                      and itytur=2 or 5)
                                 - 2: Spalart-Shur correction (default when
                                      irccor=1 and iturb=60 or 70) */
  int           idirsm;       /* turbulent diffusion model for second moment
                                 closure
                                 - 0: scalar diffusivity (Shir model)
                                 - 1: tensorial diffusivity (Daly and Harlow
                                      model, default model) */
  int           iclkep;       /* clipping of k and epsilon
                                 - 0: absolute value clipping
                                 - 1: coupled clipping based on physical
                                      relationships */
  int           igrhok;       /* take (2/3 rho grad k) in the momentum
                                 equation
                                 - 1: true
                                 - 0: false (default) */
  int           igrake;       /* buoyant term in k-epsilon
                                 - 1: true (default if rho is variable)
                                 - 0: false */
  int           igrari;       /* buoyant term in Rij-epsilon
                                 - 1: true (default if rho is variable)
                                 - 0: false */
  int           ikecou;       /* partially coupled version of
                                 k-epsilon (only for iturb=20)
                                 - 1: true (default)
                                 - 0: false */
  int           reinit_turb;  /* Advanced re-init for EBRSM and k-omega models
                                 - 1: true (default)
                                 - 0: false */
  int           irijco;       /* coupled solving of Rij
                                 - 1: true
                                 - 0: false (default) */
  int           irijnu;       /* pseudo eddy viscosity in the matrix of momentum
                                 equation to partially implicit div( rho R )
                                 - 1: true
                                 - 0: false (default) */
  int           irijrb;       /* accurate treatment of R at the boundary (see
                                 \ref condli)
                                 - 1: true
                                 - 0: false (default) */
  int           irijec;       /* wall echo term of R
                                 - 1: true
                                 - 0: false (default) */
  int           idifre;       /* whole treatment of the diagonal part of the
                                 diffusion tensor of R and epsilon
                                 - 1: true (default)
                                 - 0: simplified treatment */
  int           iclsyr;       /* partial implicitation of symmetry BCs of R
                                 - 1: true (default)
                                 - 0: false */
  int           iclptr;       /* partial implicitation of wall BCs of R
                                 - 1: true
                                 - 0: false (default) */
  double        xlomlg;       /* mixing length */

} cs_turb_rans_model_t;

/* LES turbulence model descriptor */
/*---------------------------------*/

typedef struct {

  int           idries;       /* Van Driest smoothing at the wall (only for
                                 itytur=4)
                                 - 1: true
                                 - 0: false */
  int           ivrtex;       /* vortex method (in LES)
                                 - 1: true
                                 - 0: false (default) */

} cs_turb_les_model_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to main turbulence model descriptor structure */

extern const cs_turb_model_t         *cs_glob_turb_model;

/* Pointer to reference values for turbulence descriptor structure */

extern const cs_turb_ref_values_t    *cs_glob_turb_ref_values;

/* Pointer to RANS turbulence model descriptor structure */

extern const cs_turb_rans_model_t    *cs_glob_turb_rans_model;

/* Pointer to LES turbulence model descriptor structure */

extern const cs_turb_les_model_t     *cs_glob_turb_les_model;

/* Constant for turbulence models */

extern const double cs_turb_xkappa;
extern const double cs_turb_vdriest;
extern const double cs_turb_cstlog;
extern const double cs_turb_cstlog_rough;
extern double cs_turb_cstlog_alpha;
extern const double cs_turb_apow;
extern const double cs_turb_bpow;
extern double cs_turb_dpow;
extern double cs_turb_cmu;
extern double cs_turb_cmu025;
extern const double cs_turb_ce1;
extern const double cs_turb_ce2;
extern const double cs_turb_ce4;
extern const double cs_turb_sigmak;
extern double cs_turb_sigmae;
extern double cs_turb_crij1;
extern double cs_turb_crij2;
extern double cs_turb_crij3;
extern const double cs_turb_crijp1;
extern const double cs_turb_crijp2;
extern const double cs_turb_cssge2;
extern const double cs_turb_cssgs1;
extern const double cs_turb_cssgs2;
extern const double cs_turb_cssgr1;
extern const double cs_turb_cssgr2;
extern const double cs_turb_cssgr3;
extern const double cs_turb_cssgr4;
extern const double cs_turb_cssgr5;
extern const double cs_turb_cebms1;
extern const double cs_turb_cebms2;
extern const double cs_turb_cebmr1, cebmr2, cebmr3, cebmr4, cebmr5;
extern double cs_turb_csrij;
extern const double cs_turb_cebme2;
extern const double cs_turb_cebmmu;
extern const double cs_turb_xcl;
extern const double cs_turb_xa1;
extern const double cs_turb_xct;
extern const double cs_turb_xceta;
extern const double cs_turb_cpale1;
extern const double cs_turb_cpale2;
extern const double cs_turb_cpale3;
extern const double cs_turb_cpale4;
extern const double cs_turb_cpalse;
extern const double cs_turb_cpalmu;
extern const double cs_turb_cpalc1;
extern const double cs_turb_cpalc2;
extern const double cs_turb_cpalct;
extern const double cs_turb_cpalcl;
extern const double cs_turb_cpalet;
extern const double cs_turb_ckwsk1;
extern const double cs_turb_ckwsk2;
extern const double cs_turb_ckwsw1;
extern const double cs_turb_ckwsw2;
extern const double cs_turb_ckwbt1;
extern const double cs_turb_ckwbt2;
extern double cs_turb_ckwgm1;
extern double cs_turb_ckwgm2;
extern const double cs_turb_ckwa1;
extern const double cs_turb_ckwc1;
extern const double cs_turb_csab1;
extern const double cs_turb_csab2;
extern const double cs_turb_csasig;
extern const double cs_turb_csav1;
extern double cs_turb_csaw1;
extern const double cs_turb_csaw2;
extern const double cs_turb_csaw3;
extern const double cs_turb_cssr1;
extern const double cs_turb_cssr2;
extern const double cs_turb_cssr3;
extern const double cs_turb_ccaze2;
extern const double cs_turb_ccazsc;
extern const double cs_turb_ccaza;
extern const double cs_turb_ccazb;
extern const double cs_turb_ccazc;
extern const double cs_turb_ccazd;
extern const double cs_turb_xlesfl;
extern const double cs_turb_ales;
extern const double cs_turb_bles;
extern double cs_turb_csmago;
extern const double cs_turb_xlesfd;
extern double cs_turb_smagmx;
extern double cs_turb_smagmn;
extern const double cs_turb_cdries;
extern const double cs_turb_cv2fa1;
extern const double cs_turb_cv2fe2;
extern const double cs_turb_cv2fmu;
extern const double cs_turb_cv2fc1;
extern const double cs_turb_cv2fc2;
extern const double cs_turb_cv2fct;
extern const double cs_turb_cv2fcl;
extern const double cs_turb_cv2fet;
extern double cs_turb_cwale;
extern const double cs_turb_xiafm;
extern const double cs_turb_etaafm;
extern const double cs_turb_c1trit;
extern const double cs_turb_c2trit;
extern const double cs_turb_c3trit;
extern const double cs_turb_c4trit;
extern const double cs_turb_cthafm;
extern const double cs_turb_cthdfm;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize turbulence model structures
 *----------------------------------------------------------------------------*/

void
cs_turb_model_init(void);

/*----------------------------------------------------------------------------
 * Set type and order of the turbulence model
 *----------------------------------------------------------------------------*/

void
cs_set_type_order_turbulence_model(void);

/*----------------------------------------------------------------------------
 * Set global pointer to turbulence model structure
 *----------------------------------------------------------------------------*/

void
cs_set_glob_turb_model(void);

/*----------------------------------------------------------------------------
 * Provide write access to turbulence model structure
 *----------------------------------------------------------------------------*/

cs_turb_model_t *
cs_get_glob_turb_model(void);

/*----------------------------------------------------------------------------
 * Compute turbulence model constants,
 * some of which may depend on the model choice.
 *----------------------------------------------------------------------------*/

void
cs_turb_compute_constants(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_turb_ref_values
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_ref_values_t *
cs_get_glob_turb_ref_values(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_turb_rans_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_rans_model_t *
cs_get_glob_turb_rans_model(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_turb_les_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_les_model_t *
cs_get_glob_turb_les_model(void);

/*----------------------------------------------------------------------------*
 * Print the turbulence model parameters to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_turb_model_log_setup(void);

/*----------------------------------------------------------------------------*
 * Print the turbulent constants to setup.log.
 *----------------------------------------------------------------------------*/

void
cs_turb_constants_log_setup(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_MODEL_H__ */

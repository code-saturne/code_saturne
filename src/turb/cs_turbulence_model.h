#ifndef __CS_TURBULENCE_MODEL_H__
#define __CS_TURBULENCE_MODEL_H__

/*============================================================================
 * Base turbulence model data.
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

typedef enum {

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

} cs_turb_model_type_t;

/*----------------------------------------------------------------------------
 * turbulence type of model
 *----------------------------------------------------------------------------*/

enum {

  CS_TURB_TYPE_NONE = 0,
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
 * turbulence model High or Low Reynolds number
 *----------------------------------------------------------------------------*/

enum {

  CS_TURB_HIGH_RE = 0,
  CS_TURB_LOW_RE = 1,
  CS_TURB_HIGH_LOW_RE = 2

};

/*----------------------------------------------------------------------------
 * hybrid models
 *----------------------------------------------------------------------------*/

enum {

  CS_HYBRID_NONE  = 0,
  CS_HYBRID_DES   = 1,
  CS_HYBRID_DDES  = 2,
  CS_HYBRID_SAS   = 3,
  CS_HYBRID_HTLES = 4

};

/* turbulence model general options descriptor */
/*---------------------------------------------*/

typedef struct {

  union {
    int           model;/*! turbulence model
                          - CS_TURB_NONE: no turbulence model (laminar flow)
                          - CS_TURB_MIXING_LENGTH: mixing length model
                          - CS_TURB_K_EPSILON: standard k-epsilon model
                          - CS_TURB_K_EPSILON_LIN_PROD: k-epsilon model with
                              Linear Production (LP) correction
                          - CS_TURB_K_EPSILON_LS: Launder-Sharma low Re
                              k-epsilon model
                          - CS_TURB_K_EPSILON_QUAD: Baglietto et al. low Re
                              k epsilon model
                          - CS_TURB_RIJ_EPSILON_LRR: Rij-epsilon (LRR)
                          - CS_TURB_RIJ_EPSILON_SSG: Rij-epsilon (SSG)
                          - CS_TURB_RIJ_EPSILON_EBRSM: Rij-epsilon (EBRSM)
                          - CS_TURB_LES_SMAGO_CONST: LES
                              (constant Smagorinsky model)
                          - CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic
                              Smagorisky model)
                          - CS_TURB_LES_WALE: LES (WALE)
                          - CS_TURB_V2F_PHI: v2f phi-model
                          - CS_TURB_V2F_BL_V2K: v2f BL-v2-k
                          - CS_TURB_K_OMEGA: k-omega SST
                          - CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model */
    int           iturb;  /*! Deprecated */
  };
  int           itytur;       /* class of turbulence model (integer value
                                 iturb/10) */
  int           hybrid_turb;  /*! Type of Hybrid Turbulence Model
                                   - CS_HYBRID_NONE:  No model
                                   - CS_HYBRID_DES:   Detached Eddy Simulation
                                   - CS_HYBRID_DDES:  Delayed Detached Eddy
                                                      Simulation
                                   - CS_HYBRID_SAM:   Scale Adaptive Model
                                   - CS_HYBRID_HTLES: Hybrid Temporal Large
                                                      Eddy Simulation */
  int           type;  /*! Type of turbulence modelling:
                          - CS_TURB_NONE: No model
                          - CS_TURB_RANS: RANS modelling
                          - CS_TURB_LES: LES modelling
                          - CS_TURB_HYBRID: RANS -- LES modelling */
  int           order; /*! Order of the turbulence model:
                          - CS_TURB_ALGEBRAIC: 0th order algebraik model
                          - CS_TURB_FIRST_ORDER: 1st order Eddy Viscosity
                                                 type models
                          - CS_TURB_SECOND_ORDER: 2nd order Differential
                                                  Reynolds Stress type models */
  int           high_low_re; /*! High or Low Reynolds number model:
                               - CS_TURB_HIGH_RE
                               - CS_TURB_LOW_RE
                               - CS_TURB_HIGH_LOW_RE */
} cs_turb_model_t;

/* Reference values for turbulence structure and associated pointer */
/*------------------------------------------------------------------*/

typedef struct {

  double        almax;        /* characteristic macroscopic length of the
                                 domain */
  double        uref;         /* characteristic flow velocity */

} cs_turb_ref_values_t;

/* RANS turbulence model descriptor */
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
                                 - 0: scalar diffusivity (Shir model, default)
                                 - 1: tensorial diffusivity (Daly and Harlow
                                      model) */
  int           iclkep;       /* clipping of k and epsilon
                                 - 0: absolute value clipping
                                 - 1: coupled clipping based on physical
                                      relationships */
  int           igrhok;       /* take (2/3 rho grad k) in the momentum
                                 equation
                                 - 1: true
                                 - 0: false (default) */
  int           has_buoyant_term;
                              /* take buoyant term in k-epsilon or Rij-epsilon
                               * models
                                 - 1: true (default if rho is variable)
                                 - 0: false
                                 Useful if and only if RANS models are activated
                                 and gravity is non-zero. */
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
                                 equation to partially implicit div(rho R)
                                 - 0: false (default)
                                 - 1: true
                                 - 2: Rusanov fluxes */
  int           irijrb;       /* accurate treatment of R at the boundary (see
                                 \ref cs_boundary_condition_set_coeffs)
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
  int           ikwcln;       /* Wall boundary condition on omega in k-omega SST
                                 0: Deprecated Neumann boundary condition
                                 1: Dirichlet boundary condition consistent
                                    with Menter's
                                    original model: w_wall = 60*nu/(beta*d**2) */

  double        xlomlg;       /* mixing length */

  int           dissip_buo_mdl;
                              /* Turbulent dissipation buoyant production model
                                 0: Default: Production term clipped to 0
                                 1: For EM-RSM */
} cs_turb_rans_model_t;

/* LES turbulence model descriptor */
/*---------------------------------*/

typedef struct {

  int           idries;       /* Van Driest smoothing at the wall (only for
                                 itytur=4)
                                 - 1: true
                                 - 0: false */

} cs_turb_les_model_t;

/* Hybrid turbulence model descriptor */
/*------------------------------------*/

typedef struct {

  int           iicc;      /* Internal Consistency Constraint applied
                              (only for hybrid_turb=4)
                                 - 1: true
                                 - 0: false */
  int           ishield;   /* Shielding function applied at the wall
                              (only for hybrid_turb=4)
                                 - 1: true
                                 - 0: false */

  cs_lnum_t     n_iter_mean; /* number of iteration for the exponential mean */
  cs_lnum_t     time_mean;   /* time for the exponential mean, automatically
                                computed if n_iter_mean > 0 */

} cs_turb_hybrid_model_t;

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

/* Pointer to hybrid turbulence model descriptor structure */

extern const cs_turb_hybrid_model_t  *cs_glob_turb_hybrid_model;

/* Constant for turbulence models */

extern double cs_turb_xkappa;
extern double cs_turb_vdriest;
extern double cs_turb_cstlog;
extern double cs_turb_cstlog_rough;
extern double cs_turb_cstlog_alpha;
extern double cs_turb_apow;
extern double cs_turb_bpow;
extern double cs_turb_dpow;
extern double cs_turb_cmu;
extern double cs_turb_cmu025;
extern double cs_turb_ce1;
extern double cs_turb_ce2;
extern double cs_turb_ce3;
extern double cs_turb_ce4;
extern double cs_turb_crij_eps;
extern double cs_turb_crij1;
extern double cs_turb_crij2;
extern double cs_turb_crij3;
extern double cs_turb_crij_c0;
extern double cs_turb_crijp1;
extern double cs_turb_crijp2;
extern double cs_turb_cssgs1;
extern double cs_turb_cssgs2;
extern double cs_turb_cssgr1;
extern double cs_turb_cssgr2;
extern double cs_turb_cssgr3;
extern double cs_turb_cssgr4;
extern double cs_turb_cssgr5;
extern double cs_turb_cebms1;
extern double cs_turb_cebms2;
extern double cs_turb_cebmr1;
extern double cs_turb_cebmr2;
extern double cs_turb_cebmr3;
extern double cs_turb_cebmr4;
extern double cs_turb_cebmr5;
extern double cs_turb_csrij;
extern double cs_turb_cebmmu;
extern double cs_turb_xcl;
extern double cs_turb_xa1;
extern double cs_turb_xct;
extern double cs_turb_xclt;
extern double cs_turb_xceta;
extern double cs_turb_cpale1;
extern double cs_turb_cpale2;
extern double cs_turb_cpale3;
extern double cs_turb_cpale4;
extern double cs_turb_cpalc1;
extern double cs_turb_cpalc2;
extern double cs_turb_cpalct;
extern double cs_turb_cpalcl;
extern double cs_turb_cpalet;
extern double cs_turb_ckwsk1;
extern double cs_turb_ckwsk2;
extern double cs_turb_ckwsw1;
extern double cs_turb_ckwsw2;
extern double cs_turb_ckwbt1;
extern double cs_turb_ckwbt2;
extern double cs_turb_ckwgm1;
extern double cs_turb_ckwgm2;
extern double cs_turb_ckwa1;
extern double cs_turb_ckwc1;
extern double cs_turb_cddes;
extern double cs_turb_csas;
extern double cs_turb_csas_eta2;
extern double cs_turb_chtles_bt0;
extern double cs_turb_cnl1;
extern double cs_turb_cnl2;
extern double cs_turb_cnl3;
extern double cs_turb_cnl4;
extern double cs_turb_cnl5;
extern double cs_turb_csab1;
extern double cs_turb_csab2;
extern double cs_turb_csasig;
extern double cs_turb_csav1;
extern double cs_turb_csaw1;
extern double cs_turb_csaw2;
extern double cs_turb_csaw3;
extern double cs_turb_cssr1;
extern double cs_turb_cssr2;
extern double cs_turb_cssr3;
extern double cs_turb_ccaze2;
extern double cs_turb_ccazsc;
extern double cs_turb_ccaza;
extern double cs_turb_ccazb;
extern double cs_turb_ccazc;
extern double cs_turb_ccazd;
extern double cs_turb_xlesfl;
extern double cs_turb_ales;
extern double cs_turb_bles;
extern double cs_turb_csmago;
extern double cs_turb_xlesfd;
extern double cs_turb_csmago_max;
extern double cs_turb_csmago_min;
extern double cs_turb_cdries;
extern double cs_turb_cv2fa1;
extern double cs_turb_cv2fe2;
extern double cs_turb_cv2fc1;
extern double cs_turb_cv2fc2;
extern double cs_turb_cv2fct;
extern double cs_turb_cv2fcl;
extern double cs_turb_cv2fet;
extern double cs_turb_cwale;
extern double cs_turb_xiafm;
extern double cs_turb_etaafm;
extern double cs_turb_c1trit;
extern double cs_turb_c2trit;
extern double cs_turb_c3trit;
extern double cs_turb_c4trit;
extern double cs_turb_cthafm;
extern double cs_turb_cthdfm;
extern double cs_turb_cthebdfm;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize additional turbulence model members of turbulence model
 * and RANS model structure
 *----------------------------------------------------------------------------*/

void
cs_turbulence_init_models(void);

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
 *
 * \param[in]       phase_id  turbulent phase id (-1 for single phase flow)
 *
 *----------------------------------------------------------------------------*/

void
cs_turb_compute_constants(int phase_id);

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

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_turb_hybrid_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_hybrid_model_t *
cs_get_glob_turb_hybrid_model(void);

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
/*!
 * \brief Compute characteristic length for turbulence if not already done.
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_init_ref_quantities(void);

/*----------------------------------------------------------------------------*
 * Clip turbulent fluxes
 *----------------------------------------------------------------------------*/

void
cs_clip_turbulent_fluxes(int  flux_id,
                         int  ivartt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return or estimate the value of the turbulent kinetic energy
 *        over specified elements.
 *
 * Returned values are zero for turbulence models other than RANS.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_function_k(int              location_id,
                         cs_lnum_t         n_elts,
                         const cs_lnum_t  *elt_ids,
                         void             *input,
                         void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return or estimate the value of the turbulent dissipation
 *        over specified elements.
 *
 * Returned values are zero for turbulence models other than RANS.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_function_eps(int              location_id,
                           cs_lnum_t         n_elts,
                           const cs_lnum_t  *elt_ids,
                           void             *input,
                           void             *vals);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return or estimate the value of the Reynolds stresses
 *        over specified elements.
 *
 * Returned values are zero for turbulence models other than RANS.
 *
 * This function matches the cs_eval_at_location_t function profile.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or NULL if no
 *                               filtering is required
 * \param[in, out]  input        ignored
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_function_rij(int               location_id,
                           cs_lnum_t         n_elts,
                           const cs_lnum_t  *elt_ids,
                           void             *input,
                           void             *vals);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_MODEL_H__ */

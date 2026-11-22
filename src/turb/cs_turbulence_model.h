#ifndef CS_TURBULENCE_MODEL_H
#define CS_TURBULENCE_MODEL_H

/*============================================================================
 * Base turbulence model data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

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
  CS_TURB_RIJ_EPSILON_BFH= 33,
  CS_TURB_RIJ_OMEGA = 34,
  CS_TURB_LES_SMAGO_CONST = 40,
  CS_TURB_LES_SMAGO_DYN = 41,
  CS_TURB_LES_WALE = 42,
  CS_TURB_LES_KSGS = 43,
  CS_TURB_LES_TAUSGS = 44,
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
 * Schemes for Rij discretization
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_RIJ_SCHEME_DEFAULT = 0,
  CS_RIJ_SCHEME_IMPLICIT_VISCOSITY = 1,
  CS_RIJ_SCHEME_RUSANOV = 2,
  CS_RIJ_SCHEME_GODUNOV = 3

} cs_rij_scheme_t;

/*----------------------------------------------------------------------------
 * Schemes for Rij/Rit source time stepping
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_TURB_RIJ_SOURCE_TS_CONVEXP = -1,  /* Convection full explicit */
  CS_TURB_RIJ_SOURCE_TS_IMEX = 0,      /* Standard implicit/explicit */
  CS_TURB_RIJ_SOURCE_TS_EXPONENTIAL = 1, /* Exact exponential integration */
  CS_TURB_RIJ_SOURCE_TS_VAR_TAU = 2    /* Variable tau(t) trajectory */

} cs_turb_rij_source_ts_t;

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

/*! turbulence model general options descriptor */
/*---------------------------------------------*/

struct cs_turb_model_t {

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Default constructor with suppression of deprecated declarations,
   *        needed to hide warnings simply induced by including this header
   *        file.
   */
  /*--------------------------------------------------------------------------*/

  /*! \cond DOXYGEN_SHOULD_SKIP_THIS */
  DISABLE_WARNING_PUSH
  DISABLE_WARNING(-Wdeprecated-declarations)
  /*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

  cs_turb_model_t(){}

  /*! \cond DOXYGEN_SHOULD_SKIP_THIS */
  DISABLE_WARNING_POP
  /*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

  int  model{-999};/*!< turbulence model
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
                     - CS_TURB_RIJ_OMEGA: Rij-omega (Wilcox 2006)
                     - CS_TURB_LES_SMAGO_CONST: LES
                         (constant Smagorinsky model)
                     - CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic
                         Smagorisky model)
                     - CS_TURB_LES_WALE: LES (WALE)
                     - CS_TURB_V2F_PHI: v2f phi-model
                     - CS_TURB_V2F_BL_V2K: v2f BL-v2-k
                     - CS_TURB_K_OMEGA: k-omega SST
                     - CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model */

  [[deprecated("Use model instead")]]
  int& iturb{model};    /*!< Deprecated, use `model` instead */

  int  itytur{-999};    /*!< class of turbulence model (integer value
                        model/10) */
  int  hybrid_turb{0};  /*!< Type of Hybrid Turbulence Model
                            - CS_HYBRID_NONE:  No model
                            - CS_HYBRID_DES:   Detached Eddy Simulation
                            - CS_HYBRID_DDES:  Delayed Detached Eddy
                                               Simulation
                            - CS_HYBRID_SAS:   Scale Adaptive Model
                            - CS_HYBRID_HTLES: Hybrid Temporal Large
                                               Eddy Simulation */
  int  type{-1};  /*!< Type of turbulence modelling:
                     - CS_TURB_NONE: No model
                     - CS_TURB_RANS: RANS modelling
                     - CS_TURB_LES: LES modelling
                     - CS_TURB_HYBRID: RANS -- LES modelling */
  int  order{-1}; /*!< Order of the turbulence model:
                     - CS_TURB_ALGEBRAIC: 0th order algebraik model
                     - CS_TURB_FIRST_ORDER: 1st order Eddy Viscosity
                                            type models
                     - CS_TURB_SECOND_ORDER: 2nd order Differential
                                             Reynolds Stress type models */
  int  high_low_re{-1}; /*!< High or Low Reynolds number model:
                          - CS_TURB_HIGH_RE
                          - CS_TURB_LOW_RE
                          - CS_TURB_HIGH_LOW_RE */
};

/*! Reference values for turbulence structure and associated pointer */
/*------------------------------------------------------------------*/

typedef struct {

  double        almax; /*!< characteristic macroscopic length of the
                          domain */
  double        uref;  /*!< characteristic flow velocity */

} cs_turb_ref_values_t;

/*! RANS turbulence model descriptor */
/*----------------------------------*/

typedef struct cs_turb_rans_model_t {

  int     irccor;       /*!< activation of rotation/curvature correction for
                           an eddy viscosity turbulence models
                           - 0: false
                           - 1: true */
  int     itycor;       /*!< type of rotation/curvature correction for an
                           eddy viscosity turbulence models
                           - 1: Cazalbou correction (default when irccor=1
                                and itytur=2 or 5)
                           - 2: Spalart-Shur correction (default when
                                irccor=1 and iturb=60 or 70) */
  int     idirsm;       /*!< turbulent diffusion model for second moment
                           closure
                           - 0: scalar diffusivity (Shir model, default)
                           - 1: tensorial diffusivity (Daly and Harlow
                                model) */
  int     igrhok;       /*!< take (2/3 rho grad k) in the momentum
                           equation
                           - 1: true
                           - 0: false (default) */

  int     has_buoyant_term; /*!< take buoyant term in k-epsilon or Rij-epsilon
                             * models
                               - 1: true (default if rho is variable)
                               - 0: false
                               Useful if and only if RANS models are activated
                               and gravity is non-zero. */
  int     ikecou;       /*!< partially coupled version of
                           k-epsilon (only for iturb=20)
                           - 1: true (default)
                           - 0: false */
  int     reinit_turb;  /*!< Advanced re-init for EBRSM and k-omega models
                           - 1: true (default)
                           - 0: false */
  int     irijco;       /*!< coupled solving of Rij
                           - 1: true
                           - 0: false (default) */
  int     rij_discretization_scheme{0}; /*!< Rij discretization scheme:
                                           - 0: CS_RIJ_SCHEME_DEFAULT
                                                (false, default)
                                           - 1: CS_RIJ_SCHEME_IMPLICIT_VISCOSITY
                                                (true)
                                           - 2: CS_RIJ_SCHEME_RUSANOV (Rusanov)
                                           - 3: CS_RIJ_SCHEME_GODUNOV
                                                (exact Godunov) */
  [[deprecated("renamed to rij_discretization_scheme")]] \
  int&    irijnu{rij_discretization_scheme}; /*!< Deprecated */
  int     verbosity;    /*!< verbosity level for Rij / Rit schemes
                           - 0: silent
                           - 1: moderate (default)
                           - >= 2: active logs on all cells */
  int     source_time_stepping{0};
                         /*!< treatment of the Rotta-Monin source step
                           (return-to-isotropy, thermal relaxation,
                           Boussinesq buoyancy) for the Reynolds-stress,
                           turbulent-heat-flux, temperature-variance,
                           and epsilon equations, once mean-gradient
                           production and diffusion have been handled
                           elsewhere. Requires:
                           rij_discret_sch == CS_RIJ_SCHEME_GODUNOV
                           and, for now, pure Rotta closure
                           (crij2 == 0).
                           - CS_TURB_RIJ_SOURCE_TS_CONVEXP (-1)
                           - CS_TURB_RIJ_SOURCE_TS_IMEX (0, default)
                           - CS_TURB_RIJ_SOURCE_TS_EXPONENTIAL (1)
                           - CS_TURB_RIJ_SOURCE_TS_VAR_TAU (2) */
  int     irijrb;       /*!< accurate treatment of R at the boundary (see
                           \ref cs_boundary_conditions_set_coeffs)
                           - 1: true
                           - 0: false (default) */
  int     irijec;       /*!< wall echo term of R
                           - 1: true
                           - 0: false (default) */
  int     iclsyr;       /*!< partial implicitation of symmetry BCs of R
                           - 1: true (default)
                           - 0: false */
  int     iclptr;       /*!< partial implicitation of wall BCs of R
                           - 1: true
                           - 0: false (default) */
  int     ikwcln;       /*!< Wall boundary condition on omega in k-omega SST
                           0: Deprecated Neumann boundary condition
                           1: Dirichlet boundary condition consistent
                              with Menter's
                              original model: w_wall = 60*nu/(beta*d**2) */

  double  xlomlg;       /*!< mixing length */

  int     dissip_buo_mdl; /*!< Turbulent dissipation buoyant production model
                               0: Default: Production term clipped to 0
                               1: For EM-RSM */

  /* BFH second-moment closure options */
  int     bfh_slow_model;      /*!< Slow redistribution model for BFH:
                                    0: Rotta (default)
                                    1: Chung & Kim
                                    2: Speziale & Sarkar
                                    3: Universal 2C-3C */
  double  bfh_c2_prime_inf;    /*!< Asymptotic constant for Chung & Kim
                                    (default: 15.5) */
  double  bfh_aw;              /*!< Wall damping coefficient Aw
                                    (default: 150.0) */
  double  bfh_pw;              /*!< Wall damping power exponent pw
                                    (default: 2.0) */
  int     bfh_use_beta2;       /*!< Use beta2 damping on slow term:
                                    0: false (default)
                                    1: true */
  int     bfh_ck_damping;      /*!< Reynolds-number function of the
                                    Chung & Kim slow term:
                                    1: 1/beta(R_t) of the original model
                                       (default), calibrated on homogeneous
                                       decaying turbulence
                                    0: the tau/T2 activation of the power
                                       mean, shared with the Rotta branch */
  int     bfh_beta2_form;      /*!< Form of the rapid coefficient beta2:
                                    0: 9/10 lambda_min / tr(R) (default,
                                       90 percent of the stability bound)
                                    1: 3/5 lambda_min
                                       / (lambda_max + lambda_min) */
  int     bfh_use_durbin_time; /*!< Use Durbin time scale T_drb for Rotta:
                                    0: false (default: T = k/eps or dt_wall)
                                    1: true: T_R = T_drb, T_eps = k/eps */
  double  bfh_t2_ct;           /*!< Constant C_T2 of the Rotta redistribution
                                    time scale T2 = (k/eps) [1 + (Re_s/Re_t)^p]
                                    with Re_s = C_T2^2, used when
                                    bfh_use_durbin_time = 1
                                    (default: -1, i.e. cs_turb_xct) */
  double  bfh_t2_pow;          /*!< Exponent p of the Rotta redistribution
                                    time scale (default: 0.5). It fixes the
                                    near-wall asymptotics: with p >= 1 the
                                    isotropic Rotta source vanishes at least
                                    as y^4 at the wall */
  double  bfh_t2_n;            /*!< Order n of the power mean blending the
                                    integral and near-wall time scales,
                                    T2 = tau (1 + X^n)^(1/n) with
                                    X = (Re_s/Re_t)^p. It fixes the sharpness
                                    of the transition only, not the
                                    asymptotics: n = 1 is the harmonic blend
                                    of the two rates, n -> infinity gives
                                    max(tau, T_v) (default: 1) */
  int     bfh_c2_mod;          /*!< Reynolds-number function carried by the
                                    quadratic constant C'_2 of the Chung &
                                    Kim slow term, independently of the
                                    function carried by the overall
                                    amplitude (bfh_ck_damping):
                                    - 0: 1 / beta(R_t) (Chung & Kim 1995)
                                    - 1: tau / T2 sigmoid of the power mean
                                    - 2: no damping, C'_2 constant
                                    Mode 2 makes a negative C'_2 meaningful,
                                    which reverses the transfer driven by the
                                    quadratic block while the Schumann
                                    two-component constraint A_1 = 2 +
                                    A_2 (II_b - 1/3) still holds, so
                                    realizability is preserved for either
                                    sign (default: 0) */
  double  bfh_ck_cr_add;       /*!< Additional linear return to isotropy
                                    added to the Chung & Kim slow term,
                                    - dC_R F (eps / k) dev(R), where
                                    F = 27 det R / (tr R)^3 is the flatness
                                    parameter. F vanishes on the
                                    two-component surface, so the extra
                                    return does not break realizability,
                                    while it lifts the linear coefficient
                                    above the value A_1 = 2 at which the
                                    Schumann constraint locks it in the
                                    fully three-component state
                                    (default: 0) */
  double  bfh_wall_split;      /*!< Unequal tangential split of the wall
                                    tensor. The redistribution
                                    - (5/2) R_nn P of the trace-free wall
                                    model shares equally between the two
                                    tangential directions; it is replaced by
                                    - R_nn [(5/2 + d) t.t + (5/2 - d) s.s],
                                    where t is the unit vector along the
                                    tangential part of R.n (the shear
                                    direction) and s = n x t. The sum of the
                                    two coefficients is unchanged, so W stays
                                    trace free, the wall asymptotics are
                                    untouched, and d > 0 removes more from
                                    the streamwise component and relatively
                                    less from the spanwise one
                                    (default: 0, equal split) */
  double  bfh_wall_split_eta;  /*!< Regularisation of the shear-aligned frame
                                    of the wall tensor split: the split is
                                    weighted by
                                    |r_t|^2 / (|r_t|^2 + (eta k)^2) so that
                                    it switches off smoothly where the
                                    tangential part r_t of R.n no longer
                                    defines a direction. Since |r_t| / k is
                                    of order 10^-2 throughout the viscous
                                    sublayer, eta must stay well below that
                                    for the split to act where the budgets
                                    ask for it (default: 1e-3) */
  double  bfh_wall_stab;       /*!< Constant c_w of the implicit-only
                                    stabilisation - c_w R / tau of the wall
                                    tensor. The remainder left explicit by
                                    the exact implicitation of 2 imp_phiw N
                                    is not sign definite, and c_w covers it.
                                    The sufficient condition
                                    c_w / tau >= imp_phiw (n.R.n)(n.R^-1.n)
                                    is not a bound valid for every
                                    realizable state, so c_w may have to be
                                    raised on meshes whose first cell is
                                    much thinner than y+ = 0.1, where the
                                    2 nu / z^2 prefactor grows. A negative
                                    value selects the default, 6 for
                                    uncoupled components and 2 for coupled
                                    ones (default: -1) */

  /*! \cond DOXYGEN_SHOULD_SKIP_THIS */
  DISABLE_WARNING_PUSH
  DISABLE_WARNING(-Wdeprecated-declarations)
  /*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

  /* Constructor to bind reference and set defaults */
  cs_turb_rans_model_t()
    : irccor(0),
      itycor(-999),
      idirsm(0),
      igrhok(0),
      has_buoyant_term(1),
      ikecou(0),
      reinit_turb(1),
      irijco(1),
      rij_discretization_scheme(0),
      irijnu(rij_discretization_scheme),
      verbosity(0),
      source_time_stepping(0),
      irijrb(0),
      irijec(0),
      iclsyr(1),
      iclptr(0),
      ikwcln(1),
      xlomlg(-1e13),
      dissip_buo_mdl(0),
      bfh_slow_model(0),
      bfh_c2_prime_inf(15.5),
      bfh_aw(150.0),
      bfh_pw(2.0),
      bfh_use_beta2(0),
      bfh_ck_damping(1),
      bfh_beta2_form(0),
      bfh_use_durbin_time(0),
      bfh_t2_ct(-1.0),
      bfh_t2_pow(0.5),
      bfh_t2_n(1.0),
      bfh_c2_mod(0),
      bfh_ck_cr_add(0.0),
      bfh_wall_split(0.0),
      bfh_wall_split_eta(1.e-3),
      bfh_wall_stab(-1.0)
  {}
  /*! \cond DOXYGEN_SHOULD_SKIP_THIS */
  DISABLE_WARNING_POP
  /*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
extern double cs_turb_ce1;
extern double cs_turb_ce2;
extern double cs_turb_ce3;
extern double cs_turb_ce4;
extern double cs_turb_crij_eps;
extern double cs_turb_crij1;
extern double cs_turb_crij2;
extern double cs_turb_crij3;
extern double cs_turb_crij_c0;
extern double cs_turb_crij_ct;
extern double cs_turb_crijp1;
extern double cs_turb_crijp2;
extern double cs_turb_rij_omega_alpha;
extern double cs_turb_rij_omega_beta0;
extern double cs_turb_rij_omega_sigma_d0;
extern double cs_turb_rij_omega_alpha_b_star;
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

/*----------------------------------------------------------------------------*/
/*
 * \brief Return name string associated with a turbulence model value
 *
 * \param[in]  id  model type
 *
 * \return pointer to turbulence models description name.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_turbulence_model_name(cs_turb_model_type_t  id);

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

extern "C" void
cs_clip_turbulent_fluxes(int  flux_id,
                         int  ivartt);

/*----------------------------------------------------------------------------*/
/*
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

extern "C" void
cs_turbulence_function_k(int              location_id,
                         cs_lnum_t         n_elts,
                         const cs_lnum_t  *elt_ids,
                         void             *input,
                         void             *vals);

/*----------------------------------------------------------------------------*/
/*
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

extern "C" void
cs_turbulence_function_eps(int              location_id,
                           cs_lnum_t         n_elts,
                           const cs_lnum_t  *elt_ids,
                           void             *input,
                           void             *vals);

/*----------------------------------------------------------------------------*/
/*
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

extern "C" void
cs_turbulence_function_rij(int               location_id,
                           cs_lnum_t         n_elts,
                           const cs_lnum_t  *elt_ids,
                           void             *input,
                           void             *vals);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENCE_MODEL_H */

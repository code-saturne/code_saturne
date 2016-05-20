/*============================================================================
 * Base turbulence model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_model.c
        Base turbulence model data.
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_model_t

  \brief Turbulence model general options descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_model_t::iturb
        turbulence model
        - 0: no turbulence model (laminar flow)
        - 10: mixing length model
        - 20: standard \f$ k-\varepsilon \f$ model
        - 21: \f$ k-\varepsilon \f$ model with Linear Production (LP) correction
        - 30: \f$ R_{ij}-\epsilon \f$ (LRR)
        - 31: \f$ R_{ij}-\epsilon \f$ (SSG)
        - 32: \f$ R_{ij}-\epsilon \f$ (EBRSM)
        - 40: LES (constant Smagorinsky model)
        - 41: LES ("classical" dynamic Smagorisky model)
        - 42: LES (WALE)
        - 50: v2f phi-model
        - 51: v2f \f$ BL-v^2-k \f$
        - 60: \f$ k-\omega \f$ SST
        - 70: Spalart-Allmaras model
  \var  cs_turb_model_t::itytur
        class of turbulence model (integer value iturb/10)
  \var  cs_turb_model_t::nvarcl
        number of variable plus number of turbulent fluxes (used by the
        boundary conditions)
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_rans_model_t

  \brief RANS turbulence model descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_rans_model_t::irccor
        activation of rotation/curvature correction for an eddy viscosity
        turbulence models
        - 0: false
        - 1: true
  \var  cs_turb_rans_model_t::itycor
        type of rotation/curvature correction for an eddy viscosity turbulence
        models
        - 1: Cazalbou correction (default when irccor=1 and itytur=2 or 5)
        - 2: Spalart-Shur correction (default when irccor=1 and iturb=60 or 70)
  \var  cs_turb_rans_model_t::idirsm
        turbulent diffusion model for second moment closure
        - 0: scalar diffusivity (Shir model)
        - 1: tensorial diffusivity (Daly and Harlow model, default model)
  \var  cs_turb_rans_model_t::iclkep
        clipping of k and epsilon
        - 0: absolute value clipping
        - 1: coupled clipping based on physical relationships
  \var  cs_turb_rans_model_t::igrhok
        take \f$ 2/3 \rho \grad k \f$ in the momentum equation
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::igrake
        buoyant term in \f$ k- \varepsilon \f$
        - 1: true (default if \f$ \rho \f$ is variable)
        - 0: false
  \var  cs_turb_rans_model_t::igrari
        buoyant term in \f$ R_{ij}- \varepsilon \f$
        - 1: true (default if \f$ \rho \f$ is variable)
        - 0: false
  \var  cs_turb_rans_model_t::ikecou
        partially coupled version of \f$ k-\varepsilon \f$ (only for iturb=20)
        - 1: true (default)
        - 0: false
  \var  cs_turb_rans_model_t::reinit_turb
        Advanced re-init for EBRSM and k-omega models
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijco
        coupled solving of Rij
        - 1: true
        - 0: false (default)
        \var  cs_turb_rans_model_t::irijnu
        pseudo eddy viscosity in the matrix of momentum equation to partially
        implicit \f$ \divv \left( \rho \tens{R} \right) \f$
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijrb
        accurate treatment of \f$ \tens{R} \f$ at the boundary (see \ref condli)
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijec
        wall echo term of \f$ \tens{R} \f$
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::idifre
        whole treatment of the diagonal part of the diffusion tensor of \f$
        \tens{R} \f$ and \f$ \varepsilon \f$
        - 1: true (default)
        - 0: simplified treatment
  \var  cs_turb_rans_model_t::iclsyr
        partial implicitation of symmetry BCs of \f$ \tens{R} \f$
        - 1: true (default)
        - 0: false
  \var  cs_turb_rans_model_t::iclptr
        partial implicitation of wall BCs of \f$ \tens{R} \f$
        - 1: true
        - 0: false (default)
  \var  cs_turb_ref_values_t::almax
        characteristic macroscopic length of the domain, used for the
        initialization of the turbulence and the potential clipping (with
        \ref iclkep=1)
        - Negative value: not initialized (the code then uses the cubic root of
          the domain volume).

        Useful mainly for RANS models.
  \var  cs_turb_ref_values_t::uref
        characteristic flow velocity, used for the initialization of the
        turbulence
        - Negative value: not initialized.

        Useful mainly for RANS models and if
        the turbulence is not initialized somewhere else (restart file or
        subroutine \ref cs\_user\_initialization).
  \var  cs_turb_rans_model_t::xlomlg
        mixing length for the mixing length model

        Useful if and only if \ref iturb= 10 (mixing length).
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_les_model_t

  \brief LES turbulence model descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_les_model_t::idries
        Van Driest smoothing at the wall (only for itytur=4)
        - 1: true
        - 0: false
  \var  cs_turb_les_model_t::ivrtex
        vortex method (in LES)
        - 1: true
        - 0: false (default)
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* main turbulence model structure and associated pointer */

static cs_turb_model_t  _turb_model =
{
  .iturb  = -999,
  .itytur = -999,
  .nvarcl = 0
};

const cs_turb_model_t  *cs_glob_turb_model = &_turb_model;

/* Reference values for turbulence structure and associated pointer */

static cs_turb_ref_values_t
_turb_ref_values =
{
  .almax      = -999,
  .uref       =-1e13
};

const cs_turb_ref_values_t  *cs_glob_turb_ref_values = &_turb_ref_values;

/* RANS turbulence model structure and associated pointer */

static cs_turb_rans_model_t
_turb_rans_model =
{
  .irccor     =    0,
  .itycor     = -999,
  .idirsm     =    1,
  .iclkep     =    0,
  .igrhok     =    0,
  .igrake     =    1,
  .igrari     =    1,
  .ikecou     = -999,
  .reinit_turb=    0,
  .irijco     =    0,
  .irijnu     =    0,
  .irijrb     =    0,
  .irijec     =    0,
  .idifre     =    1,
  .iclsyr     =    1,
  .iclptr     =    0,
  .xlomlg      =-1e13
};

const cs_turb_rans_model_t  *cs_glob_turb_rans_model = &_turb_rans_model;

/* LES turbulence model structure and associated pointer */

static cs_turb_les_model_t  _turb_les_model = {-1, 0};

const cs_turb_les_model_t  *cs_glob_turb_les_model = &_turb_les_model;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*!
 * Karman constant. (= 0.42)
 *
 * Useful if and only if \ref iturb >= 10.
 *  (mixing length, \f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$,
 * LES, v2f or \f$k-\omega\f$).
 */
const double cs_turb_xkappa = 0.42;

/*!
 * Van Driest constant. (= 25.6)
 *
 * Useful if and only if \ref cs_wall_functions_t::iwallf
 * "cs_glob_wall_functions::iwallf" = 5.
 *  (Two scales log law at the wall using Van Driest mixing length expression).
 */
const double cs_turb_vdriest = 25.6;

/*!
 * Constant of logarithmic smooth law function:
 * \f$ \dfrac{1}{\kappa} \ln(y^+) + cstlog \f$
 * (\f$ cstlog = 5.2 \f$).
 *
 * Constant of the logarithmic wall function.
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
const double cs_turb_cstlog = 5.2;

/*!
 * Constant of logarithmic rough law function:
 * \f$ \dfrac{1}{\kappa} \ln(y/\xi) + cstlog_{rough} \f$
 * (\f$ cstlog_{rough} = 8.5 \f$).
 *
 * Constant of the logarithmic wall function.
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
const double cs_turb_cstlog_rough = 8.5;

/*!
 * Constant \f$ \alpha \f$ for logarithmic law function switching from rough to smooth:
 * \f$ \dfrac{1}{\kappa} \ln(y u_k/(\nu + \alpha \xi u_k)) + cstlog \f$
 * (\f$ \alpha = \exp \left( -\kappa (8.5 - 5.2) \right) \f$).
 *
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
double cs_turb_cstlog_alpha;

/*! Werner and Wengle coefficient */
const double cs_turb_apow = 8.3;

/*! Werner and Wengle coefficient */
const double cs_turb_bpow = 1.0/7.0;

/*! Werner and Wengle coefficient */
double cs_turb_dpow;

/*!
 * Constant \f$C_\mu\f$ for all the RANS turbulence models except for the
 * v2f model (see \ref cs_turb_cv2fmu for the value of \f$C_\mu\f$ in case of v2f
 * modelling). Useful if and only if \ref iturb = 20, 21, 30, 31 or 60
 * (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or \f$k-\omega\f$).
 */
const double cs_turb_cmu = 0.09;

/*! \f$ C_\mu^\frac{1}{4} \f$ */
double cs_turb_cmu025;

/*!
 * Constant \f$C_{\varepsilon 1}\f$ for all the RANS turbulence models except
 * for the v2f and the \f$k-\omega\f$ models.
 * Useful if and only if \ref iturb= 20, 21, 30 or 31 (\f$k-\varepsilon\f$
 * or \f$R_{ij}-\varepsilon\f$).
 */
const double cs_turb_ce1 = 1.44;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the \f$k-\varepsilon\f$ and
 * \f$R_{ij}-\varepsilon\f$ LRR models.
 * Useful if and only if \ref iturb = 20, 21 or 30
 * (\f$k-\varepsilon\f$ or \f$R_{ij}-\varepsilon\f$ LRR).
 */
const double cs_turb_ce2 = 1.92;

/*!
 * Coefficient of interfacial coefficient in k-eps, used in Lagrange treatment.
 *
 * Constant \f$C_{\varepsilon 4}\f$ for the interfacial term (Lagrangian module)
 * in case of two-way coupling. Useful in case of Lagrangian modelling,
 * in \f$k-\varepsilon\f$ and \f$R_{ij}-\varepsilon\f$ with two-way coupling.
 */
const double cs_turb_ce4 = 1.20;

/*!
 * Prandtl number for \f$k\f$ with \f$k-\varepsilon\f$ and v2f models.
 * Useful if and only if \ref iturb=20, 21 or 50 (\f$k-\varepsilon\f$ or v2f).
 */
const double cs_turb_sigmak = 1.0;

/*!
 * Prandtl number for \f$\varepsilon\f$.
 * Useful if and only if \ref iturb= 20, 21, 30, 31 or 50
 * (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or v2f).
 */
double cs_turb_sigmae;

/*!
 * Constant \f$C_1\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double cs_turb_crij1 = 1.80;

/*
 * Constant \f$C_2\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double cs_turb_crij2 = 0.60;

/*!
 * Constant \f$C_3\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double cs_turb_crij3 = 0.55;

/*!
 * Constant \f$C_1^\prime\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model,
 * corresponding to the wall echo terms.
 * Useful if and only if \ref iturb=30 and \ref cs_turb_rans_model_t::irijec
 * "cs_turb_rans_model_t::irijec"=1
 * (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double cs_turb_crijp1 = 0.50;

/*!
 * Constant \f$C_2^\prime\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model,
 * corresponding to the wall echo terms.
 * Useful if and only if \ref iturb=30 and \ref cs_turb_rans_model_t::irijec=1
 * (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double cs_turb_crijp2 = 0.30;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssge2 = 1.83;

/*!
 * Constant \f$C_{s1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgs1 = 1.70;

/*!
 * Constant \f$C_{s2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgs2 = -1.05;

/*!
 * Constant \f$C_{r1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgr1 = 0.90;

/*!
 * Constant \f$C_{r2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgr2 = 0.80;

/*!
 * Constant \f$C_{r3}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgr3 = 0.65;

/*!
 * constant \f$C_{r4}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgr4 = 0.625;

/*!
 * Constant \f$C_{r1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cs_turb_cssgr5 = 0.20;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_cebms1 = 1.70;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_cebms2 = 0.;

const double cs_turb_cebmr1 = 0.90;
const double cs_turb_cebmr2 = 0.80;
const double cs_turb_cebmr3 = 0.65;
const double cs_turb_cebmr4 = 0.625;
const double cs_turb_cebmr5 = 0.20;
const double cs_turb_cebmr6 = 0.6;

/*!
 * Constant \f$C_s\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_csrij;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_cebme2 = 1.83;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_cebmmu = 0.22;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_xcl = 0.122;

/*! Constant in the expression of Ce1' for the Rij-epsilon EBRSM. */
const double cs_turb_xa1 = 0.1;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_xct = 6.0;

/*! Constant of the Rij-epsilon EBRSM. */
const double cs_turb_xceta = 80.0;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpale1 = 1.44;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpale2 = 1.83;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpale3 = 2.3;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpale4 = 0.4;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalse = 1.5;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalmu = 0.22;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalc1 = 1.7;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalc2 = 0.9;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalct = 4.0;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalcl = 0.164;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cs_turb_cpalet = 75.0;

/*!
 * Constant \f$\sigma_{k1}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60.
 */
const double cs_turb_ckwsk1 = 1.0/0.85;

/*!
 * Constant \f$\sigma_{k2}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60.
 */
const double cs_turb_ckwsk2 = 1.0;

/*!
 * Constant \f$\sigma_{\omega 1}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double cs_turb_ckwsw1 = 2.0;

/*!
 * Constant \f$\sigma_{\omega 2}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double cs_turb_ckwsw2 = 1.0/0.856;

/*!
 * Constant \f$\beta_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double cs_turb_ckwbt1 = 0.075;

/*!
 * Constant \f$\beta_2\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double cs_turb_ckwbt2 = 0.0828;


/*!
 * \f$\frac{\beta_1}{C_\mu}-\frac{\kappa^2}{\sqrt{C_\mu}\sigma_{\omega 1}}\f$.
 * Constant \f$\gamma_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * \warning: \f$\gamma_1\f$ is calculated before the call to
 * \ref usipsu. Hence, if \f$\beta_1\f$, \f$C_\mu\f$, \f$\kappa\f$ or
 * \f$\sigma_{\omega 1}\f$ is modified in \ref usipsu,
 * \ref cs_turb_ckwgm1 must also be modified in accordance.
 */
double cs_turb_ckwgm1;

/*!
 * \f$\frac{\beta_2}{C_\mu}-\frac{\kappa^2}{\sqrt{C_\mu}\sigma_{\omega 2}}\f$.
 * Constant \f$\gamma_2\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * \warning: \f$\gamma_2\f$ is calculated before the call to \ref usipsu. Hence,
 * if \f$\beta_2\f$, \f$C_\mu\f$, \f$\kappa\f$ or \f$\sigma_{\omega 2}\f$ is
 * modified in \ref usipsu, \ref cs_turb_ckwgm2 must also be modified
 * in accordance.
 */
double cs_turb_ckwgm2;

/*!
 * Specific constant of k-omega SST.
 * Constant \f$a_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double cs_turb_ckwa1 = 0.31;

/*!
 * Constant \f$ c_1 \f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * Specific constant of k-omega SST.
 */
const double cs_turb_ckwc1 = 10.0;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double cs_turb_csab1 = 0.1355;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double cs_turb_csab2 = 0.622;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double cs_turb_csasig = 2.0/3.0;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double cs_turb_csav1 = 7.1;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csaw1;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double cs_turb_csaw2 = 0.3;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double cs_turb_csaw3 = 2.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
const double cs_turb_cssr1 = 1.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
const double cs_turb_cssr2 = 12.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
const double cs_turb_cssr3 = 1.0;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double cs_turb_ccaze2 = 1.83;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double cs_turb_ccazsc = 0.119;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double cs_turb_ccaza = 4.3;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double cs_turb_ccazb = 5.130;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double cs_turb_ccazc = 0.453;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double cs_turb_ccazd = 0.682;

/*!
 * Constant used in the definition of LES filtering diameter:
 * \f$ \delta = \text{xlesfl} . (\text{ales} . volume)^{\text{bles}} \f$.
 */
const double cs_turb_xlesfl = 2.0;

/*!
 * Constant used to define, for each cell \f$\Omega_i\f$, the width of
 * the (implicit) filter:
 * - \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
const double cs_turb_ales = 1.0;

/*!
 * Constant used to define, for each cell \f$Omega_i\f$, the width of
 * the (implicit) filter:
 * - \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
const double cs_turb_bles = 1.0/3.0;

/*!
 * Smagorinsky constant used in the Smagorinsky model for LES.
 * The sub-grid scale viscosity is calculated by
 * \f$\displaystyle\mu_{sg}=
 * \rho C_{smago}^2\bar{\Delta}^2\sqrt{2\bar{S}_{ij}\bar{S}_{ij}}\f$
 * where \f$\bar{\Delta}\f$ is the width of the filter
 * and \f$\bar{S}_{ij}\f$ the filtered strain rate.
 *
 * Useful if and only if \ref iturb = 40.
 * \note In theory Smagorinsky constant is 0.18.
 * For a channel, 0.065 value is rather taken.
 */
double cs_turb_csmago = 0.065;

/*!
 * Ratio between explicit and explicit filter width for a dynamic model.
 * Constant used to define, for each cell \f$\Omega_i\f$, the width of the
 * explicit filter used in the framework of the LES dynamic model:
 * \f$\widetilde{\overline{\Delta}}=xlesfd\overline{\Delta}\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
const double cs_turb_xlesfd = 1.5;

/*!
 * Maximum allowed value for the variable \f$C\f$ appearing in the LES dynamic
 * model (the "square" comes from the fact that the variable of the dynamic
 * model corresponds to the square of the constant of the Smagorinsky model).
 * Any larger value yielded by the calculation procedure of the dynamic model
 * will be clipped to \f$ smagmx^2\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
double cs_turb_smagmx;

/*!
 * Van Driest constant appearing in the van Driest damping function applied to
 * the Smagorinsky constant:
 * - \f$ (1-\exp^{(-y^+/cdries}) \f$.
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
const double cs_turb_cdries = 26.0;

/*!
 * Constant \f$a_1\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fa1 = 0.05;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fe2 = 1.85;

/*!
 * Constant \f$C_\mu\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fmu = 0.22;

/*!
 * Constant \f$C_1\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fc1 = 1.4;

/*!
 * Constant \f$C_2\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fc2 = 0.3;

/*!
 * Constant \f$C_T\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fct = 6.0;

/*!
 * Constant \f$C_L\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fcl = 0.25;

/*!
 * Constant \f$C_\eta\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cs_turb_cv2fet = 110.0;

/*!
 * Constant of the WALE LES method.
 */
const double cs_turb_cwale = 0.25;

/*!
 * Coefficient of turbulent AFM flow model.
 */
const double cs_turb_xiafm = 0.7;

/*!
 * Coefficient of turbulent AFM flow model.
 */
const double cs_turb_etaafm = 0.4;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double cs_turb_c1trit = 4.15;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double cs_turb_c2trit = 0.55;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double cs_turb_c3trit= 0.5;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double cs_turb_c4trit = 0.;

/*!
 * Constant of GGDH and AFM on the thermal scalar.
 */
const double cs_turb_cthafm = 0.236;

/*!
 * Constant of GGDH and AFM on the thermal scalar.
 */
const double cs_turb_cthdfm = 0.31;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur,
                             int     **nvarcl);

void
cs_f_turb_rans_model_get_pointers(int     **irccor,
                                  int     **itycor,
                                  int     **idirsm,
                                  int     **iclkep,
                                  int     **igrhok,
                                  int     **igrake,
                                  int     **igrari,
                                  int     **ikecou,
                                  int     **reinit_turb,
                                  int     **irijco,
                                  int     **irijnu,
                                  int     **irijrb,
                                  int     **irijec,
                                  int     **idifre,
                                  int     **iclsyr,
                                  int     **iclptr);

void
cs_f_turb_les_model_get_pointers(int     **idries,
                                 int     **ivrtex);

void
cs_f_turb_reference_values(double  **almax,
                           double  **uref,
                           double  **xlomlg);

void
cs_f_turb_complete_constants(void);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global turbulence model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iturb  --> pointer to cs_glob_turb_model->iturb
 *   itytur --> pointer to cs_glob_turb_model->itytur
 *   nvarcl --> pointer to cs_glob_turb_model->nvarcl
 *----------------------------------------------------------------------------*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur,
                             int     **nvarcl)
{
  *iturb  = &(_turb_model.iturb);
  *itytur = &(_turb_model.itytur);
  *nvarcl = &(_turb_model.nvarcl);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the RANS turbulence functions structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   irccor --> pointer to cs_glob_turb_rans_model->irccor
 *   itycor --> pointer to cs_glob_turb_rans_model->itycor
 *   idirsm --> pointer to cs_glob_turb_rans_model->idirsm
 *   iclkep --> pointer to cs_glob_turb_rans_model->iclkep
 *   igrhok --> pointer to cs_glob_turb_rans_model->igrhok
 *   igrake --> pointer to cs_glob_turb_rans_model->igrake
 *   igrari --> pointer to cs_glob_turb_rans_model->igrari
 *   ikecou --> pointer to cs_glob_turb_rans_model->ikecou
 *   reinit_turb --> pointer to cs_glob_turb_rans_model->reinit_turb
 *   irijco --> pointer to cs_glob_turb_rans_model->irijco
 *   irijnu --> pointer to cs_glob_turb_rans_model->irijnu
 *   irijrb --> pointer to cs_glob_turb_rans_model->irijrb
 *   irijec --> pointer to cs_glob_turb_rans_model->irijec
 *   idifre --> pointer to cs_glob_turb_rans_model->idifre
 *   iclsyr --> pointer to cs_glob_turb_rans_model->iclsyr
 *   iclptr --> pointer to cs_glob_turb_rans_model->iclptr
 *----------------------------------------------------------------------------*/

void
cs_f_turb_rans_model_get_pointers(int     **irccor,
                                  int     **itycor,
                                  int     **idirsm,
                                  int     **iclkep,
                                  int     **igrhok,
                                  int     **igrake,
                                  int     **igrari,
                                  int     **ikecou,
                                  int     **reinit_turb,
                                  int     **irijco,
                                  int     **irijnu,
                                  int     **irijrb,
                                  int     **irijec,
                                  int     **idifre,
                                  int     **iclsyr,
                                  int     **iclptr)
{
  *irccor = &(_turb_rans_model.irccor);
  *itycor = &(_turb_rans_model.itycor);
  *idirsm = &(_turb_rans_model.idirsm);
  *iclkep = &(_turb_rans_model.iclkep);
  *igrhok = &(_turb_rans_model.igrhok);
  *igrake = &(_turb_rans_model.igrake);
  *igrari = &(_turb_rans_model.igrari);
  *ikecou = &(_turb_rans_model.ikecou);
  *reinit_turb= &(_turb_rans_model.reinit_turb);
  *irijco = &(_turb_rans_model.irijco);
  *irijnu = &(_turb_rans_model.irijnu);
  *irijrb = &(_turb_rans_model.irijrb);
  *irijec = &(_turb_rans_model.irijec);
  *idifre = &(_turb_rans_model.idifre);
  *iclsyr = &(_turb_rans_model.iclsyr);
  *iclptr = &(_turb_rans_model.iclptr);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the LES turbulence model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   idries --> pointer to cs_glob_turb_les_model->idries
 *   ivrtex --> pointer to cs_glob_turb_les_model->ivrtex
 *----------------------------------------------------------------------------*/

void
cs_f_turb_les_model_get_pointers(int     **idries,
                                 int     **ivrtex)
{
  *idries = &(_turb_les_model.idries);
  *ivrtex = &(_turb_les_model.ivrtex);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the RANS turbulence functions structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   almax  --> pointer to cs_glob_turb_ref_values->almax
 *   uref   --> pointer to cs_glob_turb_ref_values->uref
 *   xlomlg --> pointer to cs_glob_turb_rans_model->xlomlg
 *----------------------------------------------------------------------------*/

void
cs_f_turb_reference_values(double  **almax,
                           double  **uref,
                           double  **xlomlg)
{
  *almax  = &(_turb_ref_values.almax);
  *uref   = &(_turb_ref_values.uref);
  *xlomlg = &(_turb_rans_model.xlomlg);
}

/*----------------------------------------------------------------------------
 * Initialize tubulence constants
 *----------------------------------------------------------------------------*/

void
cs_f_turb_complete_constants(void)
{
  cs_turb_dpow   = 1/(1.+cs_turb_bpow);
  cs_turb_cmu025 = pow(cs_turb_cmu,0.25);
  cs_turb_cstlog_alpha = exp(-cs_turb_xkappa
                             * (cs_turb_cstlog_rough - cs_turb_cstlog));


  if (   cs_glob_turb_model->iturb == 30
      || cs_glob_turb_model->iturb == 31)
    cs_turb_sigmae = 1.22;
  else if (cs_glob_turb_model->iturb == 32)
    cs_turb_sigmae = 1.15;
  else
    cs_turb_sigmae = 1.30;

  if (cs_glob_turb_model->iturb == 32)
    cs_turb_csrij = 0.21;
  else
    cs_turb_csrij = 0.22;

  double xkappa2 = cs_turb_xkappa*cs_turb_xkappa;
  cs_turb_ckwgm1 =   cs_turb_ckwbt1/cs_turb_cmu
                   - xkappa2/(cs_turb_ckwsw1*sqrt(cs_turb_cmu));
  cs_turb_ckwgm2 =   cs_turb_ckwbt2/cs_turb_cmu
                   - xkappa2/(cs_turb_ckwsw2*sqrt(cs_turb_cmu));
  cs_turb_csaw1 =   cs_turb_csab1/xkappa2
                  + 1./cs_turb_csasig*(1. + cs_turb_csab2);
  cs_turb_smagmx = 10.*cs_turb_csmago;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_turb_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_model_t *
cs_get_glob_turb_model(void)
{
  return &_turb_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_turb_ref_values
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_ref_values_t *
cs_get_glob_turb_ref_values(void)
{
  return &_turb_ref_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_turb_rans_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_rans_model_t *
cs_get_glob_turb_rans_model(void)
{
  return &_turb_rans_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_turb_les_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_turb_les_model_t *
cs_get_glob_turb_les_model(void)
{
  return &_turb_les_model;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

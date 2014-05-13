/*============================================================================
 * Base turbulence model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
  \var  cs_turb_rans_model_t::almax
        characteristic macroscopic length of the domain, used for the
        initialisation of the turbulence and the potential clipping (with
        \ref iclkep=1)
        - Negative value: not initialised (the code then uses the cubic root of
          the domain volume).

        Useful if and only if \ref turb = 20, 21, 30, 31, 50 or 60
        (RANS models).
  \var  cs_turb_rans_model_t::uref
        characteristic flow velocity, used for the initialisation of the
        turbulence
        - Negative value: not initialised.

        Useful if and only if \ref iturb= 20, 21, 30, 31, 50 or 60 (RANS model)
        and the turbulence is not initialised somewhere else (restart file or
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

static cs_turb_model_t  _turb_model = {-999, -999, 0};

const cs_turb_model_t  *cs_glob_turb_model = &_turb_model;

/* RANS turbulence model structure and associated pointer */

static cs_turb_rans_model_t  _turb_rans_model = {0, -999, 0, 0, 1, 1, -999, 0, 0,
                                                 0, 1, 1, 0, -999, -1e13, -1e13};

const cs_turb_rans_model_t  *cs_glob_turb_rans_model = &_turb_rans_model;

/* LES turbulence model structure and associated pointer */

static cs_turb_les_model_t  _turb_les_model = {-1, 0};

const cs_turb_les_model_t  *cs_glob_turb_les_model = &_turb_les_model;

/*! \endcond (end ignore by Doxygen) */

/*!
 * Karman constant. (= 0.42)
 *
 * Useful if and only if \ref iturb >= 10.
 *  (mixing length, \f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$,
 * LES, v2f or \f$k-\omega\f$).
 */
const double xkappa = 0.42;

/*!
 * Constant of logarithmic law function:
 * \f$ \dfrac{1}{\kappa} \ln(y^+) + cstlog \f$
 * (\f$ cstlog = 5.2 \f$).
 *
 * Constant of the logarithmic wall function.
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
const double cstlog = 5.2;

/*! Werner and Wengle coefficient */
const double apow = 8.3;

/*! Werner and Wengle coefficient */
const double bpow = 1.0/7.0;

/*! Werner and Wengle coefficient */
double dpow;

/*!
 * Constant \f$C_\mu\f$ for all the RANS turbulence models except for the
 * v2f model (see \ref cv2fmu for the value of \f$C_\mu\f$ in case of v2f
 * modelling). Useful if and only if \ref iturb = 20, 21, 30, 31 or 60
 * (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or \f$k-\omega\f$).
 */
const double cmu = 0.09;

/*! \f$ C_\mu^\frac{1}{4} \f$ */
double cmu025;

/*!
 * Constant \f$C_{\varepsilon 1}\f$ for all the RANS turbulence models except
 * for the v2f and the \f$k-\omega\f$ models.
 * Useful if and only if \ref iturb= 20, 21, 30 or 31 (\f$k-\varepsilon\f$
 * or \f$R_{ij}-\varepsilon\f$).
 */
const double ce1 = 1.44;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the \f$k-\varepsilon\f$ and
 * \f$R_{ij}-\varepsilon\f$ LRR models.
 * Useful if and only if {\tt iturb}= 20, 21 or 30
 * (\f$k-\varepsilon\f$ or \f$R_{ij}-\varepsilon\f$ LRR).
 */
const double ce2 = 1.92;

/*!
 * Coefficient of interfacial coefficient in k-eps, used in Lagrange treatment.
 *
 * Constant \f$C_{\varepsilon 4}\f$ for the interfacial term (Lagrangian module)
 * in case of two-way coupling. Useful in case of Lagrangian modelling,
 * in \f$k-\varepsilon\f$ and \f$R_{ij}-\varepsilon\f$ with two-way coupling.
 */
const double ce4 = 1.20;

/*!
 * Prandtl number for \f$k\f$ with \f$k-\varepsilon\f$ and v2f models.
 * Useful if and only if \ref iturb=20, 21 or 50 (\f$k-\varepsilon\f$ or v2f).
 */
const double sigmak = 1.0;

/*!
 * Prandtl number for \f$\varepsilon\f$.
 * Useful if and only if \ref iturb= 20, 21, 30, 31 or 50
 * (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or v2f).
 */
double sigmae;

/*!
 * Constant \f$C_1\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double crij1 = 1.80;

/*
 * Constant \f$C_2\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double crij2 = 0.60;

/*!
 * Constant \f$C_3\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double crij3 = 0.55;

/*!
 * Constant \f$C_1^\prime\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model,
 * corresponding to the wall echo terms.
 * Useful if and only if \ref iturb=30 and \ref irijec=1
 * (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double crijp1 = 0.50;

/*!
 * Constant \f$C_2^\prime\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model,
 * corresponding to the wall echo terms.
 * Useful if and only if \ref iturb=30 and \ref irijec=1
 * (\f$R_{ij}-\varepsilon\f$ LRR).
 */
const double crijp2 = 0.30;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssge2 = 1.83;

/*!
 * Constant \f$C_{s1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgs1 = 1.70;

/*!
 * Constant \f$C_{s2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgs2 = -1.05;

/*!
 * Constant \f$C_{r1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgr1 = 0.90;

/*!
 * Constant \f$C_{r2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgr2 = 0.80;

/*!
 * Constant \f$C_{r3}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgr3 = 0.65;

/*!
 * constant \f$C_{r4}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgr4 = 0.625;

/*!
 * Constant \f$C_{r1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
const double cssgr5 = 0.20;

/*! Constant of the Rij-epsilon EBRSM. */
const double cebms1 = 1.70;

/*! Constant of the Rij-epsilon EBRSM. */
const double cebms2 = 0.;

const double cebmr1 = 0.90;
const double cebmr2 = 0.80;
const double cebmr3 = 0.65;
const double cebmr4 = 0.625;
const double cebmr5 = 0.20;
const double cebmr6 = 0.6;

/*!
 * Constant \f$C_s\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double csrij;

/*! Constant of the Rij-epsilon EBRSM. */
const double cebme2 = 1.83;

/*! Constant of the Rij-epsilon EBRSM. */
const double cebmmu = 0.22;

/*! Constant of the Rij-epsilon EBRSM. */
const double xcl = 0.122;

/*! Constant in the expression of Ce1' for the Rij-epsilon EBRSM. */
const double xa1 = 0.1;

/*! Constant of the Rij-epsilon EBRSM. */
const double xct = 6.0;

/*! Constant of the Rij-epsilon EBRSM. */
const double xceta = 80.0;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpale1 = 1.44;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpale2 = 1.83;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpale3 = 2.3;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpale4 = 0.4;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalse = 1.5;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalmu = 0.22;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalc1 = 1.7;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalc2 = 0.9;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalct = 4.0;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalcl = 0.164;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
const double cpalet = 75.0;

/*!
 * Constant \f$\sigma_{k1}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60.
 */
const double ckwsk1 = 1.0/0.85;

/*!
 * Constant \f$\sigma_{k2}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60.
 */
const double ckwsk2 = 1.0;

/*!
 * Constant \f$\sigma_{\omega 1}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double ckwsw1 = 2.0;

/*!
 * Constant \f$\sigma_{\omega 2}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double ckwsw2 = 1.0/0.856;

/*!
 * Constant \f$\beta_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double ckwbt1 = 0.075;

/*!
 * Constant \f$\beta_2\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double ckwbt2 = 0.0828;


/*!
 * \f$\frac{\beta_1}{C_\mu}-\frac{\kappa^2}{\sqrt{C_\mu}\sigma_{\omega 1}}\f$.
 * Constant \f$\gamma_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * \warning: \f$\gamma_1\f$ is calculated before the call to
 * \ref usipsu. Hence, if \f$\beta_1\f$, \f$C_\mu\f$, \f$\kappa\f$ or
 * \f$\sigma_{\omega 1}\f$ is modified in \ref usipsu,
 * \ref ckwgm1 must also be modified in accordance.
 */
double ckwgm1;

/*!
 * \f$\frac{\beta_2}{C_\mu}-\frac{\kappa^2}{\sqrt{C_\mu}\sigma_{\omega 2}}\f$.
 * Constant \f$\gamma_2\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * \warning: \f$\gamma_2\f$ is calculated before the call to \ref usipsu. Hence,
 * if \f$\beta_2\f$, \f$C_\mu\f$, \f$\kappa\f$ or \f$\sigma_{\omega 2}\f$ is
 * modified in \ref usipsu, \ref ckwgm2 must also be modified in accordance.
 */
double ckwgm2;

/*!
 * Specific constant of k-omega SST.
 * Constant \f$a_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
const double ckwa1 = 0.31;

/*!
 * Constant \f$ c_1 \f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * Specific constant of k-omega SST.
 */
const double ckwc1 = 10.0;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double csab1 = 0.1355;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double csab2 = 0.622;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double csasig = 2.0/3.0;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double csav1 = 7.1;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double csaw1;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double csaw2 = 0.3;

/*!
 * Specific constant of Spalart-Allmaras.
 */
const double csaw3 = 2.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
const double cssr1 = 1.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
const double cssr2 = 12.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
const double cssr3 = 1.0;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double ccaze2 = 1.83;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double ccazsc = 0.119;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double ccaza = 4.3;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double ccazb = 5.130;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double ccazc = 0.453;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
const double ccazd = 0.682;

/*!
 * Constant used in the definition of LES filtering diameter:
 * \f$ \delta = \text{xlesfl} . (\text{ales} . volume)^{\text{bles}} \f$.
 */
const double xlesfl = 2.0;

/*!
 * Constant used to define, for each cell \f$\Omega_i\f$, the width of
 * the (implicit) filter:
 * - \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
const double ales = 1.0;

/*!
 * Constant used to define, for each cell $\Omega_i$, the width of
 * the (implicit) filter:
 * - \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
const double bles = 1.0/3.0;

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
const double csmago = 0.065;

/*!
 * Ratio between explicit and explicit filter width for a dynamic model.
 * Constant used to define, for each cell \f$\Omega_i\f$, the width of the
 * explicit filter used in the framework of the LES dynamic model:
 * \f$\widetilde{\overline{\Delta}}=xlesfd\overline{\Delta}\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
const double xlesfd = 1.5;

/*!
 * Maximum allowed value for the variable \f$C\f$ appearing in the LES dynamic
 * model (the "square" comes from the fact that the variable of the dynamic
 * model corresponds to the square of the constant of the Smagorinsky model).
 * Any larger value yielded by the calculation procedure of the dynamic model
 * will be clipped to \f$ smagmx^2\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
double smagmx;

/*!
 * Van Driest constant appearing in the van Driest damping function applied to
 * the Smagorinsky constant:
 * - \f$ (1-\exp^{(-y^+/cdries}) \f$.
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
const double cdries = 26.0;

/*!
 * Constant \f$a_1\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fa1 = 0.05;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fe2 = 1.85;

/*!
 * Constant \f$C_\mu\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fmu = 0.22;

/*!
 * Constant \f$C_1\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fc1 = 1.4;

/*!
 * Constant \f$C_2\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fc2 = 0.3;

/*!
 * Constant \f$C_T\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fct = 6.0;

/*!
 * Constant \f$C_L\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fcl = 0.25;

/*!
 * Constant \f$C_\eta\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
const double cv2fet = 110.0;

/*!
 * Constant of the WALE LES method.
 */
const double cwale = 0.25;

/*!
 * Coefficient of turbulent AFM flow model.
 */
const double xiafm = 0.7;

/*!
 * Coefficient of turbulent AFM flow model.
 */
const double etaafm = 0.4;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double c1trit = 4.15;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double c2trit = 0.55;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double c3trit= 0.5;

/*!
 * Coefficient of turbulent DFM flow model.
 */
const double c4trit = 0.;

/*!
 * Constant of GGDH and AFM on the thermal scalar.
 */
const double cthafm = 0.236;

/*!
 * Constant of GGDH and AFM on the thermal scalar.
 */
const double cthdfm = 0.31;

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

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
 *   almax  --> pointer to cs_glob_turb_rans_model->almax
 *   uref   --> pointer to cs_glob_turb_rans_model->uref
 *   xlomlg --> pointer to cs_glob_turb_rans_model->xlomlg
 *----------------------------------------------------------------------------*/

void
cs_f_turb_reference_values(double  **almax,
                           double  **uref,
                           double  **xlomlg)
{
  *almax  = &(_turb_rans_model.almax);
  *uref   = &(_turb_rans_model.uref);
  *xlomlg = &(_turb_rans_model.xlomlg);
}

/*----------------------------------------------------------------------------
 * Initialise tubulence constants
 *----------------------------------------------------------------------------*/

void
cs_f_turb_complete_constants(void)
{
  dpow   = 1/(1.+bpow);
  cmu025 = pow(cmu,0.25);

  if (   cs_glob_turb_model->iturb == 30
      || cs_glob_turb_model->iturb == 31)
    sigmae = 1.22;
  else if (cs_glob_turb_model->iturb == 32)
    sigmae = 1.15;
  else
    sigmae = 1.30;

  if (cs_glob_turb_model->iturb == 32)
    csrij = 0.21;
  else
    csrij = 0.22;

  double xkappa2 = xkappa*xkappa;
  ckwgm1 = ckwbt1/cmu - xkappa2/(ckwsw1*sqrt(cmu));
  ckwgm2 = ckwbt2/cmu - xkappa2/(ckwsw2*sqrt(cmu));
  csaw1 = csab1/xkappa2 + 1./csasig*(1. + csab2);
  smagmx = 10.*csmago;
}

/*! \endcond (end ignore by Doxygen) */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

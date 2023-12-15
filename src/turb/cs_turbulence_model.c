/*============================================================================
 * Base turbulence model data.
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

#include "cs_assert.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_function.h"
#include "cs_gradient.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_wall_functions.h"

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
        - CS_TURB_NONE: no turbulence model (laminar flow)
        - CS_TURB_MIXING_LENGTH: mixing length model
        - CS_TURB_K_EPSILON: standard \f$ k-\varepsilon \f$ model
        - CS_TURB_K_EPSILON_LIN_PROD: \f$ k-\varepsilon \f$ model with
                                       Linear Production (LP) correction
        - CS_TURB_K_EPSILON_LS: Launder-Sharma \f$ k-\varepsilon \f$ model
        - CS_TURB_K_EPSILON_QUAD: Baglietto et al. quadratic
                                   \f$ k-\varepsilon \f$ model
        - CS_TURB_RIJ_EPSILON_LRR: \f$ R_{ij}-\epsilon \f$ (LRR)
        - CS_TURB_RIJ_EPSILON_SSG: \f$ R_{ij}-\epsilon \f$ (SSG)
        - CS_TURB_RIJ_EPSILON_EBRSM: \f$ R_{ij}-\epsilon \f$ (EBRSM)
        - CS_TURB_LES_SMAGO_CONST: LES (constant Smagorinsky model)
        - CS_TURB_LES_SMAGO_DYN: LES ("classical" dynamic Smagorisky model)
        - CS_TURB_LES_WALE: LES (WALE)
        - CS_TURB_V2F_PHI: v2f phi-model
        - CS_TURB_V2F_BL_V2K: v2f \f$ BL-v^2-k \f$
        - CS_TURB_K_OMEGA: \f$ k-\omega \f$ SST
        - CS_TURB_SPALART_ALLMARAS: Spalart-Allmaras model
  \var  cs_turb_model_t::itytur
        class of turbulence model (integer value iturb/10, deprecated)
  \var  cs_turb_model_t::hybrid_turb
        Type of hybrid turbulence model
        - 0: No model
        - 1: Detached Eddy Simulation
        - 2: Delayed Detached Eddy Simulation
        - 3: Scale Adaptive Model (Menter et al.)
        - 4. Hybrid Temporal LES
  \var  cs_turb_model_t::type
        Type of modelling
        - CS_TURB_NONE: No model
        - CS_TURB_RANS: RANS
        - CS_TURB_LES: LES
        - CS_TURB_HYBRID: Hybrid RANS-LES
  \var  cs_turb_model_t::order
        Order of the turbulence model:
        - CS_TURB_ALGEBRAIC: 0th order algebraic model
        - CS_TURB_FIRST_ORDER: 1st order Eddy Viscosity
                               type models
        - CS_TURB_SECOND_ORDER: 2nd order Differential
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
        - 1: Cazalbou correction (default when \ref irccor = 1 and
        \ref cs_turb_model_t::itytur "itytur" = 2 or 5)
        - 2: Spalart-Shur correction (default when \ref irccor = 1 and
        \ref iturb = 60 or 70)
  \var  cs_turb_rans_model_t::idirsm
        turbulent diffusion model for second moment closure
        - 0: scalar diffusivity (Shir model, default model)
        - 1: tensorial diffusivity (Daly and Harlow model)
  \var  cs_turb_rans_model_t::iclkep
        Indicates the clipping method used for \f$k\f$ and
        \f$\varepsilon\f$, for the \f$k-\epsilon\f$ and v2f models\n
        - 0: clipping in absolute value
        - 1: coupled clipping based on physical relationships\n
        Useful if and only if \ref iturb = 20, 21 or 50 (\f$k-\epsilon\f$ and
        v2f models). The results obtained with the method corresponding to
        \ref iclkep =1 showed in some cases a substantial sensitivity to the
        values of the length scale \ref almax.\n
        The option \ref iclkep = 1 is therefore not recommended, and,
        if chosen, must be used cautiously.
  \var  cs_turb_rans_model_t::igrhok
        Indicates if the term \f$\frac{2}{3}\grad \rho k\f$
        is taken into account in the velocity equation.
        - 1: true
        - 0: false in the velocity
        Useful if and only if \ref iturb = 20, 21, 50 or 60.\n
        This term may generate non-physical velocities at the wall.
        When it is not explicitly taken into account, it is
        implicitly included into the pressure.
  \var  cs_turb_rans_model_t::igrake
        Indicates if the terms related to gravity are taken
        into account in the equations of \f$k-\epsilon\f$.\n
        - 1: true (default if \f$ \rho \f$ is variable)
        - 0: false
        Useful if and only if \ref iturb = 20, 21, 50 or 60 and
        (\ref cs_physical_constants_t::gravity "gravity")
        \f$\ne\f$ (0,0,0) and the density is not uniform.
  \var  cs_turb_rans_model_t::igrari
        Indicates if the terms related to gravity are taken
        into account in the equations of \f$R_{ij}-\epsilon\f$.\n
        - 1: true (default if \f$ \rho \f$ is variable)
        - 0: false
        Useful if and only if \ref iturb = 30 or 31 and
        (\ref cs_physical_constants_t::gravity "gravity") \f$\ne\f$
        (0,0,0) (\f$R_{ij}-\epsilon\f$ model with gravity) and the
        density is not uniform.
  \var  cs_turb_rans_model_t::ikecou
        Indicates if the coupling of the source terms of
        \f$k\f$ and \f$\epsilon\f$ or \f$k\f$ and \f$\omega\f$
        is taken into account or not.
        - 1: true,
        - 0: false\n
        If \ref ikecou = 0 in \f$k-\epsilon\f$ model, the term
        in \f$\epsilon\f$ in the equation of \f$k\f$ is made implicit.\n
        \ref ikecou is initialised to 0 if \ref iturb = 21 or 60, and
        to 1 if \ref iturb = 20.\n
        \ref ikecou = 1 is forbidden when using the v2f model (\ref iturb = 50).\n
        Useful if and only if \ref iturb = 20, 21 or 60 (\f$k-\epsilon\f$ and
        \f$k-\omega\f$ models)
  \var  cs_turb_rans_model_t::reinit_turb
        Advanced re-init for EBRSM and k-omega models
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijco
        coupled solving of Rij
        - 1: true (default)
        - 0: false
  \var  cs_turb_rans_model_t::irijnu
        pseudo eddy viscosity in the matrix of momentum equation to partially
        implicit \f$ \divv \left( \rho \tens{R} \right) \f$
        - 1: true
        - 0: false (default)
        The goal is to improve the stability of the calculation.
        The usefulness of \ref irijnu = 1 has however not been
        clearly demonstrated.\n Since the system is solved in
        incremental form, this extra turbulent viscosity does
        not change the final solution for steady flows. However,
        for unsteady flows, the parameter \ref cs_var_cal_opt_t::nswrsm "nswrsm"
        should be increased.\n Useful if and only if \ref iturb = 30
        or 31 (\f$R_{ij}-\epsilon\f$ model).
  \var  cs_turb_rans_model_t::irijrb
        accurate treatment of \f$ \tens{R} \f$ at the boundary
        (see \ref cs_boundary_condition_set_coeffs)
        - 1: true
        - 0: false (default)
  \var  cs_turb_rans_model_t::irijec
        Indicates if the wall echo terms in
        \f$R_{ij}-\epsilon\f$ LRR model are taken into account:
        - 1: true,
        - 0: false (default)\n
        Useful if and only if \ref iturb = 30 (\f$R_{ij}-\epsilon\f$
        LRR).\n It is not recommended to take these terms into account:
        they have an influence only near the walls, their expression is hardly
        justifiable according to some authors and, in the configurations
        studied with code_saturne, they did not bring any improvement in the
        results.\n
        In addition, their use induces an increase in the calculation time.\n
        The wall echo terms imply the calculation of the distance to the wall
        for every cell in the domain. See \ref optcal::icdpar "icdpar" for
        potential
        restrictions due to this.
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
  \var  cs_turb_rans_model_t::ikwcln
        Wall boundary condition on omega in k-omega SST
        0: Deprecated Neumann boundary condition
        1: Dirichlet boundary condition consistent with Menter's
        original model: w_wall = 60*nu/(beta*d**2)
  \var  cs_turb_ref_values_t::almax
        characteristic macroscopic length of the domain, used for the
        initialization of the turbulence and the potential clipping (with
        \ref cs_turb_rans_model_t::iclkep "iclkep"= 1)
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
        Activates or the van Driest wall-damping for the
        Smagorinsky constant (the Smagorinsky constant
        is multiplied by the damping function
        \f$1-e^{-y^+/cdries}\f$, where \f$y^+\f$
        designates the non-dimensional distance to the
        nearest wall).
           - 1: true
           - 0: false
        The default value is 1 for the Smagorinsky model
        and 0 for the dynamic model.\n The van Driest
        wall-damping requires the knowledge of the
        distance to the nearest wall for each cell
        in the domain. Refer to keyword \ref optcal::icdpar "icdpar"
        for potential limitations.\n
        Useful if and only if \ref iturb = 40 or 41
*/

/*----------------------------------------------------------------------------*/

/*! \struct cs_turb_hybrid_model_t

  \brief Hybrid turbulence model descriptor.

  Members of this turbulence model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_turb_hybrid_model_t::iicc
        Applied or not the Internal Consistency
        Constraint (ICC) for the HTLES model,
        in order to recover the correct RANS
        behavior when the energy ratio is forced
        to one in the RANS region:
          - 1: True (default)
          - 0: False
        Useful if and only if \ref hybrid_turb=4

  \var  cs_turb_hybrid_model_t::ishield
        Applied or not the two-fold shielding
        function (\f$f_s(\xi_K,\xi_D)\f$ of HTLES,
        to properly control the RANS-to-LES
        transition in the vicinity of the wall:
          - 1: True (default)
          - 0: False
        Useful if and only if \ref hybrid_turb=4
*/

/*----------------------------------------------------------------------------*/

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
  .hybrid_turb = 0,
  .type = -1,
  .order = -1
};

const cs_turb_model_t  *cs_glob_turb_model = NULL;

/* Reference values for turbulence structure and associated pointer */

static cs_turb_ref_values_t
_turb_ref_values =
{
  .almax =  -999,
  .uref  = -1e13
};

const cs_turb_ref_values_t  *cs_glob_turb_ref_values = &_turb_ref_values;

/* RANS turbulence model structure and associated pointer */

static cs_turb_rans_model_t
_turb_rans_model =
{
  .irccor     =    0,
  .itycor     = -999,
  .idirsm     =    0,
  .iclkep     =    0,
  .igrhok     =    0,
  .igrake     =    1,
  .igrari     =    1,
  .ikecou     =    0,
  .reinit_turb=    1,
  .irijco     =    1, /* Coupled version of DRSM models */
  .irijnu     =    0,
  .irijrb     =    0,
  .irijec     =    0,
  .idifre     =    1,
  .iclsyr     =    1,
  .iclptr     =    0,
  .ikwcln     =    1,
  .xlomlg     = -1e13
};

const cs_turb_rans_model_t  *cs_glob_turb_rans_model = &_turb_rans_model;

/* LES turbulence model structure and associated pointer */

static cs_turb_les_model_t  _turb_les_model =
{
  .idries = -1,
};

const cs_turb_les_model_t  *cs_glob_turb_les_model = &_turb_les_model;

/* Hybrid turbulence model structure and associated pointer */

static cs_turb_hybrid_model_t  _turb_hybrid_model =
{
  .iicc    = 1,
  .ishield = 1,
  .n_iter_mean = -1,
  .time_mean = -1.
};

const cs_turb_hybrid_model_t  *cs_glob_turb_hybrid_model = &_turb_hybrid_model;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*!
 * Karman constant. (= 0.42)
 *
 * Useful if and only if \ref iturb >= 10.
 *  (mixing length, \f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$,
 * LES, v2f or \f$k-\omega\f$).
 */
double cs_turb_xkappa = 0.42;

/*!
 * Van Driest constant. (= 25.6)
 *
 * Useful if and only if \ref cs_wall_functions_t::iwallf
 * "cs_glob_wall_functions::iwallf" = CS_WALL_F_2SCALES_VDRIEST.
 *  (Two scales log law at the wall using Van Driest mixing length expression).
 */
double cs_turb_vdriest = 25.6;

/*!
 * Constant of logarithmic smooth law function:
 * \f$ \dfrac{1}{\kappa} \ln(y^+) + cstlog \f$
 * (\f$ cstlog = 5.2 \f$).
 *
 * Constant of the logarithmic wall function.
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
double cs_turb_cstlog = 5.2;

/*!
 * Constant of logarithmic rough law function:
 * \f$ \dfrac{1}{\kappa} \ln(y/\xi) + cstlog_{rough} \f$
 * (\f$ cstlog_{rough} = 8.5 \f$).
 *
 * Constant of the logarithmic wall function.
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
double cs_turb_cstlog_rough = 8.5;

/*!
 * Constant \f$ \alpha \f$ for logarithmic law function switching from rough
 * to smooth:
 * \f$ \dfrac{1}{\kappa} \ln(y u_k/(\nu + \alpha \xi u_k)) + cstlog \f$
 * (\f$ \alpha = \exp \left( -\kappa (8.5 - 5.2) \right) \f$).
 *
 * Useful if and only if \ref iturb >= 10 (mixing length, \f$k-\varepsilon\f$,
 * \f$R_{ij}-\varepsilon\f$, LES, v2f or \f$k-\omega\f$).
 */
double cs_turb_cstlog_alpha;

/*! Werner and Wengle coefficient */
double cs_turb_apow = 8.3;

/*! Werner and Wengle coefficient */
double cs_turb_bpow = 1.0/7.0;

/*! Werner and Wengle coefficient */
double cs_turb_dpow = -1.;

/*!
 * Constant \f$C_\mu\f$ for all the RANS turbulence models.
 * Warning: different value for v2f models. Useful only for
 * RANS models
 * (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or \f$k-\omega\f$).
 */
double cs_turb_cmu = 0.09;

/*! \f$ C_\mu^\frac{1}{4} \f$ */
double cs_turb_cmu025 = 0.547722557; /* computed more precisely later */

/*!
 * Constant \f$C_{\varepsilon 1}\f$ for all the RANS turbulence models except
 * for the v2f and the \f$k-\omega\f$ models.
 * Useful if and only if \ref iturb= 20, 21, 30 or 31 (\f$k-\varepsilon\f$
 * or \f$R_{ij}-\varepsilon\f$).
 */
double cs_turb_ce1 = 1.44;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the \f$k-\varepsilon\f$ and
 * \f$R_{ij}-\varepsilon\f$ LRR models.
 * Useful if and only if \ref iturb = 20, 21 or 30
 * (\f$k-\varepsilon\f$ or \f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_ce2 = 1.92;

/*!
 * Coefficient of interfacial coefficient in k-eps, used in Lagrange treatment.
 *
 * Constant \f$C_{\varepsilon 4}\f$ for the interfacial term (Lagrangian module)
 * in case of two-way coupling. Useful in case of Lagrangian modelling,
 * in \f$k-\varepsilon\f$ and \f$R_{ij}-\varepsilon\f$ with two-way coupling.
 */
double cs_turb_ce4 = 1.20;

/*!
 * Constant \f$C_1\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_crij1 = 1.80;

/*!
 * Constant \f$C_2\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_crij2 = 0.60;

/*!
 * Rotta constant \f$C_0\f$ for the \f$R_{ij}-\varepsilon\f$ model.
 * Useful for the Lagrangian model. The value is set from \f$C_1\f$
 * if and only if \ref iturb=CS_TURB_RIJ_EPSILON_LRR
 * (\f$R_{ij}-\varepsilon\f$ LRR) and \f$C_2=0\f$.
 */
double cs_turb_crij_c0 = 3.5;

/*!
 * Constant \f$C_3\f$ for the \f$R_{ij}-\varepsilon\f$ models.
 * Value is 0.55 for SSG and LRR, 0.6 for EBRSM.
 */
double cs_turb_crij3 = 0.55;

/*!
 * Constant \f$C_1^\prime\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model,
 * corresponding to the wall echo terms.
 * Useful if and only if \ref iturb=30 and \ref cs_turb_rans_model_t::irijec
 * "cs_turb_rans_model_t::irijec"=1
 * (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_crijp1 = 0.50;

/*!
 * Constant \f$C_2^\prime\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model,
 * corresponding to the wall echo terms.
 * Useful if and only if \ref iturb=30 and \ref cs_turb_rans_model_t::irijec=1
 * (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_crijp2 = 0.30;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssge2 = 1.83;

/*!
 * Constant \f$C_{s1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgs1 = 1.70;

/*!
 * Constant \f$C_{s2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgs2 = -1.05;

/*!
 * Constant \f$C_{r1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgr1 = 0.90;

/*!
 * Constant \f$C_{r2}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgr2 = 0.80;

/*!
 * Constant \f$C_{r3}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgr3 = 0.65;

/*!
 * constant \f$C_{r4}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgr4 = 0.625;

/*!
 * Constant \f$C_{r1}\f$ for the \f$R_{ij}-\varepsilon\f$ SSG model.
 * Useful if and only if \ref iturb=31 (\f$R_{ij}-\varepsilon\f$ SSG).
 */
double cs_turb_cssgr5 = 0.20;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_cebms1 = 1.70;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_cebms2 = 0.;

double cs_turb_cebmr1 = 0.90;
double cs_turb_cebmr2 = 0.80;
double cs_turb_cebmr3 = 0.65;
double cs_turb_cebmr4 = 0.625;
double cs_turb_cebmr5 = 0.20;

/*!
 * Constant \f$C_s\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_csrij;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_cebme2 = 1.83;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_cebmmu = 0.22;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_xcl = 0.122;

/*! Constant in the expression of Ce1' for the Rij-epsilon EBRSM. */
double cs_turb_xa1 = 0.1;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_xct = 6.0;

/*! Constant of the Rij-epsilon EBRSM. */
double cs_turb_xceta = 80.0;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpale1 = 1.44;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpale2 = 1.83;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpale3 = 2.3;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpale4 = 0.4;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpalc1 = 1.7;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpalc2 = 0.9;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpalct = 4.0;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpalcl = 0.164;

/*! Specific constant of v2f "BL-v2k" (or phi-alpha). */
double cs_turb_cpalet = 75.0;

/*!
 * Constant \f$\sigma_{k1}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60.
 */
double cs_turb_ckwsk1 = 1.0/0.85;

/*!
 * Constant \f$\sigma_{k2}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60.
 */
double cs_turb_ckwsk2 = 1.0;

/*!
 * Constant \f$\sigma_{\omega 1}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
double cs_turb_ckwsw1 = 2.0;

/*!
 * Constant \f$\sigma_{\omega 2}\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
double cs_turb_ckwsw2 = 1.0/0.856;

/*!
 * Constant \f$\beta_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
double cs_turb_ckwbt1 = 0.075;

/*!
 * Constant \f$\beta_2\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
double cs_turb_ckwbt2 = 0.0828;

/*!
 * \f$\frac{\beta_1}{C_\mu}-\frac{\kappa^2}{\sqrt{C_\mu}\sigma_{\omega 1}}\f$.
 * Constant \f$\gamma_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * \warning: \f$\gamma_1\f$ is calculated before the call to
 * \ref usipsu. Hence, if \f$\beta_1\f$, \f$C_\mu\f$, \f$\kappa\f$ or
 * \f$\sigma_{\omega 1}\f$ is modified in \ref usipsu,
 * \ref cs_turb_ckwgm1 must also be modified in accordance.
 */
double cs_turb_ckwgm1 = -1.;

/*!
 * \f$\frac{\beta_2}{C_\mu}-\frac{\kappa^2}{\sqrt{C_\mu}\sigma_{\omega 2}}\f$.
 * Constant \f$\gamma_2\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * \warning: \f$\gamma_2\f$ is calculated before the call to \ref usipsu. Hence,
 * if \f$\beta_2\f$, \f$C_\mu\f$, \f$\kappa\f$ or \f$\sigma_{\omega 2}\f$ is
 * modified in \ref usipsu, \ref cs_turb_ckwgm2 must also be modified
 * in accordance.
 */
double cs_turb_ckwgm2 = -1.;

/*!
 * Specific constant of k-omega SST.
 * Constant \f$a_1\f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 */
double cs_turb_ckwa1 = 0.31;

/*!
 * Constant \f$ c_1 \f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST).
 * Specific constant of k-omega SST.
 */
double cs_turb_ckwc1 = 10.0;

/*!
 * Constant \f$ C_{DDES} \f$ for the \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST) and hybrid_turb=1.
 */
double cs_turb_cddes = 0.65;

/*!
 * Constant \f$ C_{SAS}\f$ for the hybrid \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST) and hybrid_turb=3.
 */
double cs_turb_csas = 0.11;

/*! constant \f$ C_{DDES}\f$ for the hybrid \f$k-\omega\f$ SST model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST) and hybrid_turb=3.
 */
double cs_turb_csas_eta2 = 3.51;

/*!
 * Constant \f$ \beta_0 \f$ for the HTLES model.
 * Useful if and only if \ref iturb=60 (\f$k-\omega\f$ SST)
 * or if \ref iturb=51 (\f$BL-v^2-k\f$)
 * and hybrid_turb=4.
 */
double cs_turb_chtles_bt0 = 0.48;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csab1 = 0.1355;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csab2 = 0.622;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csasig = 2.0/3.0;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csav1 = 7.1;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csaw1 = -1.;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csaw2 = 0.3;

/*!
 * Specific constant of Spalart-Allmaras.
 */
double cs_turb_csaw3 = 2.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
double cs_turb_cssr1 = 1.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
double cs_turb_cssr2 = 2.0;

/*!
 * Constant of the Spalart-Shur rotation/curvature correction.
 */
double cs_turb_cssr3 = 1.0;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
double cs_turb_ccaze2 = 1.83;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
double cs_turb_ccazsc = 0.119;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
double cs_turb_ccaza = 4.3;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
double cs_turb_ccazb = 5.130;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
double cs_turb_ccazc = 0.453;

/*!
 * Constants of the Cazalbou rotation/curvature correction.
 */
double cs_turb_ccazd = 0.682;

/*!
 * Constant used in the definition of LES filtering diameter:
 * \f$ \delta = \text{xlesfl} . (\text{ales} . volume)^{\text{bles}}\f$
 * \ref cs_turb_xlesfl is a constant used to define, for
 * each cell \f$\Omega_i\f$, the width of the (implicit) filter:
 * \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$\n
 * Useful if and only if \ref iturb = 40 or 41
 */
double cs_turb_xlesfl = 2.0;

/*!
 * Constant used to define, for each cell \f$\Omega_i\f$, the width of
 * the (implicit) filter:
 * - \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
double cs_turb_ales = 1.0;

/*!
 * Constant used to define, for each cell \f$Omega_i\f$, the width of
 * the (implicit) filter:
 * - \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
double cs_turb_bles = 1.0/3.0;

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
double cs_turb_xlesfd = 1.5;

/*!
 * Maximum allowed value for the variable \f$C\f$ appearing in the LES dynamic
 * model.
 * Any larger value yielded by the calculation procedure of the dynamic model
 * will be clipped to \f$ smagmx\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
double cs_turb_csmago_max = -1.;

/*!
 * Minimum allowed value for the variable \f$C\f$ appearing in the LES dynamic
 * model.
 * Any smaller value yielded by the calculation procedure of the dynamic model
 * will be clipped to \f$ smagmn\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
double cs_turb_csmago_min = 0.;

/*!
 * Van Driest constant appearing in the van Driest damping function applied to
 * the Smagorinsky constant:
 * - \f$ (1-\exp^{(-y^+/cdries}) \f$.
 *
 * Useful if and only if \ref iturb = 40 or 41.
 */
double cs_turb_cdries = 26.0;

/*!
 * Constant \f$a_1\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fa1 = 0.05;

/*!
 * Constant \f$C_{\varepsilon 2}\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fe2 = 1.85;

/*!
 * Constant \f$C_1\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fc1 = 1.4;

/*!
 * Constant \f$C_2\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fc2 = 0.3;

/*!
 * Constant \f$C_T\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fct = 6.0;

/*!
 * Constant \f$C_L\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fcl = 0.25;

/*!
 * Constant \f$C_\eta\f$ for the v2f \f$\varphi\f$-model.
 * Useful if and only if \ref iturb=50 (v2f \f$\varphi\f$-model).
 */
double cs_turb_cv2fet = 110.0;

/*!
 * Constants for the Baglietto et al. quadratic k-epsilon model.
 * Useful if and only if \ref iturb = CS_TURB_K_EPSILON_QUAD
 */
double cs_turb_cnl1 = 0.8;
double cs_turb_cnl2 = 11.;
double cs_turb_cnl3 = 4.5;
double cs_turb_cnl4 = 1e3;
double cs_turb_cnl5 = 1.;

/*!
 * Constant of the WALE LES method.
 */
double cs_turb_cwale = 0.25;

/*!
 * Coefficient of turbulent AFM flow model.
 */
double cs_turb_xiafm = 0.7;

/*!
 * Coefficient of turbulent AFM flow model.
 */
double cs_turb_etaafm = 0.4;

/*!
 * Coefficient of turbulent DFM flow model.
 */
double cs_turb_c1trit = 4.15;

/*!
 * Coefficient of turbulent DFM flow model.
 */
double cs_turb_c2trit = 0.55;

/*!
 * Coefficient of turbulent DFM flow model.
 */
double cs_turb_c3trit= 0.5;

/*!
 * Coefficient of turbulent DFM flow model.
 */
double cs_turb_c4trit = 0.;

/*!
 * Constant of GGDH and AFM on the thermal scalar.
 */
double cs_turb_cthafm = 0.236;

/*!
 * Constant of GGDH and AFM on the thermal scalar.
 */
double cs_turb_cthdfm = 0.31;
double cs_turb_cthebdfm = 0.22;

/*!
  * constant of EB-AFM and EB-DFM (0.122*2.5, See F. Dehoux thesis)
  */
double cs_turb_xclt = 0.305;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur,
                             int     **hybrid_turb);

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
cs_f_turb_les_model_get_pointers(int     **idries);

void
cs_f_turb_hybrid_model_get_pointers(int  **iicc,
                                    int  **ishield);

void
cs_f_turb_reference_values(double  **almax,
                           double  **uref);

void
cs_f_turb_model_constants_get_pointers(double  **cmu,
                                       double  **csmago,
                                       double  **xlesfd,
                                       double  **xlesfl,
                                       double  **ales,
                                       double  **bles,
                                       double  **cdries,
                                       double  **csrij,
                                       double  **xclt);

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
 *   hybrid_turb --> pointer to cs_glob_turb_model->hybrid_turb
 *----------------------------------------------------------------------------*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur,
                             int     **hybrid_turb)
{
  *iturb  = &(_turb_model.iturb);
  *itytur = &(_turb_model.itytur);
  *hybrid_turb = &(_turb_model.hybrid_turb);
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
 *----------------------------------------------------------------------------*/

void
cs_f_turb_les_model_get_pointers(int     **idries)
{
  *idries = &(_turb_les_model.idries);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the hybrid turbulence model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iicc    --> pointer to cs_glob_turb_hybrid_model->iicc
 *   ishield --> pointer to cs_glob_turb_hybrid_model->ishield
 *----------------------------------------------------------------------------*/

void
cs_f_turb_hybrid_model_get_pointers(int  **iicc,
                                    int  ** ishield)
{
  *iicc    = &(_turb_hybrid_model.iicc);
  *ishield = &(_turb_hybrid_model.ishield);
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
                           double  **uref)
{
  *almax  = &(_turb_ref_values.almax);
  *uref   = &(_turb_ref_values.uref);
}

/*----------------------------------------------------------------------------
 * Get pointers to constants for turbulence models.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_turb_model_constants_get_pointers(double  **cmu,
                                       double  **csmago,
                                       double  **xlesfd,
                                       double  **xlesfl,
                                       double  **ales,
                                       double  **bles,
                                       double  **cdries,
                                       double  **csrij,
                                       double  **xclt)
{
  *cmu    = &cs_turb_cmu;
  *csmago= &cs_turb_csmago;
  *csmago= &cs_turb_csmago;
  *xlesfd= &cs_turb_xlesfd;
  *xlesfl= &cs_turb_xlesfl;
  *ales  = &cs_turb_ales;
  *bles  = &cs_turb_bles;
  *cdries= &cs_turb_cdries;
  *csrij = &cs_turb_csrij;
  *xclt = &cs_turb_xclt;
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return enumeration name associated with a turbulence model value
 *
 * parameters:
 *   id <-- model type
 *
 * returns:
 *   pointer to enum name.
 *----------------------------------------------------------------------------*/

static const char *
_turbulence_model_enum_name(cs_turb_model_type_t  id)
{
  const char *s = NULL;
  switch(id) {
  case CS_TURB_NONE:
    s = "CS_TURB_NONE";
    break;
  case CS_TURB_MIXING_LENGTH:
    s = "CS_TURB_MIXING_LENGTH";
    break;
  case CS_TURB_K_EPSILON:
    s = "CS_TURB_K_EPSILON";
    break;
  case CS_TURB_K_EPSILON_LIN_PROD:
    s = "CS_TURB_K_EPSILON_LIN_PROD";
    break;
  case CS_TURB_K_EPSILON_LS:
    s = "CS_TURB_K_EPSILON_LS";
    break;
  case CS_TURB_K_EPSILON_QUAD:
    s = "CS_TURB_K_EPSILON_QUAD";
    break;
  case CS_TURB_RIJ_EPSILON_LRR:
    s = "CS_TURB_RIJ_EPSILON_LRR";
    break;
  case CS_TURB_RIJ_EPSILON_SSG:
    s = "CS_TURB_RIJ_EPSILON_SSG";
    break;
  case CS_TURB_RIJ_EPSILON_EBRSM:
    s = "CS_TURB_RIJ_EPSILON_EBRSM";
    break;
  case CS_TURB_LES_SMAGO_CONST:
    s = "CS_TURB_LES_SMAGO_CONST";
    break;
  case CS_TURB_LES_SMAGO_DYN:
    s = "CS_TURB_LES_SMAGO_DYN";
    break;
  case CS_TURB_LES_WALE:
    s = "CS_TURB_LES_WALE";
    break;
  case CS_TURB_V2F_PHI:
    s = "CS_TURB_V2F_PHI";
    break;
  case CS_TURB_V2F_BL_V2K:
    s = "CS_TURB_V2F_BL_V2K";
    break;
  case CS_TURB_K_OMEGA:
    s = "CS_TURB_K_OMEGA";
    break;
  case CS_TURB_SPALART_ALLMARAS:
    s = "CS_TURB_SPALART_ALLMARAS";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Unknown cs_turb_model_type_t value: %d"), id);
  }

  return s;
}

/*----------------------------------------------------------------------------
 * Return associated with a turbulence model value
 *
 * parameters:
 *   id <-- model type
 *
 * returns:
 *   pointer to enum name.
 *----------------------------------------------------------------------------*/

static const char *
_turbulence_model_name(cs_turb_model_type_t  id)
{
  const char *s = NULL;
  switch(id) {
  case CS_TURB_NONE:
    s = _("no turbulence model");
    break;
  case CS_TURB_MIXING_LENGTH:
    s = _("mixing length model");
    break;
  case CS_TURB_K_EPSILON:
    s = _("standard k-epsilon model");
    break;
  case CS_TURB_K_EPSILON_LIN_PROD:
    s = _("k-epsilon model with Linear Production (LP) correction");
    break;
  case CS_TURB_K_EPSILON_LS:
    s = _("Launder-Sharma k-epsilon model");
    break;
  case CS_TURB_K_EPSILON_QUAD:
    s = _("Baglietto et al. quadratic k-epsilon model");
    break;
  case CS_TURB_RIJ_EPSILON_LRR:
    s = _("Rij-epsilon (LRR) model");
    break;
  case CS_TURB_RIJ_EPSILON_SSG:
    s = _("Rij-epsilon (SSG)");
    break;
  case CS_TURB_RIJ_EPSILON_EBRSM:
    s = _("Rij-epsilon (EBRSM))");
    break;
  case CS_TURB_LES_SMAGO_CONST:
    s = _("LES (constant Smagorinsky model)");
    break;
  case CS_TURB_LES_SMAGO_DYN:
    s = _("LES (classical dynamic Smagorisky model)");
    break;
  case CS_TURB_LES_WALE:
    s = _("LES (WALE)");
    break;
  case CS_TURB_V2F_PHI:
    s = _("v2f phi-model");
    break;
  case CS_TURB_V2F_BL_V2K:
    s = _("v2f BL-v2-k)");
    break;
  case CS_TURB_K_OMEGA:
    s = _("k-omega SST");
    break;
  case CS_TURB_SPALART_ALLMARAS:
    s = _("Spalart-Allmaras model");
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Unknown cs_turb_model_type_t value: %d"), id);
  }

  return s;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize turbulence model structures
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_model_init(void) {

  /* set global pointer to turbulence model */
  cs_set_glob_turb_model();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize type and order members of turbulence model structure
 */
/*----------------------------------------------------------------------------*/

void
cs_set_type_order_turbulence_model(void)
{
  _turb_model.type = CS_TURB_NONE;
  if (_turb_model.iturb == CS_TURB_MIXING_LENGTH) {
     _turb_model.type = CS_TURB_RANS;
     _turb_model.order = CS_TURB_ALGEBRAIC;
  }
  else if (   _turb_model.iturb == CS_TURB_K_EPSILON
           || _turb_model.iturb == CS_TURB_K_EPSILON_LIN_PROD
           || _turb_model.iturb == CS_TURB_K_EPSILON_LS
           || _turb_model.iturb == CS_TURB_K_EPSILON_QUAD
           || _turb_model.iturb == CS_TURB_V2F_PHI
           || _turb_model.iturb == CS_TURB_V2F_BL_V2K
           || _turb_model.iturb == CS_TURB_K_OMEGA
           || _turb_model.iturb == CS_TURB_SPALART_ALLMARAS) {
    _turb_model.type = CS_TURB_RANS;
    _turb_model.order = CS_TURB_FIRST_ORDER;
  }
  else if (   _turb_model.iturb == CS_TURB_RIJ_EPSILON_LRR
           || _turb_model.iturb == CS_TURB_RIJ_EPSILON_SSG
           || _turb_model.iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    _turb_model.type = CS_TURB_RANS;
    _turb_model.order = CS_TURB_SECOND_ORDER;
  }
  else if (   _turb_model.iturb == CS_TURB_LES_SMAGO_CONST
           || _turb_model.iturb == CS_TURB_LES_SMAGO_DYN
           || _turb_model.iturb == CS_TURB_LES_WALE) {
    _turb_model.type = CS_TURB_LES;
    _turb_model.order = CS_TURB_ALGEBRAIC;
  }

  else {
    _turb_model.iturb = 0;
    _turb_model.itytur = CS_TURB_TYPE_NONE;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide write access to turbulence model structure
 */
/*----------------------------------------------------------------------------*/

cs_turb_model_t *
cs_get_glob_turb_model(void)
{
  return &_turb_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global pointer to turbulence model structure
 *
 * This global pointer provides a read-only access to the structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_set_glob_turb_model(void)
{
  /* If not set yet, points to the locally defined structure */
  if (cs_glob_turb_model == NULL)
    cs_glob_turb_model = &_turb_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute turbulence model constants,
 *        some of which may depend on the model choice.
 *
 * \param[in]       phase_id  turbulent phase id (-1 for single phase flow)
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_compute_constants(int phase_id)
{
  cs_turb_dpow   = 1./(1.+cs_turb_bpow);

  if (cs_glob_turb_model->itytur == 5)
    cs_turb_cmu = 0.22;

  cs_turb_cmu025 = pow(cs_turb_cmu, 0.25);
  cs_turb_cstlog_alpha = exp(-cs_turb_xkappa
                             * (cs_turb_cstlog_rough - cs_turb_cstlog));

  int k_turb_schmidt = cs_field_key_id("turbulent_schmidt");

  cs_field_pointer_ensure_init();

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_phi = CS_F_(phi);
  cs_field_t *f_eps = CS_F_(eps);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_phi = CS_FI_(phi, phase_id);
    f_eps = CS_FI_(eps, phase_id);
  }

  if (f_k != NULL)
    cs_field_set_key_double(f_k, k_turb_schmidt, 1.);

  if (f_phi != NULL)
    cs_field_set_key_double(f_phi, k_turb_schmidt, 1.);

  if (   cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR
      || cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG)
    cs_field_set_key_double(f_eps, k_turb_schmidt, 1.22);
  else if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    cs_field_set_key_double(f_eps, k_turb_schmidt, 1.15);
    cs_turb_crij3 = 0.6;
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K)
    cs_field_set_key_double(f_eps, k_turb_schmidt, 1.5);
  else
    cs_field_set_key_double(f_eps, k_turb_schmidt, 1.30);

  if (cs_glob_turb_rans_model->idirsm == 0)
    cs_turb_csrij = 0.11;
  else {
    if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
      cs_turb_csrij = 0.21;
    else
      cs_turb_csrij = 0.22;
  }

  if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA){
    /* SST DDES */
    if (cs_glob_turb_model->hybrid_turb == 2)
      cs_turb_cddes = 0.65;
    else if (cs_glob_turb_model->hybrid_turb == 1)
      cs_turb_cddes = 0.61;
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    cs_turb_cddes = 0.60;
  }

  if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA){
    /* SST HTLES */
    cs_turb_chtles_bt0 = 0.48;
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    /* BL-v2/k HTLES */
    cs_turb_chtles_bt0 = 0.42;
  }

  double xkappa2 = cs_turb_xkappa*cs_turb_xkappa;
  cs_turb_ckwgm1 =   cs_turb_ckwbt1/cs_turb_cmu
                   - xkappa2/(cs_turb_ckwsw1*sqrt(cs_turb_cmu));
  cs_turb_ckwgm2 =   cs_turb_ckwbt2/cs_turb_cmu
                   - xkappa2/(cs_turb_ckwsw2*sqrt(cs_turb_cmu));
  cs_turb_csaw1 =   cs_turb_csab1/xkappa2
                  + 1./cs_turb_csasig*(1. + cs_turb_csab2);
  cs_turb_csmago_max = cs_turb_csmago*cs_turb_csmago;
  cs_turb_csmago_min = 0.;

  /* In case of Rotta model (ie LRR + Cr2 = 0) compute
   * automatically the C0 constant */
  if ((cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) &&
      (CS_ABS(cs_turb_crij2) < 1.e-12))
    cs_turb_crij_c0 = (cs_turb_crij1-1.0)*2.0/3.0;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_turb_ref_values
 *
 * needed to initialize structure with GUI
 */
/*----------------------------------------------------------------------------*/

cs_turb_ref_values_t *
cs_get_glob_turb_ref_values(void)
{
  return &_turb_ref_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_turb_rans_model
 *
 * needed to initialize structure with GUI
 */
/*----------------------------------------------------------------------------*/

cs_turb_rans_model_t *
cs_get_glob_turb_rans_model(void)
{
  return &_turb_rans_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_turb_les_model
 *
 * needed to initialize structure with GUI
 */
/*----------------------------------------------------------------------------*/

cs_turb_les_model_t *
cs_get_glob_turb_les_model(void)
{
  return &_turb_les_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_turb_hybrid_model
 *
 * needed to initialize structure with GUI
 */
/*----------------------------------------------------------------------------*/

cs_turb_hybrid_model_t *
cs_get_glob_turb_hybrid_model(void)
{
  return &_turb_hybrid_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the turbulence model parameters to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_model_log_setup(void)
{
  if (cs_glob_turb_model == NULL)
    return;

  const cs_turb_model_t *turb_model = cs_glob_turb_model;
  const cs_wall_functions_t *wall_fns = cs_get_glob_wall_functions();

  cs_var_cal_opt_t var_cal_opt;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  cs_log_printf(CS_LOG_SETUP, _("\n"
                                "Turbulence model options\n"
                                "------------------------\n\n"
                                "  Continuous phase:\n\n"));

  if (turb_model->type == CS_TURB_RANS)
    cs_log_printf(CS_LOG_SETUP,
                  _("    RANS model (type = CS_TURB_RANS)\n"));
  else if (turb_model->type == CS_TURB_LES)
    cs_log_printf(CS_LOG_SETUP,
                  _("   LES model (type = CS_TURB_LES)\n"));
  else if (turb_model->order == CS_TURB_ALGEBRAIC)
    cs_log_printf(CS_LOG_SETUP,
                  _("   Algebraic model (order = CS_TURB_ALGEBRAIC)\n"));
  else if (turb_model->order == CS_TURB_FIRST_ORDER)
    cs_log_printf(CS_LOG_SETUP,
                  _("   First order model (order = CS_TURB_FIRST_ORDER)\n"));
  else if (turb_model->order == CS_TURB_SECOND_ORDER)
    cs_log_printf(CS_LOG_SETUP,
                  _("   Second order model (order = CS_TURB_SECOND_ORDER)\n"));

  cs_log_printf(CS_LOG_SETUP,
                _("\n    %s\n"
                  "      (iturb = %s)\n\n"),
                _turbulence_model_name(turb_model->iturb),
                _turbulence_model_enum_name(turb_model->iturb));

  const char *iwallf_value_str[]
    = {N_("Disabled"),
       N_("One scale power law, forbidden for k-epsilon"),
       N_("One scale log law"),
       N_("Two scales log law"),
       N_("Scalable wall function"),
       N_("Two scales Van Driest"),
       N_("Two scales smooth/rough"),
       N_("All y+")};

  if (wall_fns->iwallf >= 0)
    cs_log_printf(CS_LOG_SETUP,
                  _("    iwallf                      (wall function:\n"
                    "                                 %s)\n"),
                  iwallf_value_str[wall_fns->iwallf]);

  const char *iwalfs_value_str[]
    = {N_("Arparci and Larsen"),
       N_("Van Driest"),
       N_("Louis (atmo flows)"),
       N_("Monin Obukhov (atmo flows)"),
       N_("smooth/rough")};

  if (wall_fns->iwalfs >= 0)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    iwalfs                      (Scalar wall function:\n"
         "                                 %s)\n"),
       iwalfs_value_str[wall_fns->iwalfs]);

  cs_log_printf(CS_LOG_SETUP,
                _("    ypluli:      %14.5e (Limit Y+)\n\n"),
                wall_fns->ypluli);

  const char *igrhok_value_str[]
    = {N_("0 (ignore Grad(rho k) in velocity equation)"),
       N_("1 (Grad(rho k) in velocity equation)")};

  cs_log_printf(CS_LOG_SETUP,
                _("    igrhok:        %s\n"),
                _(igrhok_value_str[cs_glob_turb_rans_model->igrhok]));

  if (turb_model->iturb == CS_TURB_MIXING_LENGTH) {

    cs_log_printf(CS_LOG_SETUP,
                  _("    xlomlg:      %14.5e (Characteristic length)\n"),
                  cs_glob_turb_rans_model->xlomlg);

  }
  else if (   turb_model->iturb == CS_TURB_K_EPSILON
           || turb_model->iturb == CS_TURB_K_EPSILON_LIN_PROD
           || turb_model->iturb == CS_TURB_K_EPSILON_LS
           || turb_model->iturb == CS_TURB_K_EPSILON_QUAD) {

    cs_log_printf
      (CS_LOG_SETUP,
       _("    uref:        %14.5e (Characteristic velocity)\n"
         "    iclkep:      %14d (k-epsilon clipping model)\n"
         "    ikecou:      %14d (k-epsilon coupling mode)\n"
         "    igrake:      %14d (Account for gravity)\n"),
       cs_glob_turb_ref_values->uref,
       cs_glob_turb_rans_model->iclkep,
       cs_glob_turb_rans_model->ikecou,
       cs_glob_turb_rans_model->igrake);

    if (   cs_glob_turb_rans_model->ikecou == 0
        && cs_glob_time_step_options->idtvar >= 0) {
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxve = var_cal_opt.relaxv;
      cs_log_printf(CS_LOG_SETUP,
                    _("    relaxv:      %14.5e for k (Relaxation)\n"
                      "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
                    relaxvk, relaxve);
    }
    else
      cs_log_printf(CS_LOG_SETUP, _("\n"));

  }
  else if (   turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR
           || turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG
           || turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

    cs_log_printf(CS_LOG_SETUP,
                  _("    uref:        %14.5e (Characteristic velocity)\n"
                    "    reinit_turb: %14d (Advanced re-init)\n"
                    "    irijco:      %14d (Coupled resolution)\n"
                    "    irijnu:      %14d (Matrix stabilization)\n"
                    "    irijrb:      %14d (Reconstruct at boundaries)\n"
                    "    igrari:      %14d (Account for gravity)\n"
                    "    iclsyr:      %14d (Symmetry implicitation)\n"
                    "    iclptr:      %14d (Wall implicitation)\n"
                    "    ikwcln:      %14d (Wall boundary condition"
                                           "on omega in k-omega SST)\n"),
                  cs_glob_turb_ref_values->uref,
                  cs_glob_turb_rans_model->reinit_turb,
                  cs_glob_turb_rans_model->irijco,
                  cs_glob_turb_rans_model->irijnu,
                  cs_glob_turb_rans_model->irijrb,
                  cs_glob_turb_rans_model->igrari,
                  cs_glob_turb_rans_model->iclsyr,
                  cs_glob_turb_rans_model->iclptr,
                  cs_glob_turb_rans_model->ikwcln);

    int idirsm = cs_glob_turb_rans_model->idirsm;
    if (idirsm < 0 || idirsm > 1)
      idirsm = 2;

    const char *s_turb_diff_model[3]
      = {"scalar diffusivity: Shir model",
         "tensorial diffusivity: Daly-Harlow model",
         "<unknown>"};

    cs_log_printf(CS_LOG_SETUP,
                  _("    idirsm:      %14d (%s)\n"),
                  idirsm, s_turb_diff_model[idirsm]);

    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
      cs_log_printf(CS_LOG_SETUP,
                    _("    irijec:      %14d (Wall echo terms)\n"
                      "    idifre:      %14d (Handle diffusion tensor)\n"),
                    cs_glob_turb_rans_model->irijec,
                    cs_glob_turb_rans_model->idifre);
    }
    else if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_log_printf(CS_LOG_SETUP,
                    _("    reinit_turb: %14d (turbulence reinitialization)\n"),
                    cs_glob_turb_rans_model->reinit_turb);
    }

  }
  else if (turb_model->type == CS_TURB_LES) {
    cs_log_printf(CS_LOG_SETUP,
                  _("    csmago:      %14.5e (Smagorinsky constant)\n"
                    "    cwale:       %14.5e (WALE model constant)\n"
                    "    xlesfl:      %14.5e (Filter with in a cell is)\n"
                    "    ales:        %14.5e (written as)\n"
                    "    bles:        %14.5e (xlesfl*(ales*volume)**(bles))\n"
                    "    idries:      %14d (=1 Van Driest damping)\n"
                    "    cdries:      %14.5e (Van Driest constant)\n"
                    "    xlesfd:      %14.5e (Ratio between the explicit)\n"
                    "                                (filter and LES filter)\n"
                    "                                (recommended value: 1.5)\n"
                    "    smagmx:      %14.5e (Max Smagorinsky in the)\n"
                    "                                (dynamic model case)\n"),
                  cs_turb_csmago, cs_turb_cwale, cs_turb_xlesfl,
                  cs_turb_ales, cs_turb_bles, cs_glob_turb_les_model->idries,
                  cs_turb_cdries, cs_turb_xlesfd, cs_turb_csmago_max);

  }
  else if (turb_model->iturb == CS_TURB_V2F_PHI) {

    cs_log_printf(CS_LOG_SETUP,
                  _("    uref:        %14.5e (Characteristic velocity)\n"
                    "    iclkep:      %14d (k-epsilon clipping model)\n"
                    "    ikecou:      %14d (k-epsilon coupling mode)\n"
                    "    igrake:      %14d (Account for gravity)\n"),
                  cs_glob_turb_ref_values->uref,
                  cs_glob_turb_rans_model->iclkep,
                  cs_glob_turb_rans_model->ikecou,
                  cs_glob_turb_rans_model->igrake);

    if (   cs_glob_turb_rans_model->ikecou == 0
        && cs_glob_time_step_options->idtvar >= 0) {

      cs_real_t relaxvk, relaxve;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      relaxve = var_cal_opt.relaxv;
      cs_log_printf(CS_LOG_SETUP,
                    _("    relaxv:      %14.5e for k (Relaxation)\n"
                      "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
                    relaxvk, relaxve);

    }
    else
      cs_log_printf(CS_LOG_SETUP,_("\n"));

  }
  else if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {

    const char *hybrid_turb_value_str[]
      = {N_("CS_HYBRID_NONE (no RANS-LES hybrid model)"),
         N_("CS_HYBRID_DES (RANS-LES hybrid model)"),
         N_("CS_HYBRID_DDES  (RANS-LES hybrid model)"),
         N_("CS_HYBRID_SAS (Scale Adpative Model)"),
         N_("CS_HYBRID_HTLES (Hybrid Temporal LES)")};

    cs_log_printf(CS_LOG_SETUP,
                  _("    uref:        %14.5e (Characteristic velocity)\n"
                    "    iclkep:      %14d (k-epsilon clipping model)\n"
                    "    ikecou:      %14d (k-epsilon coupling mode)\n"
                    "    hybrid_turb: %s\n"
                    "    igrake:      %14d (Account for gravity)\n"),
                  cs_glob_turb_ref_values->uref,
                  cs_glob_turb_rans_model->iclkep,
                  cs_glob_turb_rans_model->ikecou,
                  hybrid_turb_value_str[turb_model->hybrid_turb],
                  cs_glob_turb_rans_model->igrake);

    if (   cs_glob_turb_rans_model->ikecou == 0
        && cs_glob_time_step_options->idtvar >= 0) {

      cs_real_t relaxvk, relaxve;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      relaxve = var_cal_opt.relaxv;
      cs_log_printf(CS_LOG_SETUP,
                    _("    relaxv:      %14.5e for k (Relaxation)\n"
                      "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
                    relaxvk, relaxve);

    }
    else
      cs_log_printf(CS_LOG_SETUP,_("\n"));

  }
  else if (turb_model->iturb == CS_TURB_K_OMEGA) {

    const char *hybrid_turb_value_str[]
      = {N_("CS_HYBRID_NONE (no RANS-LES hybrid model)"),
         N_("CS_HYBRID_DES (RANS-LES hybrid model)"),
         N_("CS_HYBRID_DDES  (RANS-LES hybrid model)"),
         N_("CS_HYBRID_SAS (Scale Adpative Model)"),
         N_("CS_HYBRID_HTLES (Hybrid Temporal LES)")};

    cs_log_printf(CS_LOG_SETUP,
                  _("    uref:        %14.5e (Characteristic velocity)\n"
                    "    ikecou:      %14d (k-epsilon coupling mode)\n"
                    "    reinit_turb: %14d (Advanced re-init)\n"
                    "    hybrid_turb: %s\n"
                    "    igrake:      %14d (Account for gravity)\n"),
                  cs_glob_turb_ref_values->uref,
                  cs_glob_turb_rans_model->ikecou,
                  cs_glob_turb_rans_model->reinit_turb,
                  hybrid_turb_value_str[turb_model->hybrid_turb],
                  cs_glob_turb_rans_model->igrake);

    if (   cs_glob_turb_rans_model->ikecou == 0
        && cs_glob_time_step_options->idtvar >= 0) {

      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(omg), key_cal_opt_id, &var_cal_opt);
      cs_real_t relaxvo = var_cal_opt.relaxv;
      cs_log_printf(CS_LOG_SETUP,
                    _("    relaxv:      %14.5e for k (Relaxation)\n"
                      "    relaxv:      %14.5e for omega (Relaxation)\n"),
                    relaxvk, relaxvo);
    }
    else
      cs_log_printf(CS_LOG_SETUP,_("\n"));

  }
  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {

    cs_field_get_key_struct(CS_F_(nusa), key_cal_opt_id, &var_cal_opt);
    cs_log_printf(CS_LOG_SETUP,
                  _("    uref:        %14.5e (Characteristic velocity)\n"
                    "    relaxv:      %14.5e for nu (Relaxation)\n"),
                  cs_glob_turb_ref_values->uref,
                  var_cal_opt.relaxv);

  }

  if (   turb_model->type == CS_TURB_RANS
      && turb_model->order == CS_TURB_FIRST_ORDER){

    const char *irccor_value_str[]
      = {N_("0 (no rotation/curvature correction)"),
         N_("1 (rotation/curvature correction)")};
    cs_log_printf(CS_LOG_SETUP,
                  _("    irccor:        %s\n"),
                  irccor_value_str[cs_glob_turb_rans_model->irccor]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the turbulent constants to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_constants_log_setup(void)
{
  if (cs_glob_turb_model == NULL)
    return;

  const cs_turb_model_t *turb_model = cs_glob_turb_model;

  cs_log_printf(CS_LOG_SETUP,
     _("\n"
       "  Constants:\n\n"
       "    xkappa:      %14.5e (Von Karman constant)\n"
       "    cstlog:      %14.5e (U+=Log(y+)/kappa +cstlog)\n"
       "    apow:        %14.5e (U+=apow (y+)**bpow (W&W law))\n"
       "    bpow:        %14.5e (U+=apow (y+)**bpow (W&W law))\n\n"),
       cs_turb_xkappa, cs_turb_cstlog, cs_turb_apow, cs_turb_bpow);

  if (turb_model->iturb != CS_TURB_NONE)
    cs_log_printf(CS_LOG_SETUP,
                  _("  %s constants:\n"),
                  _turbulence_model_name(turb_model->iturb));

  if (   turb_model->iturb == CS_TURB_K_EPSILON
      || turb_model->iturb == CS_TURB_K_EPSILON_LIN_PROD
      || turb_model->iturb == CS_TURB_K_EPSILON_LS
      || turb_model->iturb == CS_TURB_K_EPSILON_QUAD)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    ce1:         %14.5e (Cepsilon 1: production coef.)\n"
         "    ce2:         %14.5e (Cepsilon 2: dissipat.  coef.)\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_ce1, cs_turb_ce2, cs_turb_cmu);

  else if (turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    ce1:         %14.5e (Cepsilon 1: production coef.)\n"
         "    ce2:         %14.5e (Cepsilon 2: dissipat.  coef.)\n"
         "    crij1:       %14.5e (Slow term coefficient)\n"
         "    crij2:       %14.5e (Fast term coefficient)\n"
         "    crij3:       %14.5e (Gravity term coefficient)\n"
         "    csrij:       %14.5e (Rij diffusion coeff.)\n"
         "    crijp1:      %14.5e (Slow coeff. for wall echo)\n"
         "    crijp2:      %14.5e (Fast coeff. for wall echo)\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_ce1, cs_turb_ce2, cs_turb_crij1, cs_turb_crij2,
         cs_turb_crij3, cs_turb_csrij, cs_turb_crijp1,
         cs_turb_crijp2, cs_turb_cmu);

  else if (turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    cssgs1:      %14.5e (Cs1 coeff.)\n"
         "    cssgs2:      %14.5e (Cs2 coeff.)\n"
         "    cssgr1:      %14.5e (Cr1 coeff.)\n"
         "    cssgr2:      %14.5e (Cr2 coeff.)\n"
         "    cssgr3:      %14.5e (Cr3 coeff.)\n"
         "    cssgr4:      %14.5e (Cr4 coeff.)\n"
         "    cssgr5:      %14.5e (Cr5 coeff.)\n"
         "    csrij:       %14.5e (Rij Cs diffusion coeff.)\n"
         "    crij3:       %14.5e (Gravity term coeff.)\n"
         "    ce1:         %14.5e (Ceps1 coeff.)\n"
         "    cssge2:      %14.5e (Ceps2 coeff.)\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_cssgs1, cs_turb_cssgs2, cs_turb_cssgr1,
         cs_turb_cssgr2, cs_turb_cssgr3, cs_turb_cssgr4,
         cs_turb_cssgr5, cs_turb_csrij, cs_turb_crij3,
         cs_turb_ce1, cs_turb_cssge2,
         cs_turb_cmu);

  else if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    cs_turb_crij3 = 0.6;
    cs_log_printf
      (CS_LOG_SETUP,
       _("    cebms1:      %14.5e (Cs1 coeff.)\n"
         "    cebmr1:      %14.5e (Cr1 coeff.)\n"
         "    cebmr2:      %14.5e (Cr2 coeff.)\n"
         "    cebmr3:      %14.5e (Cr3 coeff.)\n"
         "    cebmr4:      %14.5e (Cr4 coeff.)\n"
         "    cebmr5:      %14.5e (Cr5 coeff.)\n"
         "    csrij:       %14.5e (Rij Cs diffusion coeff.)\n"
         "    crij3:       %14.5e (Gravity term coeff.)\n"
         "    cebme2:      %14.5e (Coef Ceps2)\n"
         "    ce1:         %14.5e (Coef Ceps1)\n"
         "    xa1:         %14.5e (Coef A1)\n"
         "    xceta:       %14.5e (Coef Ceta)\n"
         "    xct:         %14.5e (Coef CT)\n"),
         cs_turb_cebms1, cs_turb_cebmr1, cs_turb_cebmr2,
         cs_turb_cebmr3, cs_turb_cebmr4, cs_turb_cebmr5,
         cs_turb_csrij, cs_turb_crij3,  cs_turb_cebme2,
         cs_turb_ce1, cs_turb_xa1,
         cs_turb_xceta, cs_turb_xct);

  }

  else if (turb_model->iturb == CS_TURB_V2F_PHI)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    cv2fa1:      %14.5e (a1 to calculate Cepsilon1)\n"
         "    cv2fe2:      %14.5e (Cepsilon 2: dissip. coeff.)\n"
         "    cmu   :      %14.5e (Cmu constant)\n"
         "    cv2fct:      %14.5e (CT constant)\n"
         "    cv2fcl:      %14.5e (CL constant)\n"
         "    cv2fet:      %14.5e (C_eta constant)\n"
         "    cv2fc1:      %14.5e (C1 constant)\n"
         "    cv2fc2:      %14.5e (C2 constant)\n"),
         cs_turb_cv2fa1, cs_turb_cv2fe2,
         cs_turb_cmu, cs_turb_cv2fct,
         cs_turb_cv2fcl, cs_turb_cv2fet, cs_turb_cv2fc1,
         cs_turb_cv2fc2);

  else if (turb_model->iturb == CS_TURB_V2F_BL_V2K)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    cpale1:      %14.5e (Cepsilon 1 : Prod. coeff.)\n"
         "    cpale2:      %14.5e (Cepsilon 2 : Diss. coeff.)\n"
         "    cpale3:      %14.5e (Cepsilon 3 : E term coeff.)\n"
         "    cpale4:      %14.5e (Cepsilon 4 : Mod Diss. coef.)\n"
         "    cmu   :      %14.5e (Cmu constant)\n"
         "    cpalct:      %14.5e (CT constant)\n"
         "    cpalcl:      %14.5e (CL constant)\n"
         "    cpalet:      %14.5e (C_eta constant)\n"
         "    cpalc1:      %14.5e (C1 constant)\n"
         "    cpalc2:      %14.5e (C2 constant)\n"),
         cs_turb_cpale1, cs_turb_cpale2, cs_turb_cpale3,
         cs_turb_cpale4,
         cs_turb_cmu, cs_turb_cpalct, cs_turb_cpalcl,
         cs_turb_cpalet, cs_turb_cpalc1, cs_turb_cpalc2);

  else if (turb_model->iturb == CS_TURB_K_OMEGA)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    ckwsk1:      %14.5e (sigma_k1 constant)\n"
         "    ckwsk2:      %14.5e (sigma_k2 constant)\n"
         "    ckwsw1:      %14.5e (sigma_omega1 constant)\n"
         "    ckwsw2:      %14.5e (sigma_omega2 constant)\n"
         "    ckwbt1:      %14.5e (beta1 constant)\n"
         "    ckwbt2:      %14.5e (beta2 constant)\n"
         "    ckwgm1:      %14.5e (gamma1 constant)\n"
         "    ckwgm2:      %14.5e (gamma2 constant)\n"
         "    ckwa1:       %14.5e (a1 constant to compute mu_t)\n"
         "    ckwc1:       %14.5e (c1 const. for prod. limiter)\n"
         "    cmu:         %14.5e (Cmu (or Beta*) constant for)\n"
         "                          omega/epsilon conversion)\n"),
       cs_turb_ckwsk1, cs_turb_ckwsk2, cs_turb_ckwsw1,
       cs_turb_ckwsw2, cs_turb_ckwbt1, cs_turb_ckwbt2,
       cs_turb_ckwgm1, cs_turb_ckwgm2, cs_turb_ckwa1,
       cs_turb_ckwc1, cs_turb_cmu);

  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS)
    cs_log_printf(CS_LOG_SETUP,
                  _("    csab1:        %14.5e (b1 constant)\n"
                    "    csab2:        %14.5e (b2 constant)\n"
                    "    csasig:       %14.5e (sigma constant)\n"
                    "    csav1:        %14.5e (v1 constant)\n"
                    "    csaw1:        %14.5e (w1 constant)\n"
                    "    csaw2:        %14.5e (w2 constant)\n"
                    "    csaw3:        %14.5e (w3 constant)\n"),
                  cs_turb_csab1, cs_turb_csab2, cs_turb_csasig,
                  cs_turb_csav1, cs_turb_csaw1, cs_turb_csaw2,
                  cs_turb_csaw3);

  int iokss = 0, iokcaz = 0;

  if (cs_glob_turb_rans_model != NULL) {
    if (cs_glob_turb_rans_model->irccor == 1) {
      if (cs_glob_turb_rans_model->itycor == 1)
        iokcaz = 1;
      else if (cs_glob_turb_rans_model->itycor == 2)
        iokss = 1;
    }
  }

  if (iokcaz > 0)
    cs_log_printf(CS_LOG_SETUP,
                  _("   Rotation/curvature correction (Cazalbou)\n"
                    "    ccaze2:       %14.5e (Coef Ce2^0)\n"
                    "    ccazsc:       %14.5e (Coef Csc)\n"
                    "    ccaza:        %14.5e (Coef a)\n"
                    "    ccazb:        %14.5e (Coef b)\n"
                    "    ccazc:        %14.5e (Coef c)\n"
                    "    ccazd:        %14.5e (Coef d)\n"),
                  cs_turb_ccaze2, cs_turb_ccazsc, cs_turb_ccaza,
                  cs_turb_ccazb, cs_turb_ccazc, cs_turb_ccazd);

  if (iokss > 0)
    cs_log_printf(CS_LOG_SETUP,
                  _("   Rotation/curvature correction (Spalart-Shur)\n"
                    "    cssr1:       %14.5e (Coef c_r1)\n"
                    "    cssr2:       %14.5e (Coef c_r2)\n"
                    "    cssr3:       %14.5e (Coef c_r3)\n"),
                  cs_turb_cssr1, cs_turb_cssr2, cs_turb_cssr3);

  const int kturt = cs_field_key_id("turbulent_flux_model");
  const int n_fields = cs_field_n_fields();
  bool has_dfm = false;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE))
      continue;

    const int turb_flux_model = cs_field_get_key_int(f, kturt);
    const int turb_flux_model_type = turb_flux_model / 10;
    has_dfm = has_dfm || (turb_flux_model_type == 3);
  }

  if (has_dfm)
    cs_log_printf(CS_LOG_SETUP,
                  _("  Differential Flux Model (DFM)\n"
                    "    c1trit:      %14.5e (Coef c_1)\n"
                    "    c2trit:      %14.5e (Coef c_2)\n"
                    "    c3trit:      %14.5e (Coef c_3)\n"
                    "    c4trit:      %14.5e (Coef c_4)\n"),
                  cs_turb_c1trit, cs_turb_c2trit, cs_turb_c3trit,
                  cs_turb_c4trit);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute characteristic length for turbulence if not already done.
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_init_ref_quantities(void)
{
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  if (_turb_ref_values.almax < 0)
    _turb_ref_values.almax = pow(mq->tot_vol, 1./3.);

  if (   _turb_model.type == CS_TURB_RANS
      || _turb_model.type == CS_TURB_HYBRID) {
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n"
         " Characteristic values for turbulence:\n"
         " -------------------------------------\n\n"
         "   almax: %11.3g (caracteristic length for initialization)\n"),
       _turb_ref_values.almax);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping for the turbulence flux vector
 *
 * \param[in]       flux_id       turbulent flux index
 * \param[in]       variance_id   scalar variance index
 *
 *----------------------------------------------------------------------------*/

void
cs_clip_turbulent_fluxes(int  flux_id,
                         int  variance_id)
{
  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const cs_real_6_t *cvar_rij = (const cs_real_6_t *)CS_F_(rij)->val;
  const cs_real_t *cvar_tt
    = (const cs_real_t *)cs_field_by_id(variance_id)->val;
  cs_real_3_t *cvar_rit = (cs_real_3_t *) cs_field_by_id(flux_id)->val;
  cs_real_3_t *cvar_clip_rit = NULL;

  cs_field_t *field_rit = cs_field_by_id(flux_id);

  /* Local variables */
  const cs_real_t tol_jacobi = 1.0e-12;
  const cs_real_t l_threshold = 1.0e12;
  cs_lnum_t iclip = 0;
  cs_real_t flsq, maj;
  cs_real_33_t rij;
  cs_real_33_t eigvect = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  cs_real_3_t eigval;
  cs_real_3_t rot_rit;
  cs_real_3_t rit;

  /* Get clippings field for DFM */
  cs_lnum_t kclipp, clip_rit_id;
  kclipp = cs_field_key_id("clipping_id");
  clip_rit_id = cs_field_get_key_int(field_rit, kclipp);
  if (clip_rit_id >= 0) {
    cvar_clip_rit = (cs_real_3_t *) cs_field_by_id(clip_rit_id)->val;
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cvar_clip_rit[cell_id][0] = 0;
      cvar_clip_rit[cell_id][1] = 0;
      cvar_clip_rit[cell_id][2] = 0;
    }
  }

  cs_real_t rit_min_prcoord[3] = {l_threshold,l_threshold,l_threshold};
  cs_real_t rit_max_prcoord[3] = {-l_threshold,-l_threshold,-l_threshold};

  cs_lnum_t iclip_tab_max[3] = {0,0,0};
  cs_lnum_t iclip_tab_min[3] = {0,0,0};

  for (int cell_id = 0; cell_id < n_cells; cell_id++) {
    rij[0][0] = cvar_rij[cell_id][0];
    rij[1][1] = cvar_rij[cell_id][1];
    rij[2][2] = cvar_rij[cell_id][2];
    rij[0][1] = cvar_rij[cell_id][3];
    rij[1][0] = cvar_rij[cell_id][3];
    rij[0][2] = cvar_rij[cell_id][5];
    rij[2][0] = cvar_rij[cell_id][5];
    rij[1][2] = cvar_rij[cell_id][4];
    rij[2][1] = cvar_rij[cell_id][4];

    rit[0] = cvar_rit[cell_id][0];
    rit[1] = cvar_rit[cell_id][1];
    rit[2] = cvar_rit[cell_id][2];

    cs_math_33_eig_val_vec(rij,tol_jacobi,eigval,eigvect);
    cs_math_33_3_product(eigvect,cvar_rit[cell_id],rot_rit);

    for (int i = 0; i < 3; i++) {
      rit_min_prcoord[i] = CS_MIN(rit_min_prcoord[i],rot_rit[i]);
      rit_max_prcoord[i] = CS_MAX(rit_max_prcoord[i],rot_rit[i]);
    }

    for (int i = 0; i < 3; i++) {
      flsq = pow(rot_rit[i],2);
      maj = eigval[i]*cvar_tt[cell_id];
      if ((flsq > 1.0e-12) && (flsq > maj)) {
        rot_rit[i] = rot_rit[i]*sqrt(maj/flsq);
        iclip = 1;
        if (rot_rit[i] > 0)
          iclip_tab_max[i] += 1;
        else
          iclip_tab_min[i] += 1;
      }
    }

    if (iclip > 0) {
      cs_math_33t_3_product(eigvect,rot_rit,cvar_rit[cell_id]);
      cvar_clip_rit[cell_id][0] = cvar_rit[cell_id][0] - rit[0];
      cvar_clip_rit[cell_id][1] = cvar_rit[cell_id][1] - rit[1];
      cvar_clip_rit[cell_id][2] = cvar_rit[cell_id][2] - rit[2];
    }
  }

  cs_lnum_t iclip_max = iclip_tab_max[0] + iclip_tab_max[1]
                      + iclip_tab_max[2];
  cs_lnum_t iclip_min = iclip_tab_min[0] + iclip_tab_min[1]
                      + iclip_tab_min[2];

  /* Save clippings for log */
  cs_log_iteration_clipping_field(flux_id,
                                  iclip_min,
                                  iclip_max,
                                  rit_min_prcoord,
                                  rit_max_prcoord,
                                  iclip_tab_min,
                                  iclip_tab_max);
}

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
                         void             *vals)
{
  CS_UNUSED(input);

  cs_real_t *k = vals;

  cs_assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_turb_model_t *tm = cs_glob_turb_model;

  const cs_real_t *val_k = (CS_F_(k) != NULL) ? CS_F_(k)->val : NULL;

  if (val_k != NULL) {
    if (elt_ids != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = elt_ids[i];
        k[i] = val_k[c_id];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        k[i] = val_k[i];
      }
    }
  }
  else if (tm->itytur == 3) {
    const cs_real_6_t *rij = (const cs_real_6_t *)CS_F_(rij)->val;
    if (elt_ids != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = elt_ids[i];
        k[i] = (rij[c_id][0] + rij[c_id][1] + rij[c_id][2]) * 0.5;
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        k[i] = (rij[i][0] + rij[i][1] + rij[i][2]) * 0.5;
      }
    }
  }
  else {
    if (cs_log_default_is_active()) {
      cs_log_warning(_("%s: cannot simply determine k from other variables\n"
                       "with turbulence model %s.\n"),
                     __func__,
                     _turbulence_model_enum_name(tm->iturb));
    }

    for (cs_lnum_t i = 0; i < n_elts; i++)
      k[i] = 0;
  }
}

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
                           void             *vals)
{
  CS_UNUSED(input);

  cs_real_t *eps = vals;

  cs_assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_turb_model_t *tm = cs_glob_turb_model;

  const cs_real_t *val_eps = (CS_F_(eps) != NULL) ? CS_F_(eps)->val : NULL;

  if (val_eps != NULL) {
    if (elt_ids != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = elt_ids[i];
        eps[i] = val_eps[c_id];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        eps[i] = val_eps[i];
      }
    }
  }
  else if (tm->iturb == CS_TURB_K_OMEGA) {
    const cs_real_t *val_k = CS_F_(omg)->val;
    const cs_real_t *val_omg = CS_F_(omg)->val;
    if (elt_ids != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = elt_ids[i];
        eps[i] = cs_turb_cmu * val_k[c_id] * val_omg[c_id];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        eps[i] = cs_turb_cmu * val_k[i] * val_omg[i];
      }
    }
  }
  else {
    if (cs_log_default_is_active()) {
      cs_log_warning(_("%s: cannot simply determine k from other variables\n"
                       "with turbulence model %s.\n"),
                     __func__,
                     _turbulence_model_enum_name(tm->iturb));
    }

    for (cs_lnum_t i = 0; i < n_elts; i++)
      eps[i] = 0;
  }
}

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
                           void             *vals)
{
  CS_UNUSED(input);

  cs_real_6_t *rij = vals;

  cs_assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_turb_model_t *tm = cs_glob_turb_model;

  const cs_real_6_t *val_rij = NULL;
  if (CS_F_(rij) != NULL)
    val_rij = (const cs_real_6_t *)CS_F_(rij)->val;

  /* Rij already present */

  if (val_rij != NULL) {
    if (elt_ids != NULL) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t c_id = elt_ids[i];
        for (cs_lnum_t j = 0; j < 6; j++)
          rij[i][j] = val_rij[c_id][j];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        for (cs_lnum_t j = 0; j < 6; j++)
          rij[i][j] = val_rij[i][j];
      }
    }
  }

  /* Rij estimated from turbulent viscosity and velocity gradient */

  else if (CS_F_(k) != NULL) {

    const cs_real_t d2o3 = 2./3.;
    const cs_real_t *cpro_mu_t = CS_F_(mu_t)->val;
    const cs_real_t *cpro_rho = CS_F_(rho)->val;
    const cs_real_t *cvar_k = CS_F_(k)->val;

    const cs_mesh_t *m = cs_glob_mesh;
    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
    cs_halo_type_t halo_type
      = (m->cell_cells_idx != NULL) ? CS_HALO_EXTENDED : CS_HALO_STANDARD;

    const cs_field_t *f_vel = CS_F_(vel);

#   pragma omp parallel for if(n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {

      cs_lnum_t c_id = (elt_ids != NULL) ? elt_ids[i] : i;

      cs_real_t gradv[6][3];

      cs_gradient_vector_cell(m,
                              fvq,
                              c_id,
                              halo_type,
                              (const cs_real_3_t *)f_vel->bc_coeffs->a,
                              (const cs_real_33_t *)f_vel->bc_coeffs->b,
                              (const cs_real_3_t *)f_vel->val,
                              NULL,
                              gradv);

      cs_real_t divu = gradv[0][0] + gradv[1][1] + gradv[2][2];
      cs_real_t nut = cpro_mu_t[c_id] / cpro_rho[c_id];
      cs_real_t xdiag = d2o3*(cvar_k[c_id]+ nut*divu);

      rij[i][0] =  xdiag - 2.*nut*gradv[0][0];
      rij[i][1] =  xdiag - 2.*nut*gradv[1][1];
      rij[i][2] =  xdiag - 2.*nut*gradv[2][2];
      rij[i][3] = -nut*(gradv[1][0]+gradv[0][1]);
      rij[i][4] = -nut*(gradv[2][1]+gradv[1][2]);
      rij[i][5] = -nut*(gradv[2][0]+gradv[0][2]);

    }

  }

  else {
    if (cs_log_default_is_active()) {
      cs_log_warning(_("%s: cannot simply determine Rij from other variables\n"
                       "with turbulence model %s.\n"),
                     __func__,
                     _turbulence_model_enum_name(tm->iturb));
    }

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t j = 0; j < 6; j++)
        rij[i][j] = 0;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*============================================================================
 * Base turbulence model data.
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

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_parameters.h"
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
        - 0: scalar diffusivity (Shir model)
        - 1: tensorial diffusivity (Daly and Harlow model, default model)
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
        - 1: true
        - 0: false (default)
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
        accurate treatment of \f$ \tens{R} \f$ at the boundary (see \ref condli)
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
        studied with Code_Saturne, they did not bring any improvement in the results.\n
        In addition, their use induces an increase in the calculation time.\n
        The wall echo terms imply the calculation of the distance to the wall
        for every cell in the domain. See \ref optcal::icdpar "icdpar" for potential
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
  \var  cs_turb_les_model_t::ivrtex
        Activates or not the generation of synthetic turbulence at the
        different inlet boundaries with the LES model (generation of
        unsteady synthetic eddies).\n
        - 1: true
        - 0: false (default)
        Useful if \ref iturb =40, 41 or 42\n
        This keyword requires the completion of the routine  \ref usvort
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
  .reinit_turb=    1,
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
double cs_turb_dpow = -1.;

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
double cs_turb_sigmae = -1.;

/*!
 * Constant \f$C_1\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_crij1 = 1.80;

/*
 * Constant \f$C_2\f$ for the \f$R_{ij}-\varepsilon\f$ LRR model.
 * Useful if and only if \ref iturb=30 (\f$R_{ij}-\varepsilon\f$ LRR).
 */
double cs_turb_crij2 = 0.60;

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
double cs_turb_csaw1 = -1.;

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
 * \f$ \delta = \text{xlesfl} . (\text{ales} . volume)^{\text{bles}}\f$
 * \ref cs_turb_xlesfl is a constant used to define, for
 * each cell \f$\Omega_i\f$, the width of the (implicit) filter:
 * \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$\n
 * Useful if and only if \ref iturb = 40 or 41
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
 * model.
 * Any larger value yielded by the calculation procedure of the dynamic model
 * will be clipped to \f$ smagmx\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
double cs_turb_smagmx = -1.;

/*!
 * Minimum allowed value for the variable \f$C\f$ appearing in the LES dynamic
 * model.
 * Any smaller value yielded by the calculation procedure of the dynamic model
 * will be clipped to \f$ smagmn\f$.
 *
 * Useful if and only if \ref iturb = 41.
 */
double cs_turb_smagmn = 0.;

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
                             int     **itytur);

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
 *----------------------------------------------------------------------------*/

void
cs_f_turb_model_get_pointers(int     **iturb,
                             int     **itytur)
{
  *iturb  = &(_turb_model.iturb);
  *itytur = &(_turb_model.itytur);
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_turb_model
 *
 * needed to initialize structure with GUI
 */
/*----------------------------------------------------------------------------*/

cs_turb_model_t *
cs_get_glob_turb_model(void)
{
  return &_turb_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute turbulence model constants,
 *        some of which may depend on the model choice.
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_compute_constants(void)
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
  cs_turb_smagmx = cs_turb_csmago*cs_turb_csmago;
  cs_turb_smagmn = 0.;

  /* LRR constants */
  cs_turb_crij1 = 1.80;
  cs_turb_crij2 = 0.60;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_turb_ref_values
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
 * \brief Provide acces to cs_glob_turb_rans_model
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
 * \brief Provide acces to cs_glob_turb_les_model
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
 * \brief Print the turbulence model parameters to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_turb_model_log_setup(void)
{
  cs_var_cal_opt_t var_cal_opt;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "Turbulence model options\n"
       "------------------------\n\n"));

  cs_log_printf
    (CS_LOG_SETUP,
     _("  Continuous phase:\n\n"
       "    iturb :      %14d (Turbulence model)\n"
       "    iwallf:      %14d (wall function)\n"
       "                                (0: disabled)\n"
       "                                (1: one scale power law\n"
       "                                (forbidden for k-epsilon))\n"
       "                                (2: one scale log law)\n"
       "                                (3: two scales log law)\n"
       "                                (4: scalable wall function)\n"
       "                                (5: two scales V. Driest)\n"
       "                                (6: two scales smooth/rough)\n"
       "    iwallt:      %14d (Exch. coeff. correlation)\n"
       "                                (0: not activated)\n"
       "                                (1: activated)\n"
       "    ypluli:      %14.5e (Limit Y+)\n"
       "    igrhok:      %14d (1: computed Grad(rho k)\n\n"),
       cs_glob_turb_model->iturb,
       cs_glob_wall_functions->iwallf,
       cs_glob_wall_functions->iwallt,
       cs_glob_wall_functions->ypluli,
       cs_glob_turb_rans_model->igrhok);

  if (cs_glob_turb_model->iturb == 10) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Mixing length       (iturb = 10)\n"
         "    xlomlg:      %14.5e (Characteristic length)\n"),
         cs_glob_turb_rans_model->xlomlg);
  } else if (cs_glob_turb_model->iturb == 20) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   k-epsilon           (iturb = 20)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    iclkep:      %14d (k-epsilon clipping model)\n"
         "    ikecou:      %14d (k-epsilon coupling mode)\n"
         "    igrake:      %14d (Account for gravity)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->iclkep,
         cs_glob_turb_rans_model->ikecou,
         cs_glob_turb_rans_model->igrake);

    if (cs_glob_turb_rans_model->ikecou == 0 &&
        cs_glob_time_step_options->idtvar >= 0) {
      cs_real_t relaxvk, relaxve;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      relaxve = var_cal_opt.relaxv;
      cs_log_printf
        (CS_LOG_SETUP,
         _("    relaxv:      %14.5e for k (Relaxation)\n"
           "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
           relaxvk, relaxve);
    } else {
      cs_log_printf(CS_LOG_SETUP,_("\n"));
    }
  } else if (cs_glob_turb_model->iturb == 21) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Linear production k-epsilon (iturb = 21)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    iclkep:      %14d (k-epsilon clipping model)\n"
         "    ikecou:      %14d (k-epsilon coupling mode)\n"
         "    igrake:      %14d (Account for gravity)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->iclkep,
         cs_glob_turb_rans_model->ikecou,
         cs_glob_turb_rans_model->igrake);

    if (cs_glob_turb_rans_model->ikecou == 0 &&
        cs_glob_time_step_options->idtvar >= 0) {
      cs_real_t relaxvk, relaxve;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      relaxve = var_cal_opt.relaxv;
      cs_log_printf
        (CS_LOG_SETUP,
         _("    relaxv:      %14.5e for k (Relaxation)\n"
           "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
           relaxvk, relaxve);
    } else {
      cs_log_printf(CS_LOG_SETUP,_("\n"));
    }
  } else if (cs_glob_turb_model->iturb == 30) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Rij-epsilon         (iturb = 30)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    irijco:      %14d (Coupled resolution)\n"
         "    irijnu:      %14d (Matrix stabilization)\n"
         "    irijrb:      %14d (Reconstruct at boundaries)\n"
         "    irijec:      %14d (Wall echo terms)\n"
         "    idifre:      %14d (Handle diffusion tensor)\n"
         "    igrari:      %14d (Account for gravity)\n"
         "    iclsyr:      %14d (Symmetry implicitation)\n"
         "    iclptr:      %14d (Wall implicitation)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->irijco,
         cs_glob_turb_rans_model->irijnu,
         cs_glob_turb_rans_model->irijrb,
         cs_glob_turb_rans_model->irijec,
         cs_glob_turb_rans_model->idifre,
         cs_glob_turb_rans_model->igrari,
         cs_glob_turb_rans_model->iclsyr,
         cs_glob_turb_rans_model->iclptr);
  } else if (cs_glob_turb_model->iturb == 31) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   SSG Rij-epsilon     (iturb = 31)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    irijco:      %14d (Coupled resolution)\n"
         "    irijnu:      %14d (Matrix stabilization)\n"
         "    irijrb:      %14d (Reconstruct at boundaries)\n"
         "    igrari:      %14d (Account for gravity)\n"
         "    iclsyr:      %14d (Symmetry implicitation)\n"
         "    iclptr:      %14d (Wall implicitation)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->irijco,
         cs_glob_turb_rans_model->irijnu,
         cs_glob_turb_rans_model->irijrb,
         cs_glob_turb_rans_model->igrari,
         cs_glob_turb_rans_model->iclsyr,
         cs_glob_turb_rans_model->iclptr);
  } else if (cs_glob_turb_model->iturb == 32) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Rij-epsilon EBRSM     (iturb = 32)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    reinit_                      (Reinitialization of the\n"
         "     turb:       %14d  turbulence)\n"
         "    irijco:      %14d (Coupled resolution)\n"
         "    irijnu:      %14d (Matrix stabilization)\n"
         "    irijrb:      %14d (Reconstruct at boundaries)\n"
         "    igrari:      %14d (Account for gravity)\n"
         "    iclsyr:      %14d (Symmetry implicitation)\n"
         "    iclptr:      %14d (Wall implicitation)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->reinit_turb,
         cs_glob_turb_rans_model->irijco,
         cs_glob_turb_rans_model->irijnu,
         cs_glob_turb_rans_model->irijrb,
         cs_glob_turb_rans_model->igrari,
         cs_glob_turb_rans_model->iclsyr,
         cs_glob_turb_rans_model->iclptr);
  } else if (cs_glob_turb_model->itytur == 4) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   LES                 (iturb = 40, 41, 42)\n"
         "                            (Sub-grid scale model)\n"
         "                            (40 Smagorinsky model)\n"
         "                            (41 Dynamic model)\n"
         "                            (42 WALE model)\n"
         "    csmago:      %14.5e (Smagorinsky constant)\n"
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
         "                                (dynamic model case)\n"
         "    ivrtex:      %14d (Use of the vortex method)\n"),
         cs_turb_csmago, cs_turb_cwale, cs_turb_xlesfl,
         cs_turb_ales, cs_turb_bles, cs_glob_turb_les_model->idries,
         cs_turb_cdries, cs_turb_xlesfd, cs_turb_smagmx,
         cs_glob_turb_les_model->ivrtex);

  } else if (cs_glob_turb_model->iturb == 50) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   v2f phi-model       (iturb = 50)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    iclkep:      %14d (k-epsilon clipping model)\n"
         "    ikecou:      %14d (k-epsilon coupling mode)\n"
         "    igrake:      %14d (Account for gravity)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->iclkep,
         cs_glob_turb_rans_model->ikecou,
         cs_glob_turb_rans_model->igrake);
    if (cs_glob_turb_rans_model->ikecou == 0 &&
        cs_glob_time_step_options->idtvar >= 0) {
      cs_real_t relaxvk, relaxve;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      relaxve = var_cal_opt.relaxv;
      cs_log_printf
        (CS_LOG_SETUP,
         _("    relaxv:      %14.5e for k (Relaxation)\n"
           "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
           relaxvk, relaxve);
    } else {
      cs_log_printf(CS_LOG_SETUP,_("\n"));
    }
  } else if (cs_glob_turb_model->iturb == 51) {
    cs_log_printf
     (CS_LOG_SETUP,
      _("   v2f BL-v2/k         (iturb = 51)\n"
        "    almax:       %14.5e (Characteristic length)\n"
        "    uref:        %14.5e (Characteristic velocity)\n"
        "    iclkep:      %14d (k-epsilon clipping model)\n"
        "    ikecou:      %14d (k-epsilon coupling mode)\n"
        "    igrake:      %14d (Account for gravity)\n"),
        cs_glob_turb_ref_values->almax,
        cs_glob_turb_ref_values->uref,
        cs_glob_turb_rans_model->iclkep,
        cs_glob_turb_rans_model->ikecou,
        cs_glob_turb_rans_model->igrake);
    if (cs_glob_turb_rans_model->ikecou == 0 &&
        cs_glob_time_step_options->idtvar >= 0) {
      cs_real_t relaxvk, relaxve;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &var_cal_opt);
      relaxve = var_cal_opt.relaxv;
      cs_log_printf
        (CS_LOG_SETUP,
         _("    relaxv:      %14.5e for k (Relaxation)\n"
           "    relaxv:      %14.5e for epsilon (Relaxation)\n"),
           relaxvk, relaxve);
    } else {
      cs_log_printf(CS_LOG_SETUP,_("\n"));
    }
  } else if (cs_glob_turb_model->iturb == 60) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   k-omega SST         (iturb = 60)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    ikecou:      %14d (k-epsilon coupling mode)\n"
         "    igrake:      %14d (Account for gravity)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         cs_glob_turb_rans_model->ikecou,
         cs_glob_turb_rans_model->igrake);
    if (cs_glob_turb_rans_model->ikecou == 0 &&
        cs_glob_time_step_options->idtvar >= 0) {
      cs_real_t relaxvk, relaxvo;
      cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &var_cal_opt);
      relaxvk = var_cal_opt.relaxv;
      cs_field_get_key_struct(CS_F_(omg), key_cal_opt_id, &var_cal_opt);
      relaxvo = var_cal_opt.relaxv;
      cs_log_printf
        (CS_LOG_SETUP,
         _("    relaxv:      %14.5e for k (Relaxation)\n"
           "    relaxv:      %14.5e for omega (Relaxation)\n"),
           relaxvk, relaxvo);
    } else {
      cs_log_printf(CS_LOG_SETUP,_("\n"));
    }
  } else if (cs_glob_turb_model->iturb == 70) {
    cs_field_get_key_struct(CS_F_(nusa), key_cal_opt_id, &var_cal_opt);
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Spalart-Allmaras    (iturb = 70)\n"
         "    almax:       %14.5e (Characteristic length)\n"
         "    uref:        %14.5e (Characteristic velocity)\n"
         "    relaxv:      %14.5e for nu (Relaxation)\n"),
         cs_glob_turb_ref_values->almax,
         cs_glob_turb_ref_values->uref,
         var_cal_opt.relaxv);
  }

  if (cs_glob_turb_model->itytur == 2
   || cs_glob_turb_model->itytur == 5
   || cs_glob_turb_model->iturb == 60
   || cs_glob_turb_model->iturb == 70) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Rotation/curvature correction\n"
         "    irccor:      %14d (0: desactivated)\n"
         "                                (1: activated)\n"),
         cs_glob_turb_rans_model->irccor);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the turbulent constants to setup.log.
 *
 *----------------------------------------------------------------------------*/

void
cs_turb_constants_log_setup(void)
{
  cs_log_printf
    (CS_LOG_SETUP,
     _("\nConstants\n\n"
       "    xkappa:      %14.5e (Von Karman constant)\n"
       "    cstlog:      %14.5e (U+=Log(y+)/kappa +cstlog)\n"
       "    apow:        %14.5e (U+=apow (y+)**bpow (W&W law))\n"
       "    bpow:        %14.5e (U+=apow (y+)**bpow (W&W law))\n\n"),
       cs_turb_xkappa, cs_turb_cstlog, cs_turb_apow, cs_turb_bpow);

  if (cs_glob_turb_model->iturb == 20) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   k-epsilon           (iturb = 20)\n"
         "    ce1:         %14.5e (Cepsilon 1: production coef.)\n"
         "    ce2:         %14.5e (Cepsilon 2: dissipat.  coef.)\n"
         "    sigmak:      %14.5e (Prandtl relative to k)\n"
         "    sigmae:      %14.5e (Prandtl relative to epsilon )\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_ce1, cs_turb_ce2, cs_turb_sigmak,
         cs_turb_sigmae, cs_turb_cmu);
  } else if (cs_glob_turb_model->iturb == 21) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Linear production k-epsilon (iturb = 21)\n"
         "    ce1:         %14.5e (Cepsilon 1: production coef.)\n"
         "    ce2:         %14.5e (Cepsilon 2: dissipat.  coef.)\n"
         "    sigmak:      %14.5e (Prandtl relative to k)\n"
         "    sigmae:      %14.5e (Prandtl relative to epsilon )\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_ce1, cs_turb_ce2, cs_turb_sigmak,
         cs_turb_sigmae, cs_turb_cmu);
  } else if (cs_glob_turb_model->iturb == 30) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Rij-epsilon         (iturb = 30)\n"
         "    ce1:         %14.5e (Cepsilon 1: production coef.)\n"
         "    ce2:         %14.5e (Cepsilon 2: dissipat.  coef.)\n"
         "    crij1:       %14.5e (Slow term coefficient)\n"
         "    crij2:       %14.5e (Fast term coefficient)\n"
         "    crij3:       %14.5e (Gravity term coefficient)\n"
         "    sigmae:      %14.5e (sigma_eps coeff.)\n"
         "    csrij:       %14.5e (Rij diffusion coeff.)\n"
         "    crijp1:      %14.5e (Slow coeff. for wall echo)\n"
         "    crijp2:      %14.5e (Fast coeff. for wall echo)\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_ce1, cs_turb_ce2, cs_turb_crij1, cs_turb_crij2,
         cs_turb_crij3, cs_turb_sigmae, cs_turb_csrij, cs_turb_crijp1,
         cs_turb_crijp2, cs_turb_cmu);
  } else if (cs_glob_turb_model->iturb == 31) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   SSG Rij-epsilon     (iturb = 31)\n"
         "    cssgs1:      %14.5e (Cs1 coeff.)\n"
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
         "    sigmae:      %14.5e (sigma_eps coeff.)\n"
         "    cmu:         %14.5e (Cmu constant)\n"),
         cs_turb_cssgs1, cs_turb_cssgs2, cs_turb_cssgr1,
         cs_turb_cssgr2, cs_turb_cssgr3, cs_turb_cssgr4,
         cs_turb_cssgr5, cs_turb_csrij, cs_turb_crij3,
         cs_turb_ce1, cs_turb_cssge2, cs_turb_sigmae,
         cs_turb_cmu);
  } else if (cs_glob_turb_model->iturb == 32) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   EBRSM Rij-epsilon     (iturb = 32)\n"
         "    cebms1:      %14.5e (Cs1 coeff.)\n"
         "    cebmr1:      %14.5e (Cr1 coeff.)\n"
         "    cebmr2:      %14.5e (Cr2 coeff.)\n"
         "    cebmr3:      %14.5e (Cr3 coeff.)\n"
         "    cebmr4:      %14.5e (Cr4 coeff.)\n"
         "    cebmr5:      %14.5e (Cr5 coeff.)\n"
         "    csrij:       %14.5e (Rij Cs diffusion coeff.)\n"
         "    cebmr6:      %14.5e (Gravity term coeff.)\n"
         "    cebme2:      %14.5e (Coef Ceps2)\n"
         "    ce1:         %14.5e (Coef Ceps1)\n"
         "    sigmae:      %14.5e (Coef sigma_eps)\n"
         "    xa1:         %14.5e (Coef A1)\n"
         "    sigmak:      %14.5e (Coef sigma_k)\n"
         "    xceta:       %14.5e (Coef Ceta)\n"
         "    xct:         %14.5e (Coef CT)\n"),
         cs_turb_cebms1, cs_turb_cebmr1, cs_turb_cebmr2,
         cs_turb_cebmr3, cs_turb_cebmr4, cs_turb_cebmr5,
         cs_turb_csrij, cs_turb_cebmr6, cs_turb_cebme2,
         cs_turb_ce1, cs_turb_sigmae, cs_turb_xa1,
         cs_turb_sigmak, cs_turb_xceta, cs_turb_xct);
  } else if (cs_glob_turb_model->iturb == 50) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   v2f phi-model       (iturb = 50)\n"
         "    cv2fa1:      %14.5e (a1 to calculate Cepsilon1)\n"
         "    cv2fe2:      %14.5e (Cepsilon 2: dissip. coeff.)\n"
         "    sigmak:      %14.5e (Prandtl relative to k)\n"
         "    sigmae:      %14.5e (Prandtl relative to epsilon)\n"
         "    cv2fmu:      %14.5e (Cmu constant)\n"
         "    cv2fct:      %14.5e (CT constant)\n"
         "    cv2fcl:      %14.5e (CL constant)\n"
         "    cv2fet:      %14.5e (C_eta constant)\n"
         "    cv2fc1:      %14.5e (C1 constant)\n"
         "    cv2fc2:      %14.5e (C2 constant)\n"),
         cs_turb_cv2fa1, cs_turb_cv2fe2, cs_turb_sigmak,
         cs_turb_sigmae, cs_turb_cv2fmu, cs_turb_cv2fct,
         cs_turb_cv2fcl, cs_turb_cv2fet, cs_turb_cv2fc1,
         cs_turb_cv2fc2);
  } else if (cs_glob_turb_model->iturb == 51) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   v2f BL-v2/k         (iturb = 51)\n"
         "    cpale1:      %14.5e (Cepsilon 1 : Prod. coeff.)\n"
         "    cpale2:      %14.5e (Cepsilon 2 : Diss. coeff.)\n"
         "    cpale3:      %14.5e (Cepsilon 3 : E term coeff.)\n"
         "    cpale4:      %14.5e (Cepsilon 4 : Mod Diss. coef.)\n"
         "    sigmak:      %14.5e (Prandtl relative to k)\n"
         "    cpalse:      %14.5e (Prandtl relative to epsilon)\n"
         "    cpalmu:      %14.5e (Cmu constant)\n"
         "    cpalct:      %14.5e (CT constant)\n"
         "    cpalcl:      %14.5e (CL constant)\n"
         "    cpalet:      %14.5e (C_eta constant)\n"
         "    cpalc1:      %14.5e (C1 constant)\n"
         "    cpalc2:      %14.5e (C2 constant)\n"),
         cs_turb_cpale1, cs_turb_cpale2, cs_turb_cpale3,
         cs_turb_cpale4, cs_turb_sigmak, cs_turb_cpalse,
         cs_turb_cpalmu, cs_turb_cpalct, cs_turb_cpalcl,
         cs_turb_cpalet, cs_turb_cpalc1, cs_turb_cpalc2);
  } else if (cs_glob_turb_model->iturb == 60) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   k-omega SST         (iturb = 60)\n"
         "    ckwsk1:      %14.5e (sigma_k1 constant)\n"
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
  } else if (cs_glob_turb_model->iturb == 70) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Spalart-Allmaras    (iturb = 70)\n"
         "    csab1:        %14.5e (b1 constant)\n"
         "    csab2:        %14.5e (b2 constant)\n"
         "    csasig:       %14.5e (sigma constant)\n"
         "    csav1:        %14.5e (v1 constant)\n"
         "    csaw1:        %14.5e (w1 constant)\n"
         "    csaw2:        %14.5e (w2 constant)\n"
         "    csaw3:        %14.5e (w3 constant)\n"),
         cs_turb_csab1, cs_turb_csab2, cs_turb_csasig,
         cs_turb_csav1, cs_turb_csaw1, cs_turb_csaw2,
         cs_turb_csaw3);
  }

  int iokss = 0, iokcaz = 0;

  if (cs_glob_turb_rans_model->irccor == 1) {
    if (cs_glob_turb_rans_model->itycor == 1)
      iokcaz = 1;
    else if (cs_glob_turb_rans_model->itycor == 2)
      iokss = 1;
  }

  if (iokcaz > 0) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Rotation/curvature correction (Cazalbou)\n"
         "    ccaze2:       %14.5e (Coef Ce2^0)\n"
         "    ccazsc:       %14.5e (Coef Csc)\n"
         "    ccaza:        %14.5e (Coef a)\n"
         "    ccazb:        %14.5e (Coef b)\n"
         "    ccazc:        %14.5e (Coef c)\n"
         "    ccazd:        %14.5e (Coef d)\n"),
         cs_turb_ccaze2, cs_turb_ccazsc, cs_turb_ccaza,
         cs_turb_ccazb, cs_turb_ccazc, cs_turb_ccazd);
  }

  if (iokss > 0) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Rotation/curvature correction (Spalart-Shur)\n"
         "    cssr1:       %14.5e (Coef c_r1)\n"
         "    cssr2:       %14.5e (Coef c_r2)\n"
         "    cssr3:       %14.5e (Coef c_r3)\n"),
         cs_turb_cssr1, cs_turb_cssr2, cs_turb_cssr3);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

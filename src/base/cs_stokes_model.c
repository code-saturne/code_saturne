/*============================================================================
 * Stokes equation model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_stokes_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_stokes_model.c
        Stokes equation model data.

  \struct cs_stokes_model_t

  \brief Stokes equation model descriptor.

  Members of these Stokes equation model descriptor are publicly accessible, to
  allow for concise syntax, as it is expected to be used in many places.

  \var  cs_stokes_model_t::ivisse
        <a name="ivisse"></a>
        Indicates whether the source terms in transposed gradient
        and velocity divergence should be taken into account in the
        momentum equation. In the compressible module, these terms
        also account for the volume viscosity (cf. \ref ppincl::viscv0 "viscv0" and
        \ref ppincl::iviscv "iviscv")
        \f$\partial_i \left[(\kappa -2/3\,(\mu+\mu_t))\partial_k U_k  \right]
        +     \partial_j \left[ (\mu+\mu_t)\partial_i U_j \right]\f$:
        - 0: not taken into account,
        - 1: taken into account.
  \var  cs_stokes_model_t::irevmc
        reconstruction of the velocity field with the updated pressure option
        - 0: standard gradient of pressure increment (default)
  \var  cs_stokes_model_t::iprco
        compute the pressure step thanks to the continuity equation
        - 1: true (default)
        - 0: false
  \var  cs_stokes_model_t::arak
        <a name="arak"></a>
        Arakawa multiplicator for the Rhie and Chow filter (1 by default).\n\n
        Please refer to the
        <a href="../../theory.pdf#arak"><b>Rhie and Chow filter</b></a> section
        of the theory guide for more informations.
  \var  cs_stokes_model_t::ipucou
        indicates the algorithm for velocity/pressure coupling:
        - 0: standard algorithm,
        - 1: reinforced coupling in case calculation with long time steps\n
        Always useful (it is seldom advised, but it can prove very useful,
        for instance, in case of flows with weak convection effects and
        highly variable viscosity).
  \var  cs_stokes_model_t::iccvfg
        indicates whether the dynamic field should be frozen or not:
           - 1: true
           - 0: false (default)\n
        In such a case, the values of velocity, pressure and the
        variables related to the potential turbulence model
        (\f$k\f$, \f$R_{ij}\f$, \f$\varepsilon\f$, \f$\varphi\f$,
        \f$\bar{f}\f$, \f$\omega\f$, turbulent viscosity) are kept
        constant over time and only the equations for the scalars
        are solved.\n Also, if \ref iccvfg = 1, the physical properties
        modified in \ref cs_user_physical_properties will keep being
        updated. Beware of non-consistencies if these properties would
        normally affect the dynamic field (modification of density for
        instance).\n Useful if and only if \ref dimens::nscal "nscal"
        \f$>\f$ 0 and the calculation is a restart.
  \var  cs_stokes_model_t::idilat
        algorithm to take into account the density variation in time
        - 1: dilatable steady algorithm (default)
        - 2: dilatable unsteady algorithm
        - 3: low-Mach algorithm
        - 4: algorithm for fire
  \var  cs_stokes_model_t::epsdp
        parameter of diagonal pressure strengthening
  \var  cs_stokes_model_t::itbrrb
        accurate treatment of the wall temperature
        (reconstruction of wall temperature)
        - 1: true
        - 0: false (default)
        (see \ref condli, useful in case of coupling with syrthes)
  \var  cs_stokes_model_t::iphydr
        <a name="iphydr"></a>
        improve static pressure algorithm
        Take into account the balance or imbalance between the pressure
        gradient and source terms (as gravity and head losses)
        - 1: impose the equilibrium of the static part of the pressure with
          any external force, even head losses
        - 0: no treatment (default)
        - 2: hydrostatic pressure computation with a apriori momentum equation
             to obtain a hydrostatic pressure taking into account the imbalance
             between the pressure gradient and the gravity source term.\n\n
        When the density effects are important, the choice of \ref iphydr = 1
        allows to improve the interpolation of the pressure and correct the
        non-physical velocities which may appear in highly stratified areas
        or near horizontal walls (thus avoiding the use of
        \ref cs_var_cal_opt_t::extrag "extrag"
        if the non-physical velocities are due only to gravity effects).\n
        The improved algorithm also allows eradicating the velocity oscillations
        which tend to appear at the frontiers of areas with high head losses.\n
        In the case of a stratified flow, the calculation cost is higher when
        the improved algorithm is used (about 30\% depending on the case)
        because the hydrostatic pressure must be recalculated at the outlet
        boundary conditions: see \ref icalhy.\n
        On meshes of insufficient quality, in order to
        improve the convergence, it may be useful to increase the number of
        iterations for the reconstruction of the pressure right-hand side,
        i.e. \ref cs_var_cal_opt_t::nswrsm "nswrsm".\n If head losses are
        present just along an outlet boundary, it is necessary to specify
        \ref icalhy = 0 in order to deactivate the recalculation of the
        hydrostatic pressure at the boundary, which may otherwise cause
        instabilities. Please refer to the
        <a href="../../theory.pdf#iphydr"><b>handling of the hydrostatic pressure</b></a>
        section of the theory guide for more informations.
  \var  cs_stokes_model_t::igprij
        improve static pressure algorithm
        - 1: take -div(rho R) in the static pressure
          treatment IF iphydr=1
        - 0: no treatment (default)
  \var  cs_stokes_model_t::igpust
        improve static pressure algorithm
        - 1: take user momentum source terms in the static pressure
          treatment IF iphydr=1
        - 0: no treatment (default)
  \var  cs_stokes_model_t::iifren
        indicates the presence of a Bernoulli boundary face (automatically
        computed)
        - 0: no face
        - 1: at least one face
  \var  cs_stokes_model_t::icalhy
        compute the hydrostatic pressure in order to compute the Dirichlet
        conditions on the pressure at outlets
        - 1: calculation of the hydrostatic pressure at the outlet boundary
        - 0: no calculation of the hydrostatic pressure at the outlet boundary (default)
        This option is automatically specified depending on the choice of
        \ref iphydr and the value of gravity (\ref icalhy = 1 if  \ref iphydr = 1
        and gravity is different from 0; otherwise \ref icalhy = 0). The
        activation of this option generates an additional calculation cost
        (about 30\% depending on the case).\n If head losses are present
        just along an outlet boundary, it is necessary to specify \ref icalhy = 0
        in order to deactivate the recalculation of the hydrostatic pressure
        at the boundary, which may otherwise cause instabilities
  \var  cs_stokes_model_t::irecmf
        use interpolated face diffusion coefficient instead of cell diffusion coefficient
        for the mass flux reconstruction for the non-orthogonalities
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

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main Stokes equation model descriptor structure and associated pointer:
 * Default Options */

static cs_stokes_model_t  _stokes_model = {
  .ivisse = 1,
  .irevmc = 0,
  .iprco  = 1,
  .arak   = 1.0,
  .ipucou = 0,
  .iccvfg = 0,
  .idilat = 1,
  .epsdp  = 1.e-12,
  .itbrrb = 0,
  .iphydr = 0,
  .igprij = 0,
  .igpust = 1,
  .iifren = 0,
  .icalhy = -1,
  .irecmf = 0};

const cs_stokes_model_t  *cs_glob_stokes_model = &_stokes_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_stokes_options_get_pointers(int     **ivisse,
                                 int     **irevmc,
                                 int     **iprco,
                                 double  **arak,
                                 int     **ipucou,
                                 int     **iccvfg,
                                 int     **idilat,
                                 double  **epsdp,
                                 int     **itbrrb,
                                 int     **iphydr,
                                 int     **igprij,
                                 int     **igpust,
                                 int     **iifren,
                                 int     **icalhy,
                                 int     **irecmf);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global Stokes model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ivisse  --> pointer to cs_glob_stokes_model->ivisse
 *   irevmc  --> pointer to cs_glob_stokes_model->irevmc
 *   iprco   --> pointer to cs_glob_stokes_model->iprco
 *   arak    --> pointer to cs_glob_stokes_model->arak
 *   ipucou  --> pointer to cs_glob_stokes_model->ipucou
 *   iccvfg  --> pointer to cs_glob_stokes_model->iccvfg
 *   idilat  --> pointer to cs_glob_stokes_model->idilat
 *   epsdp   --> pointer to cs_glob_stokes_model->epsdp
 *   itbrrb  --> pointer to cs_glob_stokes_model->itbrrb
 *   iphydr  --> pointer to cs_glob_stokes_model->iphydr
 *   igprij  --> pointer to cs_glob_stokes_model->igprij
 *   igpust  --> pointer to cs_glob_stokes_model->igpust
 *   iifren  --> pointer to cs_glob_stokes_model->iifren
 *   icalhy  --> pointer to cs_glob_stokes_model->icalhy
 *   irecmf  --> pointer to cs_glob_stokes_model->irecmf
 *----------------------------------------------------------------------------*/

void
cs_f_stokes_options_get_pointers(int     **ivisse,
                                 int     **irevmc,
                                 int     **iprco,
                                 double  **arak,
                                 int     **ipucou,
                                 int     **iccvfg,
                                 int     **idilat,
                                 double  **epsdp,
                                 int     **itbrrb,
                                 int     **iphydr,
                                 int     **igprij,
                                 int     **igpust,
                                 int     **iifren,
                                 int     **icalhy,
                                 int     **irecmf)
{
  *ivisse = &(_stokes_model.ivisse);
  *irevmc = &(_stokes_model.irevmc);
  *iprco  = &(_stokes_model.iprco);
  *arak   = &(_stokes_model.arak);
  *ipucou = &(_stokes_model.ipucou);
  *iccvfg = &(_stokes_model.iccvfg);
  *idilat = &(_stokes_model.idilat);
  *epsdp  = &(_stokes_model.epsdp );
  *itbrrb = &(_stokes_model.itbrrb);
  *iphydr = &(_stokes_model.iphydr);
  *igprij = &(_stokes_model.igprij);
  *igpust = &(_stokes_model.igpust);
  *iifren = &(_stokes_model.iifren);
  *icalhy = &(_stokes_model.icalhy);
  *irecmf = &(_stokes_model.irecmf);
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_stokes_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_stokes_model_t *
cs_get_glob_stokes_model(void)
{
  return &_stokes_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the stokes model parameters to setup log.
 */
/*----------------------------------------------------------------------------*/

void
cs_stokes_model_log_setup(void)
{
  if (cs_glob_field_pointers == NULL)
    return;

  cs_var_cal_opt_t var_cal_opt;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  cs_field_t *f_pot = NULL;
  if (cs_glob_physical_model_flag[CS_GROUNDWATER] > 0)
    f_pot = CS_F_(head);
  else
    f_pot = CS_F_(p);

  if (f_pot == NULL)
    return;

  const char *f_pot_label = cs_field_get_label(f_pot);

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "Secondary viscosity\n"
       "-------------------\n\n"
       "   Continuous phase:\n\n"
       "    ivisse:      %14d (1: accounted for)\n\n"),
     cs_glob_stokes_model->ivisse);

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "Stokes model\n"
       "------------\n\n"
       "    idilat:      %14d (1: without unsteady term\n"
       "                                    in the continuity equation\n"
       "                                 2: with unsteady term in\n"
       "                                    the continuity equation\n"
       "                                 3 : with unsteady term in\n"
       "                                     the continuity equation\n"
       "                                     and a thermo pressure\n"
       "                                     constant in the domain\n"
       "                                 4 : with unsteady term in\n"
       "                                and  the continuity equation\n"
       "                                 5   for fire modelling)\n"

       "    iporos:      %14d (0: without porous media\n"
       "                                 1: with porous media \n"
       "                                 2: with tensorial porous media\n"
       "                                 3: with intergal formulation\n"
       "                                    including fluid volumes and\n"
       "                                    fluid surfaces)\n"
       "    iphydr:      %14d (1: account for explicit\n"
       "                                    balance between pressure\n"
       "                                    gradient, gravity source\n"
       "                                    terms, and head losses\n"
       "                                  2: compute a hydrostatic\n"
       "                                     pressure which is balance\n"
       "                                     in balance with buoyancy)\n"
       "    icalhy:      %14d (1: compute hydrostatic\n"
       "                                    pressure for dirichlet\n"
       "                                    conditions for pressure\n"
       "                                    on outlet)\n"
       "    iprco :      %14d (1: pressure-continuity)\n"
       "    ipucou:      %14d (1: reinforced u-p coupling)\n"
       "    nterup:      %14d (n: n sweeps on navsto for\n"
       "                                    velocity/pressure coupling)\n"),
     cs_glob_stokes_model->idilat,
     cs_glob_porous_model,
     cs_glob_stokes_model->iphydr,
     cs_glob_stokes_model->icalhy,
     cs_glob_stokes_model->iprco,
     cs_glob_stokes_model->ipucou,
     cs_glob_piso->nterup);

  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "   Continuous phase:\n\n"
       "    irevmc:      %14d (Velocity reconstruction mode)\n"),
     cs_glob_stokes_model->irevmc);

  if (cs_glob_time_step_options->idtvar >= 0) {
    cs_field_get_key_struct(f_pot, key_cal_opt_id, &var_cal_opt);
    cs_log_printf
      (CS_LOG_SETUP,
       _("    relaxv:      %14.5e for %s (relaxation)\n"
         "    arak:        %14.5e (Arakawa factor)\n"),
       var_cal_opt.relaxv, f_pot_label, cs_glob_stokes_model->arak);
  } else {
    cs_field_get_key_struct(CS_F_(vel), key_cal_opt_id, &var_cal_opt);
    cs_log_printf
      (CS_LOG_SETUP,
       _("    arak:        %14.5e (Arakawa factor)\n"),
       var_cal_opt.relaxv * cs_glob_stokes_model->arak);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

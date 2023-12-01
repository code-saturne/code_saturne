/*============================================================================
 * Velocity-pressure coupling model and parameters.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_mesh_location.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_velocity_pressure.c
        Velocity-pressure coupling model and parameters.
*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_velocity_pressure_model_t

  \brief Stokes equation model descriptor.

  Members of these Stokes equation model descriptor are publicly accessible, to
  allow for concise syntax, as it is expected to be used in many places.

  \var  cs_velocity_pressure_model_t::ivisse
        <a name="ivisse"></a>
        Indicates whether the source terms in transposed gradient
        and velocity divergence should be taken into account in the
        momentum equation. In the compressible module, these terms
        also account for the volume viscosity (cf. \ref ppincl::viscv0 "viscv0"
        and \ref ppincl::iviscv "iviscv")
        \f$\partial_i \left[(\kappa -2/3\,(\mu+\mu_t))\partial_k U_k  \right]
        +     \partial_j \left[ (\mu+\mu_t)\partial_i U_j \right]\f$:
        - 0: not taken into account,
        - 1: taken into account.

  \var  cs_velocity_pressure_model_t::idilat
        algorithm to take into account the density variation in time
        - 0: Boussinesq approximation (rho constant except in the buoyant
             term where \f$\Delta \rho \vect{g} = - \rho \beta \Delta T \vect{g} \f$
        - 1: dilatable steady algorithm (default)
        - 2: dilatable unsteady algorithm
        - 3: low-Mach algorithm
        - 4: algorithm for fire

  \var  cs_velocity_pressure_model_t::fluid_solid
        Has a solid zone where dynamics must be killed?
        - false (default)
        - true
  \var  cs_velocity_pressure_model_t::iprcdo
        Indicate discretization method for presssure
        - 0: Use Legacy Finite Volume method
        - 1: Use CDO method, if CDO/FV is coupled

*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_velocity_pressure_param_t

  \brief Inner velocity/pressure iteration options descriptor.

  Members of this structure are publicly accessible, to allow for
  concise  syntax, as they are expected to be used in many places.

  \var  cs_velocity_pressure_param_t::iphydr
        <a name="iphydr"></a>
        Improved pressure interpolation scheme.
        Take into account the balance or imbalance between the pressure
        gradient and source terms (as gravity and head losses)
        - 1: impose the equilibrium of the static part of the pressure with
          any external force, even head losses (default)
        - 2: hydrostatic pressure computation with an apriori momentum equation
             to obtain a hydrostatic pressure taking into account the imbalance
             between the pressure gradient and the gravity source term.
        - 0: no treatment\n\n
        When the density effects are important, the choice of \ref iphydr = 1
        allows to improve the interpolation of the pressure and correct the
        non-physical velocities which may appear in highly stratified areas
        or near horizontal walls.\n
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
        section of the theory guide for more information.\n
        The iphydr = 2 option is a legacy treatment to improve the computation
        of the pressure gradient for buoyant/stratified flows. In most cases,
        iphydr = 2 is equivalent to iphydr = 1, but for the following situations,
        iphydr = 2 can yield better results:
        - multiple inlet/outlets with different altitudes
        - outlets normal to the gravity
        Note that iphydr = 2 is less general than iphydr = 1: only gravity forces
        are taken into account.\n

  \var  cs_velocity_pressure_param_t::icalhy
        compute the hydrostatic pressure in order to compute the Dirichlet
        conditions on the pressure at outlets
        - 1: calculation of the hydrostatic pressure at the outlet boundary
        - 0: no calculation of the hydrostatic pressure at the outlet boundary
             (default)
        This option is automatically specified depending on the choice of
        \ref iphydr and the value of gravity (\ref icalhy = 1 if  \ref iphydr = 1
        and gravity is different from 0; otherwise \ref icalhy = 0). The
        activation of this option generates an additional calculation cost
        (about 30\% depending on the case).\n If head losses are present
        just along an outlet boundary, it is necessary to specify \ref icalhy = 0
        in order to deactivate the recalculation of the hydrostatic pressure
        at the boundary, which may otherwise cause instabilities

  \var  cs_velocity_pressure_param_t::iprco
        compute the pressure step thanks to the continuity equation
        - 1: true (default)
        - 0: false

  \var  cs_velocity_pressure_param_t::ipredfl
        Switch on mass flux prediction before momentum solving to be fully
        conservative in momentum over time for variable density flows.
        \deprecated Will be removed in a future version.

  \var  cs_velocity_pressure_param_t::irevmc
        reconstruction of the velocity field with the updated pressure option
        - 0: standard gradient of pressure increment (default)

  \var  cs_velocity_pressure_param_t::iifren
        indicates the presence of a Bernoulli boundary face (automatically
        computed)
        - 0: no face
        - 1: at least one face

  \var  cs_velocity_pressure_param_t::irecmf
        use interpolated face diffusion coefficient instead of cell diffusion
        coefficient for the mass flux reconstruction for the
        non-orthogonalities
        - 1: true
        - 0: false (default)

  \var  cs_velocity_pressure_param_t::igprij
        improve static pressure algorithm
        - 1: take -div(rho R) in the static pressure
          treatment IF iphydr=1
        - 0: no treatment (default)

  \var  cs_velocity_pressure_param_t::igpust
        Improved pressure interpolation scheme:
        - 1: take user momentum source terms in the static pressure
          treatment IF iphydr=1
        - 0: no treatment (default)

  \var  cs_velocity_pressure_param_t::igrdpp
        For the compressible algorithm, indicate whether the pressure
        should be updated after the solution of the acoustic equation.
        - 1: true (default)
        - 0: false

  \var  cs_velocity_pressure_param_t::ipucou
        indicates the algorithm for velocity/pressure coupling:
        - 0: standard algorithm,
        - 1: reinforced coupling in case calculation with long time steps\n
        Always useful (it is seldom advised, but it can prove very useful,
        for instance, in case of flows with weak convection effects and
        highly variable viscosity).

  \var  cs_velocity_pressure_param_t::itpcol
        Time scheme option:
        - 0: staggered time scheme. On the time grids, the velocity is
             half a time step behind the density and the buoyant scalar.
             (See the thesis of \cite Pierce:2004)
        - 1: collocated time scheme. On the time grids, the velocity is
             at the same location as the density and the buoyant scalar.
             (See \cite Ma:2019)

  \var  cs_velocity_pressure_param_t::arak
        <a name="arak"></a>
        Arakawa multiplicator for the Rhie and Chow filter (1 by default).\n\n
        Please refer to the
        <a href="../../theory.pdf#arak"><b>Rhie and Chow filter</b></a> section
        of the theory guide for more informations.

  \var  cs_velocity_pressure_param_t::rcfact
        <a name="rcfact"></a>
        Factor of the Rhie and Chow filter:\n
        - 0: dt (by default),\n
        - 1: 1/A_u.\n

  \var  cs_velocity_pressure_param_t::staggered
        <a name="staggered"></a>
        1D zone simulator option:\n
         - 0: colocated.\n
         - 1: staggered.\n

  \var  cs_velocity_pressure_param_t::nterup
        number of iterations on the pressure-velocity coupling on Navier-Stokes

  \var  cs_velocity_pressure_param_t::epsup
        relative precision for the convergence test of the iterative process on
        pressure-velocity coupling

  \var  cs_velocity_pressure_param_t::xnrmu
        norm  of the increment \f$ \vect{u}^{k+1} - \vect{u}^k \f$ of the
        iterative process on pressure-velocity coupling

  \var  cs_velocity_pressure_param_t::xnrmu0
        norm of \f$ \vect{u}^0 \f$

  \var  cs_velocity_pressure_param_t::epsdp
        parameter of diagonal pressure strengthening
*/

/*----------------------------------------------------------------------------*/

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

static cs_velocity_pressure_model_t  _velocity_pressure_model = {
  .ivisse = 1,
  .idilat = 1,
  .fluid_solid = false,
  .n_buoyant_scal = 0,
  .iprcdo = 0,
};

const cs_velocity_pressure_model_t  *cs_glob_velocity_pressure_model
  = &_velocity_pressure_model;

/* Velocity/pressure inner iterations structure and associated pointer */

static cs_velocity_pressure_param_t  _velocity_pressure_param =
{
  .iphydr = 1,
  .icalhy = -1,
  .iprco  = 1,
  .ipredfl = 0,
  .irevmc = 0,
  .iifren = 0,
  .irecmf = 0,
  .igprij = 0,
  .igpust = 1,
  .igrdpp = 1,
  .ipucou = 0,
  .itpcol = -1,
  .arak   = 1.0,
  .rcfact = 0,
  .staggered = 0,
  .nterup = 1,
  .epsup = 1e-5,
  .xnrmu = 0.,
  .xnrmu0 = 0.,
  .epsdp  = 1.e-12,
};

const cs_velocity_pressure_param_t  *cs_glob_velocity_pressure_param
  = &_velocity_pressure_param;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_velocity_pressure_model_get_pointers(int     **ivisse,
                                          int     **idilat,
                                          bool    **fluid_solid,
                                          int     **n_buoyant_scal,
                                          int     **iprcdo);
void
cs_f_velocity_pressure_param_get_pointers(int     **iphydr,
                                          int     **icalhy,
                                          int     **iprco,
                                          int     **ipredfl,
                                          int     **irevmc,
                                          int     **iifren,
                                          int     **irecmf,
                                          int     **igprij,
                                          int     **igpust,
                                          int     **ipucou,
                                          int     **itpcol,
                                          double  **arak,
                                          int     **rcfact,
                                          int     **staggered,
                                          int     **nterup,
                                          double  **epsup,
                                          double  **xnrmu,
                                          double  **xnrmu0,
                                          double  **epsdp);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global Stokes model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_velocity_pressure_model_get_pointers(int     **ivisse,
                                          int     **idilat,
                                          bool    **fluid_solid,
                                          int     **n_buoyant_scal,
                                          int     **iprcdo)
{
  *ivisse = &(_velocity_pressure_model.ivisse);
  *idilat = &(_velocity_pressure_model.idilat);
  *fluid_solid = &(_velocity_pressure_model.fluid_solid);
  *n_buoyant_scal = &(_velocity_pressure_model.n_buoyant_scal);
  *iprcdo = &(_velocity_pressure_model.iprcdo);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global velocity_pressure structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_velocity_pressure_param_get_pointers(int     **iphydr,
                                          int     **icalhy,
                                          int     **iprco,
                                          int     **ipredfl,
                                          int     **irevmc,
                                          int     **iifren,
                                          int     **irecmf,
                                          int     **igprij,
                                          int     **igpust,
                                          int     **ipucou,
                                          int     **itpcol,
                                          double  **arak,
                                          int     **rcfact,
                                          int     **staggered,
                                          int     **nterup,
                                          double  **epsup,
                                          double  **xnrmu,
                                          double  **xnrmu0,
                                          double  **epsdp)
{
  *iphydr = &(_velocity_pressure_param.iphydr);
  *icalhy = &(_velocity_pressure_param.icalhy);
  *iprco  = &(_velocity_pressure_param.iprco);
  *ipredfl = &(_velocity_pressure_param.ipredfl);
  *irevmc = &(_velocity_pressure_param.irevmc);
  *iifren = &(_velocity_pressure_param.iifren);
  *irecmf = &(_velocity_pressure_param.irecmf);
  *igprij = &(_velocity_pressure_param.igprij);
  *igpust = &(_velocity_pressure_param.igpust);
  *ipucou = &(_velocity_pressure_param.ipucou);
  *itpcol = &(_velocity_pressure_param.itpcol);
  *arak   = &(_velocity_pressure_param.arak);
  *rcfact = &(_velocity_pressure_param.rcfact);
  *staggered = &(_velocity_pressure_param.staggered);
  *nterup = &(_velocity_pressure_param.nterup);
  *epsup  = &(_velocity_pressure_param.epsup);
  *xnrmu  = &(_velocity_pressure_param.xnrmu);
  *xnrmu0 = &(_velocity_pressure_param.xnrmu0);
  *epsdp  = &(_velocity_pressure_param.epsdp );
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
 * \brief Provide read/write access to cs_glob_velocity_pressure_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_velocity_pressure_model_t *
cs_get_glob_velocity_pressure_model(void)
{
  return &_velocity_pressure_model;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cs_glob_velocity_pressure_param.
 *
 * Needed to initialize structure with GUI.
 *
 * \return  velocity_pressure information structure
 */
/*----------------------------------------------------------------------------*/

cs_velocity_pressure_param_t *
cs_get_glob_velocity_pressure_param(void)
{
  return &_velocity_pressure_param;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Count and set number of buoyant scalars.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_set_n_buoyant_scalars(void)
{
  const int n_fields = cs_field_n_fields();
  const int key_sca = cs_field_key_id("scalar_id");
  const int key_buo = cs_field_key_id("is_buoyant");

  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (   f->type & CS_FIELD_VARIABLE
        && cs_field_get_key_int(f, key_sca) > -1) {
      if (cs_field_get_key_int(f, key_buo)) {
        _velocity_pressure_model.n_buoyant_scal += 1;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 *!
 * \brief Set `fluid_solid` flag if solid zones are present.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_set_solid(void)
{
  const int n_zones = cs_volume_zone_n_zones();

  for (int id = 0; id < n_zones; id++) {
    const cs_zone_t  *z = cs_volume_zone_by_id(id);
    if (z->type & CS_VOLUME_ZONE_SOLID) {
      /* Activate the solid flag */
      _velocity_pressure_model.fluid_solid = true;
      break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the volocity-pressure model parameters to setup log.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_model_log_setup(void)
{
  if (cs_glob_field_pointers == NULL)
    return;

  const cs_velocity_pressure_model_t *vp_model = cs_glob_velocity_pressure_model;

  cs_field_t *f_p = NULL;
  if (cs_glob_physical_model_flag[CS_GROUNDWATER] > 0)
    f_p = CS_F_(head);
  else
    f_p = CS_F_(p);

  if (f_p == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "Velocity-pressure model\n"
                 "-----------------------\n"));

  const char *ivisse_value_str[] = {N_("0 (ignored)"),
                                    N_("1 (taken into account)")};

  cs_log_printf(CS_LOG_SETUP,
                _("\n  Viscous term of transposed velocity gradient:\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    ivisse:        %s\n\n"),
                _(ivisse_value_str[vp_model->ivisse]));

  cs_log_printf(CS_LOG_SETUP,
                _("\n  Variable density / dilatable model:\n"));
  const char *idilat_value_str[]
    = {N_("0 (Boussinesq approximation)"),
       N_("1 (without unsteady term in the continuity equation)"),
       N_("2 (with unsteady term in the continuity equation)"),
       N_("3 (with unsteady term in the continuity equationnn\n"
          "                   "
          "   and a thermo pressure constant in the domain)"),
       N_("4 (with unsteady term in the continuity equation)"),
       N_("5 (for fire modelling)")};
  cs_log_printf(CS_LOG_SETUP,
                _("    idilat:        %s\n"),
                _(idilat_value_str[vp_model->idilat]));

  cs_log_printf(CS_LOG_SETUP,
                _("\n  Porosity model:\n"));
  const char *iporos_value_str[]
    = {N_("0 (without porous media)"),
       N_("1 (with porous media)"),
       N_("2 (with tensorial porous media)"),
       N_("3 (with integral formulation\n"
          "                   "
          "   including fluid volumes and fluid surfaces)")};
  cs_log_printf(CS_LOG_SETUP,
                _("    iporos:        %s\n"),
                _(iporos_value_str[cs_glob_porous_model]));

  if (vp_model->fluid_solid)
    cs_log_printf
      (CS_LOG_SETUP,
       _("\n"
         "  Fluid-solid mode (disable dynamics in the solid part)\n\n"));

  if (vp_model->iprcdo)
    cs_log_printf
      (CS_LOG_SETUP,
       _("\n"
         "  Pressure correction equation is solved by CDO\n\n"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print Velocity-pressure parameters to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_velocity_pressure_param_log_setup(void)
{
  cs_field_t *f_p = NULL;
  if (cs_glob_physical_model_flag[CS_GROUNDWATER] > 0)
    f_p = CS_F_(head);
  else
    f_p = CS_F_(p);

  if (f_p == NULL)
    return;

  const char *f_p_label = cs_field_get_label(f_p);

  const cs_velocity_pressure_param_t *vp_param
    = cs_glob_velocity_pressure_param;

  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "Velocity-pressure parameters\n"
                 "----------------------------\n\n"));

  cs_log_printf(CS_LOG_SETUP,
                ("    nterup:        %d (number of U-P sub iterations)\n"),
                vp_param->nterup);

  const char *iphydr_value_str[]
    = {N_("0 (no treatment (default) for the improvement of\n"
          "                   "
          "   static pressure algorithm)"),
       N_("1 (account for explicit balance between pressure\n"
          "                   "
          "   gradient, gravity source terms and head losses)"),
       N_("2 (compute a hydrostatic pressure which is\n"
          "                   "
          "   in balance with buoyancy)")};
  cs_log_printf(CS_LOG_SETUP,
                _("    iphydr:        %s\n"),
                _(iphydr_value_str[vp_param->iphydr]));

  /* Sub options of "iphydr=1" */
  if (vp_param->iphydr == 1) {

    const char *icalhy_value_str[]
      = {N_("0 ((default)\n"
            "                   "
            "   do not compute hydrostatic pressure for dirichlet\n"
            "                   "
            "   conditions for pressure on outlet)"),
         N_("1 (compute hydrostatic pressure for dirichlet\n"
            "                   "
            "   conditions for pressure on outlet)")};

    cs_log_printf(CS_LOG_SETUP,
                  _("    icalhy:        %s\n"),
                  _(icalhy_value_str[vp_param->icalhy]));

    const char *igpust_value_str[]
      = {N_("0 (no treatment for the improvment of static\n"
            "                   "
            "   pressure algorithm)"),
         N_("1 (take user momentum source terms into account\n"
            "                   "
            "   in the hydrostatic pressure computation)")};

    cs_log_printf(CS_LOG_SETUP,
                  _("    igpust:        %s\n"),
                  _(igpust_value_str[vp_param->igpust]));

    const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();
    if (turb_model != NULL) {
      if (turb_model->order == CS_TURB_SECOND_ORDER){
        const char *igprij_value_str[]
          = {N_("0 (do not take into account div(rho R) terms in the\n"
                "                   "
                "   hydrostatic pressure computation)"),
             N_("1 (take div(rho R) terms into account\n"
                "                   "
                "   in the hydrostatic pressure computation)")};

      cs_log_printf(CS_LOG_SETUP,
                    _("    igprij:        %s\n"),
                    _(igprij_value_str[vp_param->igprij]));
      }
    }
  }

  const char *iprco_value_str[]
    = {N_("0 (do not compute the pressure step\n"
          "                   "
          "   using the continuity equation)\n"),
       N_("1 (compute the pressure step\n"
          "                   "
          "   using the continuity equation)")};

  const char *ipucou_value_str[]
    = {N_("0 (standard algorithm for velocity/pressure coupling)\n"),
       N_("1 (reinforced velocity/pressure coupling\n"
          "                   "
          "   in case calculation with long time steps)")};

  cs_log_printf(CS_LOG_SETUP,
                _("    iprco:         %s\n"),
                _(iprco_value_str[vp_param->iprco]));

  cs_log_printf(CS_LOG_SETUP,
                _("    ipucou:        %s\n"),
                _(ipucou_value_str[vp_param->ipucou]));

  cs_log_printf
    (CS_LOG_SETUP,
     _("    irevmc:     %5d (Velocity reconstruction mode)\n"),
     vp_param->irevmc);

  const char *itpcol_type_str[] = {N_("staggered time scheme"),
                                   N_("colocated time scheme")};

  cs_log_printf(CS_LOG_SETUP,
                _("    itpcol:        %d (%s)\n"),
                vp_param->itpcol,
                _(itpcol_type_str[vp_param->itpcol]));

  const cs_equation_param_t *eqp = NULL;

  if (cs_glob_time_step_options->idtvar >= 0) {
    eqp = cs_field_get_equation_param_const(f_p);
    cs_log_printf
      (CS_LOG_SETUP,
       _("    relaxv:      %14.5e for %s (relaxation)\n"
         "    arak:        %14.5e (Arakawa factor)\n"),
       eqp->relaxv, f_p_label, vp_param->arak);
  }
  else {
    eqp = cs_field_get_equation_param_const(CS_F_(vel));
    cs_log_printf
      (CS_LOG_SETUP,
       _("    arak:        %14.5e (Arakawa factor)\n"),
       eqp->relaxv * vp_param->arak);
  }
  cs_log_printf
      (CS_LOG_SETUP,
       _("\n  Factor of Rhie and Chow %d\n"),
       vp_param->rcfact);

  cs_log_printf
    (CS_LOG_SETUP,
     _("    staggered %d (1D staggered scheme option)\n"),
     vp_param->staggered);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

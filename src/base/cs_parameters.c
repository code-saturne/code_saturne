/*============================================================================
 * General parameters management.
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
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_at_opt_interp.h"
#include "cs_convection_diffusion.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_gradient.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_rad_transfer_fields.h"
#include "cs_turbulence_model.h"
#include "cs_time_moment.h"
#include "cs_thermal_model.h"
#include "cs_cf_model.h"
#include "cs_tree.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_parameters.c
        General parameters and options management.
*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_space_disc_t

  \brief Space discretisation options descriptor.

  Members of the space discretisation structure are publicly accessible, to
  allow for concise syntax, as they are expected to be used in many places.

  \var  cs_space_disc_t::imvisf
        <a name="imvisf"></a>
        face viscosity field interpolation
        - 1: harmonic
        - 0: arithmetic (default)
  \var  cs_space_disc_t::imrgra
        type of gradient reconstruction
        - 0: Green-Gauss with iterative handling of non-orthogonalities
        - 1: standard least squares method
        - 2: least squares method with extended neighborhood
        - 3: least squares method with reduced extended neighborhood
        - 4: Green-Gauss with least squares gradient face values
        - 5: Green-Gauss with least squares gradient over extended
             neighborhood face values
        - 6: Green-Gauss with least squares gradient over reduced
             extended neighborhood face values
        - 7: Green-Gauss with face values based on cell to vertex interpolation

  \var  cs_space_disc_t::iflxmw
        method to compute interior mass flux due to ALE mesh velocity
        - 0: based on nodes displacement
        - 1: based on cell center mesh velocity

  \var  cs_space_disc_t::itbrrb
        accurate treatment of the wall temperature
        (reconstruction of wall temperature)
        - 1: true
        - 0: false (default)
        (see \ref cs_boundary_condition_set_coeffs,
        useful in case of coupling with syrthes)
*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_time_scheme_t

  \brief Time scheme descriptor.

  Members of the time scheme structure are publicly accessible, to allow for
  concise  syntax, as they are expected to be used in many places.

  \var  cs_time_scheme_t::time_order
        Global time order of time stepping
        - 2: 2nd order
        - 1: 1st order (default)

  \var  cs_time_scheme_t::istmpf
        \anchor istmpf
        Time order of the mass flux scheme
        The chosen value for \ref istmpf will automatically
        determine the value given to the variable \ref thetfl.
        - 2: theta scheme with theta > 0 (theta=0.5 means 2nd order)
             the mass flow used in the momentum equations is extrapolated at
             n+ \ref thetfl (= n+1/2) from the values at the two former time
             steps (Adams Bashforth); the mass flow used in the equations for
             turbulence and scalars is interpolated at time n+ \ref thetfl
             (= n+1/2) from the values at the former time step and at the
             newly calculated \f$n+1\f$ time step.
        - 0: theta scheme with theta = 0 (explicit): the mass flow
             calculated at the previous time step is used in the convective
             terms of all the equations (momentum, turbulence and scalars)
        - 1: implicit scheme (default) : the mass flow calculated
             at the previous time step is used in the convective terms of the
             momentum equation, and the updated mass flow is used in the
             equations of turbulence and scalars. By default, \ref istmpf=2
             is used in the case of a second-order time scheme (if \ref ischtp=2)
             and \ref istmpf = 1 otherwise.

  \var  cs_time_scheme_t::isno2t
        \anchor isno2t
        Specifies the time scheme for the source
        terms of the momentum equation, apart from convection and
        diffusion (for instance: head loss, transposed gradient, ...).
        - 0: "standard" first-order: the terms which are linear
             functions of the solved variable are implicit and the others
             are explicit.
        - 1: second-order: the terms of the form \f$S_i\phi\f$ which are
             linear functions of the solved variable \f$\phi\f$ are expressed
             as second-order terms by interpolation (according to the formula
             \f$(S_i\phi)^{n+\theta}=S_i^n[(1-\theta)\phi^n+\theta\phi^{n+1}]\f$,
             \f$\theta\f$ being given by the value of \ref thetav associated
             with the variable \f$\phi\f$); the other terms \f$S_e\f$ are
             expressed as second-order terms by extrapolation (according to the
             formula \f$(S_e)^{n+\theta}=[(1+\theta)S_e^n-\theta S_e^{n-1}]\f$,
             \f$\theta\f$ being given by the value of \ref thetsn = 0.5).\n
        - 2: the linear terms \f$S_i\phi\f$ are treated in the same
             way as when \ref isno2t = 1; the other terms \f$S_e\f$ are
             extrapolated according to the same formula as when \ref isno2t = 1,
             but with \f$\theta\f$= \ref thetsn = 1. By default, \ref isno2t
             is initialized to 1 (second-order) when the selected time scheme
             is second-order (\ref ischtp = 2), otherwise to 0.

  \var  cs_time_scheme_t::isto2t
        \anchor isto2t
        Specifies the time scheme for
        the source terms of the turbulence equations i.e. related to
        \f$k\f$, \f$R_{ij}\f$, \f$\varepsilon\f$, \f$\omega\f$, \f$\varphi\f$,
        \f$\overline{f}\f$), apart from convection and diffusion.
        - 0: standard first-order: the terms which are linear
             functions of the solved variable are implicit and the others
             are explicit.
        - 1: second-order: the terms of the form \f$S_i\phi\f$ which are
             linear functions of the solved variable \f$\phi\f$ are
             expressed as second-order terms by interpolation (according to
             the formula
             \f$(S_i\phi)^{n+\theta}=S_i^n[(1-\theta)\phi^n+\theta\phi^{n+1}]\f$,
             \f$\theta\f$ being given by the value of \ref thetav associated
             with the variable \f$\phi\f$); the other terms \f$S_e\f$ are
             expressed as second-order terms by extrapolation (according to
             the formula
             \f$(S_e)^{n+\theta}=[(1+\theta)S_e^n-\theta S_e^{n-1}]\f$,
             \f$\theta\f$ being given by the value of \ref thetst = 0.5)
        - 2: the linear terms \f$S_i\phi\f$ are treated in the same
             \ref isto2t = 1; the other terms \f$S_e\f$ are way as when
             extrapolated according to the same formula as when
             \ref isto2t = 1, but with \f$\theta\f$= \ref thetst = 1.\n
             Due to certain specific couplings between the turbulence equations,
             \ref isto2t is allowed the value 1 or 2 only for the
             \f$R_{ij}\f$ models (\ref iturb = 30 or 31);
             hence, it is always initialised to 0.

  \var  cs_time_scheme_t::thetsn
        \anchor thetsn
        \f$ \theta_S \f$-scheme for the source terms \f$S_e\f$ in the
        Navier-Stokes equations when the source term extrapolation has
        been activated (see \ref isno2t), following the formula
        \f$(S_e)^{n+\theta}=(1+\theta)S_e^n-\theta S_e^{n-1}\f$.\n The value
        \f$theta\f$ = \ref thetsn is deduced from the value chosen for
        \ref isno2t. Generally only the value 0.5 is used.
        -  0 : second viscosity explicit
        - 1/2: second viscosity extrapolated in n+1/2
        -  1 : second viscosity extrapolated in n+1

  \var  cs_time_scheme_t::thetst
        \anchor thetst
        \f$ \theta \f$-scheme for the extrapolation of the nonlinear
        explicit source terms $S_e$ of the turbulence equations when the
        source term extrapolation has been activated (see \ref isto2t),
        following the formula
        \f$(S_e)^{n+\theta}=(1+\theta)S_e^n-\theta S_e^{n-1}\f$.\n
        The value of \f$theta\f$ is deduced from the value chosen for
        \ref isto2t. Generally, only the value 0.5 is used.
        -  0 : explicit
        - 1/2: extrapolated in n+1/2
        -  1 : extrapolated in n+1

  \var  cs_time_scheme_t::thetvi
        \anchor thetvi
        \f$ \theta \f$-scheme for the extrapolation of the physical
        property \f$\phi\f$ "total viscosity" when the extrapolation
        has been activated (see \ref time_extrapolated key word), according to
        the formula \f$\phi^{n+\theta}=(1+\theta)\phi^n-\theta \phi^{n-1}\f$.\n
        The value of \f$\theta\f$ = \ref thetvi is deduced from the value
        chosen for \ref time_extrapolated key word for the viscosity.
        Generally, only the value 0.5 is used.
        -  0 : explicit
        - 1/2: extrapolated in n+1/2
        -  1 : extrapolated in n+1

  \var  cs_time_scheme_t::thetcp
        \anchor thetcp
        \f$ \theta \f$-scheme for the extrapolation of the physical
        property \f$\phi\f$ "specific heat" when the extrapolation has
        been activated (see \ref time_extrapolated field key int), according to
        the formula \f$\phi^{n+\theta}=(1+\theta)\phi^n-\theta \phi^{n-1}\f$.\n
        The value of \f$\theta\f$ = \ref thetcp is deduced from the value chosen
        for the specific heat. Generally, only the value 0.5 is used.
        -  0 : explicit
        - 1/2: extrapolated in n+1/2
        -  1 : extrapolated in n+1

  \var  cs_time_scheme_t::iccvfg
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

*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_restart_auxiliary_t

  \brief Additional checkpoint/restart files

  \var  cs_restart_auxiliary_t::read_auxiliary
        activate reading of auxiliary restart file
  \var  cs_restart_auxiliary_t::write_auxiliary
        activate writing of auxiliary restart file
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _CS_HODGE_LEGACY_INIT \
{.inv_pty = false, \
.type = CS_HODGE_N_TYPES, \
.algo = CS_HODGE_N_ALGOS, \
.coef = 0}

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure of variable calculation options mappable to Fortran
 *----------------------------------------------------------------------------*/

typedef struct {

  int     iwarni;
  int     iconv;
  int     istat;
  int     idircl;
  int     ndircl;
  int     idiff;
  int     idifft;
  int     idften;
  int     iswdyn;
  int     ischcv;
  int     ibdtso;
  int     isstpc;
  int     nswrgr;
  int     nswrsm;
  int     imvisf;
  int     imrgra;
  int     imligr;
  int     ircflu;
  int     iwgrec;       /* gradient calculation
                           - 0: standard (default)
                           - 1: weighted (could be used with imvisf = 1) */
  int     icoupl;       /* internal coupling
                           - -1: not coupled (default)
                           -  1: coupled                                 */

  double  thetav;
  double  blencv;
  double  blend_st;
  double  epsilo;
  double  epsrsm;
  double  epsrgr;
  double  climgr;
  double  relaxv;

} cs_f_var_cal_opt_t;

/*----------------------------------------------------------------------------*/

/* Definition of user variable */

typedef struct {

  char     *name;               /* Variable name */
  char     *ref_name;           /* Name of variable referred to */
  int       dim;                /* Variable dimension */
  bool      is_variance;        /* True if the variable is a variance */

} cs_user_variable_def_t;

/* Definition of user property */

typedef struct {

  char     *name;               /* Property name */
  int       dim;                /* Property dimension */
  int       location_id;        /* Property location id */

} cs_user_property_def_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default variable compute options */

static cs_equation_param_t _equation_param_default
= {
   .name = NULL,
   .type = CS_EQUATION_N_TYPES,
   .dim = 1,
   .verbosity = 0,

   .flag = 0,
   .post_flag = 0,
   .space_scheme = CS_SPACE_SCHEME_LEGACY,
   .dof_reduction = CS_PARAM_REDUCTION_AVERAGE,
   .space_poly_degree = 0,

   .iconv  = 1,
   .istat  = 1,
   .idircl = 1,
   .ndircl = 0,
   .idiff  = 1,
   .idifft = 1,
   .idften = CS_ISOTROPIC_DIFFUSION,
   .iswdyn = -1,
   .ischcv = 1,
   .ibdtso = 1,
   .isstpc = 1,
   .nswrgr = 100,
   .nswrsm = 1,
   .imvisf = 0,
   .imrgra = -1,
   .imligr = CS_GRADIENT_LIMIT_NONE,
   .ircflu = 1,
   .iwgrec = 0,
   .icoupl = -1,
   .thetav = 1.,
   .blencv = 1.,
   .blend_st = 0.,
   .epsilo = 1.e-5,
   .epsrsm = 1.e-4,
   .epsrgr = 1.e-4,
   .climgr = 1.5,
   .relaxv = 1.,
   .b_gradient_r = 2,

   .default_bc = CS_PARAM_BC_HMG_NEUMANN,
   .n_bc_defs = 0,
   .bc_defs = NULL,

   .default_enforcement = CS_PARAM_BC_ENFORCE_ALGEBRAIC,
   .strong_pena_bc_coeff = -1,
   .weak_pena_bc_coeff = -1,

   .n_ic_defs = 0,
   .ic_defs = NULL,

   .do_lumping = false,
   .time_hodgep = _CS_HODGE_LEGACY_INIT,
   .time_property = NULL,
   .time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT,
   .theta = 1,

   .diffusion_hodgep = _CS_HODGE_LEGACY_INIT,
   .diffusion_property = NULL,
   .curlcurl_hodgep = _CS_HODGE_LEGACY_INIT,
   .curlcurl_property = NULL,
   .graddiv_hodgep = _CS_HODGE_LEGACY_INIT,
   .graddiv_property = NULL,

   .adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV,
   .adv_scheme = CS_PARAM_N_ADVECTION_SCHEMES,
   .upwind_portion = 0.,
   .adv_field = NULL,
   .adv_scaling_property = NULL,

   .reaction_hodgep = _CS_HODGE_LEGACY_INIT,
   .n_reaction_terms = 0,
   .reaction_properties = NULL,
   .n_source_terms = 0,
   .source_terms = NULL,
   .n_volume_mass_injections = 0,
   .volume_mass_injections = NULL,

   .n_enforcements = 0,
   .enforcement_params = NULL,

   .sles_param = NULL,
  };

/* Space discretisation options structure and associated pointer */

static cs_space_disc_t  _space_disc =
{
  .imvisf = 0,
  .imrgra = 4,
  .iflxmw = 0,
  .itbrrb = 0
};

const cs_space_disc_t  *cs_glob_space_disc = &_space_disc;

/* Time scheme options structure and associated pointer */

static cs_time_scheme_t  _time_scheme =
{
  .time_order = -1,
  .istmpf = -999,
  .isno2t = -999,
  .isto2t = -999,
  .thetsn = -999.0,
  .thetst = -999.0,
  .thetvi = -999.0,
  .thetcp = -999.0,
  .iccvfg = 0
};

const cs_time_scheme_t  *cs_glob_time_scheme = &_time_scheme;

/* Auxiliary checkpoint/restart file parameters */

static cs_restart_auxiliary_t  _restart_auxiliary =
{
  .read_auxiliary = 1,
  .write_auxiliary = 1
};

cs_restart_auxiliary_t  *cs_glob_restart_auxiliary = &_restart_auxiliary;

/* Definition of user variables and properties */

int                      _n_user_variables = 0;
int                      _n_user_properties = 0;
cs_user_variable_def_t  *_user_variable_defs = NULL;
cs_user_property_def_t  *_user_property_defs = NULL;

static cs_solving_info_t _solving_info =
{
  0,     /* n_it: number of iterations for the linear solver */
  0.,    /* rhs_norm: right hand side norm                   */
  0.,    /* res_norm: normed residual                        */
  0.,    /* derive: norm of the time derivative              */
  0.,    /* l2residual: L2 time residual                     */
};

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Global parameters tree structure */

cs_tree_node_t  *cs_glob_tree = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_space_disc_get_pointers(int     **imvisf,
                             int     **imrgra,
                             int     **iflxmw,
                             int     **itbrrb);

void
cs_f_time_scheme_get_pointers(int     **ischtp,
                              int     **istmpf,
                              int     **isno2t,
                              int     **isto2t,
                              double  **thetsn,
                              double  **thetst,
                              double  **thetvi,
                              double  **thetcp,
                              int     **iccvfg);

void
cs_f_restart_auxiliary_get_pointers(int  **ileaux,
                                    int  **iecaux);

void
cs_f_field_get_key_struct_var_cal_opt(int                  f_id,
                                      cs_f_var_cal_opt_t  *vcopt);

void
cs_f_field_set_key_struct_var_cal_opt(int                        f_id,
                                      const cs_f_var_cal_opt_t  *vcopt);

void *
cs_f_equation_param_from_var_cal_opt(const cs_f_var_cal_opt_t  *vcopt);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log values of the structure */

static void
_log_func_var_cal_opt(const void *t)
{
  const char fmt_i[] = N_("      %-19s  %d\n");
  const char fmt_r[] = N_("      %-19s  %-12.3g\n");
  const cs_var_cal_opt_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, fmt_i, "verbosity", _t->verbosity);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iconv ", _t->iconv );
  cs_log_printf(CS_LOG_SETUP, fmt_i, "istat ", _t->istat );
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idircl", _t->idircl);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ndircl", _t->ndircl);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idiff ", _t->idiff );
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idifft", _t->idifft);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idften", _t->idften);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iswdyn", _t->iswdyn);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ischcv", _t->ischcv);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ibdtso", _t->ibdtso);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "isstpc", _t->isstpc);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "nswrgr", _t->nswrgr);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "nswrsm", _t->nswrsm);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "imvisf", _t->imvisf);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "imrgra", _t->imrgra);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "imligr", _t->imligr);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ircflu", _t->ircflu);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iwgrec", _t->iwgrec);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "icoupl", _t->icoupl);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "thetav", _t->thetav);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "blencv", _t->blencv);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "blend_st", _t->blend_st);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "epsilo", _t->epsilo);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "epsrsm", _t->epsrsm);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "epsrgr", _t->epsrgr);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "climgr", _t->climgr);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "relaxv", _t->relaxv);

  cs_log_printf(CS_LOG_SETUP, fmt_i, "b_gradient_r", _t->b_gradient_r);
}

/* Log default values of the structure */

static void
_log_func_default_var_cal_opt(const void *t)
{
  const char fmt_i[] = "      %-19s  %-12d %s\n";
  const char fmt_r[] = "      %-19s  %-12.3g %s\n";
  const char fmt_c[] = "      %-19s  %-12s %s\n";
  const cs_equation_param_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP,"  var_cal_opt\n");

  cs_log_printf(CS_LOG_SETUP,_("    Printing\n"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "verbosity", _t->verbosity,
                _("Verbosity level."));

  cs_log_printf(CS_LOG_SETUP,"    Time stepping\n");
  cs_log_printf(CS_LOG_SETUP, fmt_i, "istat ", _t->istat,
                _("Take unsteady terms into account."));

  cs_log_printf(CS_LOG_SETUP,"    Convection/Diffusion\n");

  cs_log_printf(CS_LOG_SETUP, fmt_i, "iconv ", _t->iconv,
                _("Take convection into account."));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idiff ", _t->idiff,
                _("Take diffusion into account."));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idifft", _t->idifft,
                _("Take turbulent diffusion into account."));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idften", _t->idften,
                _("Type of diffusivity: scalar (1), orthotropic (3) "
                  "or symmetric tensor (6)"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ischcv", _t->ischcv,
                _("Type of convective scheme:"));
  cs_log_printf(CS_LOG_SETUP, fmt_c, " ", " ",
                _("  0: 2nd order with centered-gradient upwind reconstruction,"));
  cs_log_printf(CS_LOG_SETUP, fmt_c, " ", " ",
                _("  1: centered,"));
  cs_log_printf(CS_LOG_SETUP, fmt_c, " ", " ",
                _("  2: 2nd order with upwind-gradient upwind-reconstruction "
                  "(SOLU)"));
  cs_log_printf(CS_LOG_SETUP, fmt_c, " ", " ",
                _("  3: continuous blending between upwind and another scheme"));
  cs_log_printf(CS_LOG_SETUP, fmt_c, " ", " ",
                _("  4: NVD/TVD scheme"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "isstpc", _t->isstpc,
                _("0 for slope test, 1 for no slope test, 2 for min/max limiter "));
  cs_log_printf(CS_LOG_SETUP, fmt_r, "blencv", _t->blencv,
                _("[0.;1.] (1-upwind proportion (0: upwind))"));
  cs_log_printf(CS_LOG_SETUP, fmt_r, "blend_st", _t->blend_st,
                _("[0.;1.] (1-upwind proportion after slope test (0: upwind))"));

  cs_log_printf(CS_LOG_SETUP,"    Gradients calculation\n");
  cs_log_printf(CS_LOG_SETUP, fmt_i, "imrgra", _t->imrgra,
                _("Reconstruction mode"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "nswrgr", _t->nswrgr,
                _("Number of sweeps gradient reconstruction"));
  cs_log_printf(CS_LOG_SETUP, fmt_r, "epsrgr", _t->epsrgr,
                _("Gradient reconstruction precision"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "imligr", _t->imligr,
                _("< 0, 0 or 1 (gradient limitation method)"));
  cs_log_printf(CS_LOG_SETUP, fmt_r, "climgr", _t->climgr,
                _("> 1 or 1 (gradient limitation coefficient)"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iwgrec", _t->iwgrec,
                _("Gradient calculation: standard (0) or weighted (1)"));

  cs_log_printf(CS_LOG_SETUP,"    Rhs reconstruction\n");
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ircflu", _t->ircflu,
                _("0 or 1 (flux reconstruction)"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "nswrsm", _t->nswrsm,
                _("Number of sweeps rhs reconstruction"));
  cs_log_printf(CS_LOG_SETUP, fmt_r, "epsrsm", _t->epsrsm,
                _("Rhs reconstruction precision"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iswdyn", _t->iswdyn,
                _("Dynamic relaxation type"));

  cs_log_printf(CS_LOG_SETUP,"    Iterative solvers\n");
  cs_log_printf(CS_LOG_SETUP, fmt_r, "epsilo", _t->epsilo,
                _("Resolution precision"));

  cs_log_printf(CS_LOG_SETUP,"    Time-scheme\n");
  cs_log_printf(CS_LOG_SETUP, fmt_r, "thetav", _t->thetav,
                _("[0.;1.] theta-scheme for the main variables (0.5 for "
                  "Crank-Nicolson)"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ibdtso", _t->ibdtso,
                _("Backward differential scheme in time order"));
  cs_log_printf(CS_LOG_SETUP, fmt_r, "relaxv", _t->relaxv,
                _("Relaxation of variables (1 for no relaxation)"));
}

/*----------------------------------------------------------------------------
 * Copy values from a Fortran var_cal_opt structure to an
 * equation_params_t structure, whose other members are unchanged.
 *
 * parameters:
 *   vcopt  <-- associated Fortran var_cal_opt structure
 *   eqp    <-> associated cs_equation_params_t structure
 *----------------------------------------------------------------------------*/

static void
_var_cal_opt_to_equation_params(const cs_f_var_cal_opt_t  *vcopt,
                                cs_equation_param_t       *eqp)
{
  eqp->verbosity = vcopt->iwarni;
  eqp->iconv  = vcopt->iconv;
  eqp->istat  = vcopt->istat;
  eqp->idircl = vcopt->idircl;
  eqp->ndircl = vcopt->ndircl;
  eqp->idiff  = vcopt->idiff;
  eqp->idifft = vcopt->idifft;
  eqp->idften = vcopt->idften;
  eqp->iswdyn = vcopt->iswdyn;
  eqp->ischcv = vcopt->ischcv;
  eqp->ibdtso = vcopt->ibdtso;
  eqp->isstpc = vcopt->isstpc;
  eqp->nswrgr = vcopt->nswrgr;
  eqp->nswrsm = vcopt->nswrsm;
  eqp->imvisf = vcopt->imvisf;
  eqp->imrgra = vcopt->imrgra;
  eqp->imligr = vcopt->imligr;
  eqp->ircflu = vcopt->ircflu;
  eqp->iwgrec = vcopt->iwgrec;
  eqp->icoupl = vcopt->icoupl;

  eqp->thetav = vcopt->thetav;
  eqp->blencv = vcopt->blencv;
  eqp->blend_st = vcopt->blend_st;
  eqp->epsilo = vcopt->epsilo;
  eqp->epsrsm = vcopt->epsrsm;
  eqp->epsrgr = vcopt->epsrgr;
  eqp->climgr = vcopt->climgr;
  eqp->relaxv = vcopt->relaxv;
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global space disc structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_space_disc_get_pointers(int     **imvisf,
                             int     **imrgra,
                             int     **iflxmw,
                             int     **itbrrb)
{
  *imvisf = &(_space_disc.imvisf);
  *imrgra = &(_space_disc.imrgra);
  *iflxmw = &(_space_disc.iflxmw);
  *itbrrb = &(_space_disc.itbrrb);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global time scheme structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *----------------------------------------------------------------------------*/

void
cs_f_time_scheme_get_pointers(int     **ischtp,
                              int     **istmpf,
                              int     **isno2t,
                              int     **isto2t,
                              double  **thetsn,
                              double  **thetst,
                              double  **thetvi,
                              double  **thetcp,
                              int     **iccvfg)
{
  *ischtp = &(_time_scheme.time_order);
  *istmpf = &(_time_scheme.istmpf);
  *isno2t = &(_time_scheme.isno2t);
  *isto2t = &(_time_scheme.isto2t);
  *thetsn = &(_time_scheme.thetsn);
  *thetst = &(_time_scheme.thetst);
  *thetvi = &(_time_scheme.thetvi);
  *thetcp = &(_time_scheme.thetcp);
  *iccvfg = &(_time_scheme.iccvfg);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global restart_auxiliary structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ileaux  --> pointer to cs_glob_restart_auxiliary->read_auxiliary
 *   iecaux  --> pointer to cs_glob_restart_auxiliary->write_auxiliary
 *----------------------------------------------------------------------------*/

void
cs_f_restart_auxiliary_get_pointers(int  **ileaux,
                                    int  **iecaux)
{
  *ileaux = &(_restart_auxiliary.read_auxiliary);
  *iecaux = &(_restart_auxiliary.write_auxiliary);
}

/*----------------------------------------------------------------------------
 * Copy values from cs_equation_params structure of a given field
 * to the matching Fortran var_cal_opt structure
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   f_id   associated field id
 *   vcopt  associated structure
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_key_struct_var_cal_opt(int                  f_id,
                                      cs_f_var_cal_opt_t  *vcopt)
{
  static int c_k_id = -1;
  if (c_k_id < 0)
    c_k_id = cs_field_key_id("var_cal_opt");

  const cs_equation_param_t *eqp
    = cs_field_get_key_struct_const_ptr(cs_field_by_id(f_id),
                                        c_k_id);

  vcopt->iwarni = eqp->verbosity;
  vcopt->iconv  = eqp->iconv;
  vcopt->istat  = eqp->istat;
  vcopt->idircl = eqp->idircl;
  vcopt->ndircl = eqp->ndircl;
  vcopt->idiff  = eqp->idiff;
  vcopt->idifft = eqp->idifft;
  vcopt->idften = eqp->idften;
  vcopt->iswdyn = eqp->iswdyn;
  vcopt->ischcv = eqp->ischcv;
  vcopt->ibdtso = eqp->ibdtso;
  vcopt->isstpc = eqp->isstpc;
  vcopt->nswrgr = eqp->nswrgr;
  vcopt->nswrsm = eqp->nswrsm;
  vcopt->imvisf = eqp->imvisf;
  vcopt->imrgra = eqp->imrgra;
  vcopt->imligr = eqp->imligr;
  vcopt->ircflu = eqp->ircflu;
  vcopt->iwgrec = eqp->iwgrec;
  vcopt->icoupl = eqp->icoupl;

  vcopt->thetav = eqp->thetav;
  vcopt->blencv = eqp->blencv;
  vcopt->blend_st = eqp->blend_st;
  vcopt->epsilo = eqp->epsilo;
  vcopt->epsrsm = eqp->epsrsm;
  vcopt->epsrgr = eqp->epsrgr;
  vcopt->climgr = eqp->climgr;
  vcopt->relaxv = eqp->relaxv;
}

/*----------------------------------------------------------------------------
 * Copy values from cs_equation_params structure of a given field
 * to the matching Fortran var_cal_opt structure
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   f_id   associated field id
 *   vcopt  associated structure
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_key_struct_var_cal_opt(int                        f_id,
                                      const cs_f_var_cal_opt_t  *vcopt)
{
  static int c_k_id = -1;
  if (c_k_id < 0)
    c_k_id = cs_field_key_id("var_cal_opt");

  cs_equation_param_t *eqp
    = cs_field_get_key_struct_ptr(cs_field_by_id(f_id), c_k_id);

  _var_cal_opt_to_equation_params(vcopt, eqp);
}

/*----------------------------------------------------------------------------
 * Return a pointer to a cs_equation_params structure initialized from
 * a Fortran var_cal_opt structure.
 *
 * Note that to avoid issues with the structure not being interoperable
 * with Fortran, and its size not being known easily in Fortran either
 * (and subject to change with code maintenance), a pointer to a static
 * variable of this function is returned. This means this function is not
 * thread safe, which should not be an issue, since the functions to which
 * this object is passed are not expected to be either, as they are quite
 * high level and usually include MPI operations on a global communicator.
 *
 * parameters:
 *   vcopt <-- Fortran var_cal_opt
 *
 * returns:
 *   pointer to matching cs_equation_params
 *----------------------------------------------------------------------------*/

void *
cs_f_equation_param_from_var_cal_opt(const cs_f_var_cal_opt_t  *vcopt)
{
  static cs_equation_param_t eqp;
  memcpy(&eqp, &_equation_param_default, sizeof(cs_equation_param_t));

  _var_cal_opt_to_equation_params(vcopt, &eqp);

  return &eqp;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cs_glob_space_disc.
 *
 * Needed to initialize structure with GUI and user C functions.
 *
 * \return  velocity_pressure information structure
 */
/*----------------------------------------------------------------------------*/

cs_space_disc_t *
cs_get_glob_space_disc(void)
{
  return &_space_disc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to cs_glob_time_scheme
 *
 * needed to initialize structure with GUI and user C functions.
 *
 * \return  time scheme information structure
 */
/*----------------------------------------------------------------------------*/

cs_time_scheme_t *
cs_get_glob_time_scheme(void)
{
  return &_time_scheme;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define general field keys.
 *
 * A recommended practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_keys(void)
{
  /* field ids of the mass or volume fluxes convecting the variable field */
  cs_field_define_key_int("inner_mass_flux_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("boundary_mass_flux_id", -1, CS_FIELD_VARIABLE);

  /* field ids of the fluxes of the variable field (it needs to be stored
     for some quantities such as void fraction flux for the VoF algo.) */
  cs_field_define_key_int("inner_flux_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("boundary_flux_id", -1, CS_FIELD_VARIABLE);

  /* field id of the variable associated to a given property */
  cs_field_define_key_int("parent_field_id", -1, 0);

  cs_field_define_key_int("variable_id", -1, 0); /* inverse of ivarfl(ivar) */
  cs_field_define_key_int("scalar_id", -1, 0);   /* inverse of isca(iscal) */

  cs_field_define_key_int("diffusion_coef_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_double("diffusion_coef_ref",
                             -1.e12*10., CS_FIELD_VARIABLE);
  cs_field_define_key_int("diffusivity_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_double("diffusivity_ref",
                             -1.e12*10., CS_FIELD_VARIABLE);

  cs_field_define_key_int("scalar_diffusivity_prev", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("turbulent_diffusivity_id", -1, CS_FIELD_VARIABLE);

  /* Only used for turbulent scalar flux model in LES */
  cs_field_define_key_int("sgs_scalar_flux_coef_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("density_id", -1, CS_FIELD_VARIABLE);

  /* Is the field buoyant? 0 if not, 1 if yes */
  cs_field_define_key_int("is_buoyant", 0, CS_FIELD_VARIABLE);

  /* Does the field behave like a temperature ? (iscacp) */
  cs_field_define_key_int("is_temperature", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("turbulent_flux_model", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("turbulent_flux_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("alpha_turbulent_flux_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_double("turbulent_flux_ctheta",
                             1., CS_FIELD_VARIABLE); /* ctheta(iscal) */

  cs_field_define_key_double("time_step_factor",
                             1., CS_FIELD_VARIABLE); /* cdtvar(ivar) */

  cs_field_define_key_double("turbulent_schmidt",
                             1., CS_FIELD_VARIABLE); /* sigmas(iscal) */
  cs_field_define_key_int("turbulent_schmidt_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("gradient_weighting_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("diffusivity_tensor", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("drift_scalar_model", 0, 0);

  cs_field_define_key_int("scalar_class", 0, 0);
  cs_field_define_key_int("first_moment_id", -1, 0); /* iscavr(iscal) */

  cs_field_define_key_int("variance_clipping", -1, 0); /* iclvfl(iscal) */
  cs_field_define_key_double("variance_dissipation", 0.8, 0); /* rvarfl(iscal) */

  cs_field_define_key_int("syrthes_coupling", 0, 0); /* icpsyr(iscal) */

  cs_field_define_key_int("source_term_prev_id", -1, CS_FIELD_VARIABLE);
  /* TODO merge with previous key word */
  cs_field_define_key_int("source_term_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("slope_test_upwind_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("clipping_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("is_clipped", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("boundary_value_id", -1, 0);

  cs_field_define_key_int("convection_limiter_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("diffusion_limiter_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("coupling_entity", -1, 0);

  /*
   * Is the field time-extrapolated?
   * -1: default automatic value
   *  0: "standard" first-order: the value calculated at
   *     the beginning of the current time step (from the
   *     variables known at the end of the previous time step) is used
   *  1: second-order: the physical property \f$\phi\f$ is
   *     extrapolated according to the formula
   *     \f$\phi^{n+\theta}=[(1+\theta)\phi^n-\theta \phi^{n-1}]\f$,
   *     \f$\theta\f$ being given by the value of 0.5
   *  2: first-order: the physical property \f$\phi\f$ is
   *     extrapolated at $n+1$ according to the same formula
   *     as when = 1 but with \f$\theta\f$ = 1
   */

  cs_field_define_key_int("time_extrapolated", -1, 0);
  cs_field_define_key_int("scalar_time_scheme", -1, 0); /* ex-isso2t */
  cs_field_define_key_double("st_exp_extrapolated", -1,
                             CS_FIELD_VARIABLE); /* ex-thetss */
  cs_field_define_key_double("diffusivity_extrapolated", -1,
                             CS_FIELD_VARIABLE); /* ex-thetvs */
  cs_field_define_key_int("measures_set_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("opt_interp_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("opt_interp_analysis_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_double("min_scalar_clipping", -1.e12, 0);
  cs_field_define_key_double("max_scalar_clipping", 1.e12, 0);

  /* Bounds of a given scalar which won't be used in clipping */

  cs_field_define_key_double("max_scalar", 1., 0);
  cs_field_define_key_double("min_scalar", 0., 0);

 /* Integer corresponding to the type of Roe-Sweby Limiter:
  * 1->minmod
  * 2->van-leer
  * 3->van-albada
  * 4->superbee */

  cs_field_define_key_int("limiter_choice", -1, CS_FIELD_VARIABLE);

  /* Structure containing the calculation options of the field variables */
  cs_field_define_key_struct
    ("var_cal_opt",
     &_equation_param_default,
     _log_func_var_cal_opt,
     _log_func_default_var_cal_opt,
     (cs_field_clear_key_struct_t *)cs_equation_clear_param,
     sizeof(cs_var_cal_opt_t),
     CS_FIELD_VARIABLE);

  /* Structure containing the solving info of the field variables
     (used for log, not setup, so set NULL setup logging function) */
  cs_field_define_key_struct("solving_info",
                             &_solving_info,
                             NULL,
                             NULL,
                             NULL,
                             sizeof(cs_solving_info_t),
                             CS_FIELD_VARIABLE);
  cs_field_key_disable_setup_log(cs_field_key_id("solving_info"));

  /* Restart options */
  cs_field_define_key_int("restart_file", CS_RESTART_DISABLED, 0);
  cs_field_define_key_int("restart_n_values", 1, 0);

  /* field units */
  cs_field_define_key_str("units", "", 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read general restart info.
 *
 * This updates the previous time step info and notebook varaibles values.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_read_restart_info(void)
{
  if (cs_restart_present()) {
    cs_restart_t *r
      = cs_restart_create("main.csc", "restart", CS_RESTART_MODE_READ);
    cs_restart_read_time_step_info(r);
    cs_restart_read_notebook_variables(r);
    cs_restart_destroy(&r);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user variable.
 *
 * Solved variables are always defined on cells.
 *
 * \param[in]  name  name of variable and associated field
 * \param[in]  dim   variable dimension
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_variable(const char  *name,
                           int          dim)
{
  BFT_REALLOC(_user_variable_defs,
              _n_user_variables + 1,
              cs_user_variable_def_t);

  BFT_MALLOC((_user_variable_defs + _n_user_variables)->name,
             strlen(name) + 1,
             char);
  strcpy((_user_variable_defs + _n_user_variables)->name, name);

  (_user_variable_defs + _n_user_variables)->dim = dim;
  (_user_variable_defs + _n_user_variables)->is_variance = false;

  if (dim > 3)
    bft_error(__FILE__, __LINE__, 0,
              _("Only user variables of dimension lower or equal to 3 are"
                "currently handled,\nbut %s is defined with dimension %d."),
              name, dim);

  _n_user_variables++;

  /* Make this immediate if fields have already been defined */

  if (cs_field_n_fields() > 0)
    cs_parameters_create_added_variables();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user variable which is a variance of another variable.
 *
 * Only variances of thermal or user-defined variables are currently handled.
 *
 * \param[in]  name           name of variance and associated field
 * \param[in]  variable_name  name of associated variable
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_variable_variance(const char  *name,
                                    const char  *variable_name)
{
  BFT_REALLOC(_user_variable_defs,
              _n_user_variables + 1,
              cs_user_variable_def_t);
  BFT_MALLOC((_user_variable_defs + _n_user_variables)->name,
             strlen(name) + 1,
             char);
  BFT_MALLOC((_user_variable_defs + _n_user_variables)->ref_name,
             strlen(variable_name) + 1,
             char);

  strcpy((_user_variable_defs + _n_user_variables)->name, name);
  strcpy((_user_variable_defs + _n_user_variables)->ref_name, variable_name);
  (_user_variable_defs + _n_user_variables)->dim = -1;
  (_user_variable_defs + _n_user_variables)->is_variance = true;

  _n_user_variables++;

  /* Make this immediate if fields have already been defined */

  if (cs_field_n_fields() > 0)
    cs_parameters_create_added_variables();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a user property.
 *
 * \param[in]  name         name of property and associated field
 * \param[in]  dim          property dimension
 * \param[in]  location_id  id of associated mesh location
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_add_property(const char  *name,
                           int          dim,
                           int          location_id)
{
  BFT_REALLOC(_user_property_defs,
              _n_user_properties + 1,
              cs_user_property_def_t);
  BFT_MALLOC((_user_property_defs + _n_user_properties)->name,
             strlen(name) + 1,
             char);

  strcpy((_user_property_defs + _n_user_properties)->name, name);
  (_user_property_defs + _n_user_properties)->dim = dim;
  (_user_property_defs + _n_user_properties)->location_id = location_id;

  _n_user_properties++;

  /* Make this immediate if fields have already been defined */

  if (cs_field_n_fields() > 0)
    cs_parameters_create_added_properties();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined user variables not added yet.
 *
 * This number is reset to 0 when \ref cs_parameters_create_added_variables
 * is called.
 *
 * \return number of defined user variables
 */
/*----------------------------------------------------------------------------*/

int
cs_parameters_n_added_variables(void)
{
  return _n_user_variables;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined user properties not added yet.
 *
 * \return number of defined user properties
 */
/*----------------------------------------------------------------------------*/

int
cs_parameters_n_added_properties(void)
{
  return _n_user_properties;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create previously added user variables.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_create_added_variables(void)
{
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_USER;

  for (int i = 0; i < _n_user_variables; i++) {

    cs_field_t *f;

    const char *name = (_user_variable_defs + i)->name;

    int cmp_id = cs_field_id_by_name(name);

    if (cmp_id > -1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error defining user variable \"%s\";\n"
                  "this name is already reserved for field with id %d."),
                name, cmp_id);

    /* Case where we define a variance */

    if ((_user_variable_defs + i)->is_variance) {

      const char *ref_name = (_user_variable_defs + i)->ref_name;
      const cs_field_t *f_ref = cs_field_by_name_try(ref_name);

      if (f_ref == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("Error defining user variance \"%s\";\n"
                    "which refers to yet undefined variable \"%s\"."),
                  name, ref_name);

      f = cs_field_create(name,
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          f_ref->dim,
                          true);
      int k_var = cs_field_key_id("first_moment_id");
      cs_field_set_key_int(f, k_var, f_ref->id);
      cs_field_lock_key(f, k_var);
      BFT_FREE((_user_variable_defs + i)->ref_name);

    }

    /* General case */

    else {

      f = cs_field_create(name,
                          field_type,
                          CS_MESH_LOCATION_CELLS,
                          (_user_variable_defs + i)->dim,
                          true);

    }

    BFT_FREE((_user_variable_defs + i)->name);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(f, cs_field_key_id("log"), 1);
    cs_field_set_key_int(f, cs_field_key_id("post_vis"), post_flag);

  }

  BFT_FREE(_user_variable_defs);
  _n_user_variables = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create auxiliary fields for some numerical algorithms.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_auxiliary_fields(void)
{
  /* Define variable diffusivities for the temperature or
     user-defined variables */

  cs_thermal_model_t *th_model = cs_get_glob_thermal_model();
  cs_cf_model_t *th_cf_model = cs_get_glob_cf_model();

  if (th_model->has_kinetic_st == 1) {
    cs_field_t *fld
      = cs_field_create("kinetic_energy_thermal_st",
                        CS_FIELD_PROPERTY,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        true);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);
  }

  /* Fields used to model the aforepresented source term */

  if (th_model->has_kinetic_st == 1) {
    cs_field_create("rho_k_prev",
                    0,
                    CS_MESH_LOCATION_CELLS,
                    1,
                    false);

  }

  if (th_model->has_kinetic_st == 1) {
    cs_field_create("inner_face_velocity",
                    0,
                    CS_MESH_LOCATION_INTERIOR_FACES,
                    3,
                    true);
  }

  if (th_model->has_kinetic_st == 1) {
    cs_field_create("boundary_face_velocity",
                    0,
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    3,
                    true);
  }

  /* If humid air equation of state, define yv,
     the mass fraction of water vapor*/

  if (th_cf_model->ieos == CS_EOS_MOIST_AIR) {
    cs_field_t *fld = cs_field_create("yv",
                                      CS_FIELD_PROPERTY | CS_FIELD_INTENSIVE,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      true);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);
  }

  if (th_cf_model->ieos != CS_EOS_NONE) {

    /* Pressure gradient */
    if (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
        || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {

      cs_field_create("pressure_gradient",
                      0,
                      CS_MESH_LOCATION_CELLS,
                      3,
                      false);
    }

    if (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
        || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {

      cs_field_create("pressure_increment_gradient",
                      0,
                      CS_MESH_LOCATION_CELLS,
                      3,
                      false);
    }

    if (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
        || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {

      /* FIXME: check relation between this and "specific_heat" field */

      cs_field_t *fld = cs_field_create("isobaric_heat_capacity",
                                        CS_FIELD_PROPERTY | CS_FIELD_INTENSIVE,
                                        CS_MESH_LOCATION_CELLS,
                                        1,
                                        false);

      cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 0);
    }

  }

  /* Temperature in case of solving the internal energy equation */
  if (th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {
    cs_field_t *fld = cs_field_create("temperature",
                                      CS_FIELD_PROPERTY | CS_FIELD_INTENSIVE,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      true);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);
  }

  /* CFL conditions */
  if (th_model->cflt) {
    cs_field_t *fld = cs_field_create("cfl_t",
                                      CS_FIELD_POSTPROCESS,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      false);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);
  }

  if (th_model->cflp) {
    cs_field_t *fld = cs_field_create("cfl_p",
                                      CS_FIELD_POSTPROCESS,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      false);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);
  }

  /* Property fields relative to radiative transfer */
  cs_rad_transfer_add_property_fields();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create previously added user properties.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_create_added_properties(void)
{
  /* Define variable diffusivities for the temperature or
     user-defined variables */

  /* Define regular user properties */

  for (int i = 0; i < _n_user_properties; i++) {

    const char *name = (_user_property_defs + i)->name;

    int cmp_id = cs_field_id_by_name(name);

    if (cmp_id > -1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error defining user property \"%s\";\n"
                  "this name is already reserved for field with id %d."),
                name, cmp_id);

    cs_field_t *fld =
      cs_field_create(name,
                      CS_FIELD_PROPERTY | CS_FIELD_USER,
                      (_user_property_defs + i)->location_id,
                      (_user_property_defs + i)->dim,
                      false);

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
    cs_field_set_key_int(fld, cs_field_key_id("post_vis"), post_flag);

    BFT_FREE((_user_property_defs + i)->name);

  }

  BFT_FREE(_user_property_defs);
  _n_user_properties = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a boundary values field for a variable field.
 *
 * \param[in, out]  f  pointer to field structure
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_values(cs_field_t  *f)
{
  cs_field_t *bf = NULL;

  /* Check we are on cells and don't already have such value */

  if (f->location_id != CS_MESH_LOCATION_CELLS)
    return bf;

  int kbf = cs_field_key_id_try("boundary_value_id");

  int bf_id = cs_field_get_key_int(f, kbf);
  if (bf_id > -1) {
    bf = cs_field_by_id(bf_id);
    return bf;
  }

  /* Currently only managed for scalars or temperature property */

  int ks = cs_field_key_id_try("scalar_id");
  if (ks < 0)
    return bf;

  int scalar_id = (f->type & CS_FIELD_VARIABLE) ?
    cs_field_get_key_int(f, ks) : -1;

  if (scalar_id < 0&& strcmp(f->name, "temperature") != 0)
    return bf;

  /* Build new field */

  char *b_name;
  size_t l = strlen("boundary_") + strlen(f->name) + 1;
  BFT_MALLOC(b_name, l, char);
  snprintf(b_name, l, "boundary_%s", f->name);

  /* Field may already have been defined */

  bf = cs_field_by_name_try(b_name);

  if (bf == NULL) {

    int type_flag =   (f->type & (CS_FIELD_INTENSIVE | CS_FIELD_EXTENSIVE))
                    | CS_FIELD_POSTPROCESS;

    bf = cs_field_create(b_name,
                         type_flag,
                         CS_MESH_LOCATION_BOUNDARY_FACES,
                         f->dim,
                         false);

    /* Set same label as parent */

    cs_field_set_key_str(bf,
                         cs_field_key_id("label"),
                         cs_field_get_label(f));

    /* Set same postprocessing and logging defaults as parent */

    int k_log = cs_field_key_id("log");
    cs_field_set_key_int(bf,
                         k_log,
                         cs_field_get_key_int(f, k_log));

    int k_vis = cs_field_key_id("post_vis");
    int f_vis = cs_field_get_key_int(f, k_vis);
    f_vis = f_vis | CS_POST_ON_LOCATION;
    cs_field_set_key_int(bf, k_vis, f_vis);

  }
  else {

    if (   f->dim != bf->dim
        || bf->location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
      bft_error(__FILE__, __LINE__, 0,
                _("Error defining variable boundary field:\n"
                  "  parent name:   \"%s\"\n"
                  "  name:          \"%s\"\n"
                  "  dimension:     %d\n\n"
                  "An incompatible field with matching name already exists:\n"
                  "  id:          %d\n"
                  "  location_id: %d\n"
                  "  dimension:   %d"),
                f->name, bf->name, f->dim,
                bf->id, bf->location_id, bf->dim);

  }

  BFT_FREE(b_name);

  cs_field_set_key_int(f, kbf, bf->id);
  cs_field_lock_key(f, kbf);

  return bf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a boundary values field for temperature, if applicable.
 *
 * When a volume temperature variable field already exists, this amounts
 * to calling \ref cs_parameters_add_boundary_values for that field.
 * When such a variable does not exist but we have an Enthalpy variables,
 * an associated temperature boundary field is returned.
 *
 * \return  pointer to boundary values field, or NULL if not applicable
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_parameters_add_boundary_temperature(void)
{
  cs_field_t *bf = NULL;

  /* Check if we already have a temperature variable field
    (temperature or enthalpy) */

  cs_field_t *f = cs_field_by_name_try("temperature");

  if (f != NULL) //FIXME it might be not a variable as in Cooling towers
    bf = cs_parameters_add_boundary_values(f);

  else {

    f = cs_field_by_name_try("enthalpy");

    if (f != NULL) {
      if (   f->location_id != CS_MESH_LOCATION_CELLS
          || (f->type & CS_FIELD_VARIABLE) == 0)
        f = NULL;
    }

    /* If we have a compatible cell enthalpy field,
       use if to define default output options */

    if (f != NULL) {

      char b_name[] = "boundary_temperature";

      bf = cs_field_by_name_try(b_name);

      if (bf == NULL) {

        int type_flag =   (f->type & (CS_FIELD_INTENSIVE | CS_FIELD_EXTENSIVE))
                           | CS_FIELD_POSTPROCESS;

        bf = cs_field_create(b_name,
                             type_flag,
                             CS_MESH_LOCATION_BOUNDARY_FACES,
                             f->dim,
                             false);

        /* Set same postprocessing and logging defaults as enthalpy */

        int k_log = cs_field_key_id("log");
        cs_field_set_key_int(bf,
                             k_log,
                             cs_field_get_key_int(f, k_log));

        int k_vis = cs_field_key_id("post_vis");
        int f_vis = cs_field_get_key_int(f, k_vis);
        f_vis = CS_MAX(f_vis, 1);
        cs_field_set_key_int(bf, k_vis, f_vis);

      }
      else {

        if (   1 != bf->dim
            || bf->location_id != CS_MESH_LOCATION_BOUNDARY_FACES)
          bft_error
            (__FILE__, __LINE__, 0,
             _("Error defining variable \"boundary_temperature\" field:\n"
               "An incompatible field with matching name already exists:\n"
               "  id:          %d\n"
               "  location_id: %d\n"
               "  dimension:   %d"),
             bf->id, bf->location_id, bf->dim);

      }

    }

  }

  return bf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if extended neighborhood is needed.
 *
 * \return  true if extended neighborhoo is needed, false otherwise.
 */
/*----------------------------------------------------------------------------*/

bool
cs_parameters_need_extended_neighborhood(void)
{
  bool need_extended_neighborhood = false;

  /* Check if a gradient requires an extended neighborhood */

  cs_gradient_type_t  gradient_type = CS_GRADIENT_GREEN_ITER;
  cs_halo_type_t  halo_type = CS_HALO_STANDARD;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  if (halo_type == CS_HALO_EXTENDED)
    need_extended_neighborhood = true;

  else {
    const int n_fields = cs_field_n_fields();
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
        if (eqp != NULL) {
          cs_gradient_type_by_imrgra(eqp->imrgra,
                                     &gradient_type,
                                     &halo_type);
          if (halo_type == CS_HALO_EXTENDED) {
            need_extended_neighborhood = true;
            break;
          }
        }
      }
    }
  }

  /* Check for other options requiring extended neighborhood */

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {
    if (cs_at_opt_interp_is_p1_proj_needed())
      need_extended_neighborhood = true;
  }

  return need_extended_neighborhood;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Complete global parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_global_complete(void)
{
  /* Set restart_file key for various fields */
  cs_restart_set_auxiliary_field_options();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Complete general equation parameter definitions.
 *
 * Also set associated field properties such as number of associated
 * time values.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_eqp_complete(void)
{
  /* Complete settings for variable fields. */

  const int ks = cs_field_key_id_try("scalar_id");
  const int kturt = cs_field_key_id_try("turbulent_flux_model");
  const int kctheta = cs_field_key_id_try("turbulent_flux_ctheta");

  const int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (eqp != NULL) {
        if (eqp->dim != f->dim)
          eqp->dim = f->dim;

        /* Additional settings for scalars */

        int scalar_id = (ks > -1) ? cs_field_get_key_int(f, ks) : -1;
        if (scalar_id > -1) {

          /* Turbulent flux constants for GGDH, AFM and DFM */

          int turb_flux_model
            = (kturt > -1) ? cs_field_get_key_int(f, kturt) : 0;
          int turb_flux_model_type = turb_flux_model / 10;

          /* AFM and GGDH on the scalar */
          if (   turb_flux_model_type == 1
              || turb_flux_model_type == 2) {
            eqp->idften = CS_ANISOTROPIC_RIGHT_DIFFUSION;
            cs_field_set_key_double(f, kctheta, cs_turb_cthafm);
          }

          /* DFM on the scalar */
          else if (turb_flux_model_type == 3) {
            eqp->idifft = 0;
            eqp->idften = CS_ISOTROPIC_DIFFUSION;
            if (turb_flux_model == 31) {
              cs_field_set_key_double(f, kctheta, cs_turb_cthebdfm);
              cs_turb_c2trit = 0.3;
            }
            else
              cs_field_set_key_double(f, kctheta, cs_turb_cthdfm);

            /* GGDH on the thermal fluxes is automatic.
               GGDH on the variance of the thermal scalar set here. */

            cs_field_t *f_v = cs_field_get_variance(f);
            if (f_v != NULL) {
              cs_equation_param_t *eqp_v = cs_field_get_equation_param(f_v);
              eqp_v->idften = CS_ANISOTROPIC_RIGHT_DIFFUSION;
              cs_field_set_key_double(f_v, kctheta, cs_turb_csrij);
            }
          }

          /* Non-GGDH, AFM or DFM cases */
          else
            cs_field_set_key_double(f, kctheta, cs_turb_csrij);
        }
      }

      /* Harmonic face viscosity interpolation */

      if (cs_glob_space_disc->imvisf) {
        eqp->imvisf = 1;
      }

      /* Set iswdyn to 2 by default if not modified for pure
         diffusion equations */
      if (eqp->iswdyn == -1 && eqp->iconv == 0)
        eqp->iswdyn = 2;

      /* Set previous values for backward n order in time. */

      if (eqp->ibdtso > 1)
        cs_field_set_n_time_vals(f, eqp->ibdtso + 1);

    }  /* end if (f->type & CS_FIELD_VARIABLE) */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Complete general output option definitions.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_output_complete(void)
{
  const int k_log = cs_field_key_id("log");
  const int k_vis = cs_field_key_id("post_vis");

  /* Complete default settings for moment fields */

  const int n_moments = cs_time_moment_n_moments();
  if (n_moments > 0) {
    for (int m_id = 0; m_id < n_moments; m_id++) {
      cs_field_t *f = cs_time_moment_get_field(m_id);
      if (f != NULL) {
        if (cs_field_is_key_set(f, k_vis) == false) {
          int flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
          cs_field_set_key_int(f, k_vis, flag);
        }
        if (cs_field_is_key_set(f, k_log) == false)
          cs_field_set_key_int(f, k_log, 1);
      }
    }
  }

  /* Complete settings for property fields; properties
     which are hidden should already have had their settings
     forced at this stage, so should not be changed here. */

  const int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_PROPERTY) {
      if (cs_field_is_key_set(f, k_vis) == false) {
        int flag = CS_POST_ON_LOCATION;
        cs_field_set_key_int(f, k_vis, flag);
      }
      if (cs_field_is_key_set(f, k_log) == false)
        cs_field_set_key_int(f, k_log, 1);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a local variable calculation options structure,
 *        with default options.
 *
 * \return  variable calculations options structure
 */
/*----------------------------------------------------------------------------*/

cs_var_cal_opt_t
cs_parameters_var_cal_opt_default(void)
{
  return _equation_param_default;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the time scheme structure to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_scheme_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "Time discretization options\n"
                 "----------------------------\n\n"));

  /* Frozen velocity field */
  if (cs_glob_time_scheme->iccvfg) {
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Frozen velocity field\n\n"));
    return;
  }

  cs_log_printf
    (CS_LOG_SETUP,
     _("  Time scheme:\n\n"
       "    time_order:  %d (order of base time stepping scheme)\n"
       "    istmpf:      %d (time order of the mass flux scheme)\n"
       "    isno2t:      %d (time scheme for the momentum source terms,\n"
       "                     apart from convection and diffusion)\n"
       "    isto2t:      %d (time scheme for the turbulence source terms)\n"),
     cs_glob_time_scheme->time_order,
     cs_glob_time_scheme->istmpf,
     cs_glob_time_scheme->isno2t,
     cs_glob_time_scheme->isto2t);

  if (cs_glob_time_scheme->isno2t > 0)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    thetsn:      %g (theta_S for Navier-Stokes source terms)\n"),
       cs_glob_time_scheme->thetsn);

  if (cs_glob_time_scheme->isto2t > 0)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    thetst:      %g (theta for turbulence explicit source terms)\n"),
       cs_glob_time_scheme->thetst);

  int key_t_ext_id = cs_field_key_id("time_extrapolated");

  {
    int iviext = 0;
    const cs_field_t *f = cs_field_by_name_try("molecular_viscosity");
    if (f != NULL)
      iviext = cs_field_get_key_int(f, key_t_ext_id);
    if (iviext > 0)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    thetvi:      %g (theta for total viscosity)\n"),
         cs_glob_time_scheme->thetvi);
  }

  {
    int icpext = 0;
    const cs_field_t *f = cs_field_by_name_try("specific_heat");
    if (f != NULL)
      icpext = cs_field_get_key_int(f, key_t_ext_id);
    if (icpext > 0)
      cs_log_printf
        (CS_LOG_SETUP,
         _("    thetvi:      %g (theta for specific heat)\n"),
         cs_glob_time_scheme->thetcp);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the space discretization structure to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_space_disc_log_setup(void)
{
  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "Space discretization options\n"
                 "----------------------------\n\n"));

  const char *imvisf_value_str[] = {N_("arithmetic"),
                                    N_("harmonic")};
  const char *halo_type_str[] = {N_("face neighbors"),
                                 N_("extended neighborhood")};

  cs_log_printf(CS_LOG_SETUP,
                _("    imvisf:    %d (%s face viscosity field interpolation)\n"),
                cs_glob_space_disc->imvisf,
                _(imvisf_value_str[cs_glob_space_disc->imvisf]));

  cs_gradient_type_t  gradient_type = CS_GRADIENT_GREEN_ITER;
  cs_halo_type_t  halo_type = CS_HALO_STANDARD;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "    imrgra:    %d (gradient reconstruction:\n"
                  "                  %s,\n"
                  "                  using %s)\n"),
                cs_glob_space_disc->imrgra,
                _(cs_gradient_type_name[gradient_type]),
                _(halo_type_str[halo_type]));

  const char *iflxmw_value_str[]
    = {N_("0 (based on mesh velocity at cell centers)"),
       N_("1 (based on nodes displacement)")};

  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "    Method to compute inner mass flux due to mesh"
                 " velocity in ALE\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    iflxmw:    %s\n"),
                _(iflxmw_value_str[cs_glob_space_disc->iflxmw]));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

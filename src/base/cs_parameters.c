/*============================================================================
 * General parameters management.
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
#include "cs_ale.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_parall.h"
#include "cs_parameters_check.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_rad_transfer_fields.h"
#include "cs_syr_coupling.h"
#include "cs_turbulence_model.h"
#include "cs_time_moment.h"
#include "cs_thermal_model.h"
#include "cs_cf_model.h"
#include "cs_sat_coupling.h"
#include "cs_tree.h"
#include "cs_turbomachinery.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_wall_functions.h"

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
             \f$\theta\f$ being given by the value of \ref theta associated
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
             \f$\theta\f$ being given by the value of \ref theta associated
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

  double  theta;
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

/* Default equation param */

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
   .blencv = 1.,
   .blend_st = 0.,
   .epsilo = 1.e-5,
   .epsrsm = 1.e-4,
   .epsrgr = 1.e-4,
   .climgr = 1.5,
   .relaxv = 1.,
   .b_gradient_r = 2,

   .default_bc = CS_BC_SYMMETRY,
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
   .adv_strategy = CS_PARAM_N_ADVECTION_STRATEGIES,
   .adv_extrapol = CS_PARAM_N_ADVECTION_EXTRAPOLATIONS,
   .upwind_portion = 0.,
   .cip_scaling_coef = -1.0,
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
   .saddle_param = NULL,

   .incremental_algo_type = CS_PARAM_N_NL_ALGOS,
   .incremental_algo_cvg =
   {.atol = -1., .rtol = -1., .dtol = -1., .n_max_iter = -1},
   .incremental_relax_factor = -1.,
   .incremental_anderson_param = {.n_max_dir = 0, .starting_iter = 0,
     .max_cond = -1., .beta = 0., .dp_type = CS_PARAM_N_DOTPROD_TYPES }

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

/* Flags used for re-initializing rho, cp and mu during restart if needed */

static int _initvi = 0;
static int _initro = 0;
static int _initcp = 0;

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
                              int     **iccvfg,
                              int     **initro);

void
cs_f_restart_auxiliary_get_pointers(int  **ileaux);

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
  cs_log_printf(CS_LOG_SETUP, fmt_r, "theta", _t->theta);
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
  cs_log_printf(CS_LOG_SETUP, fmt_r, "theta", _t->theta,
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

  eqp->theta = vcopt->theta;
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
                              int     **iccvfg,
                              int     **initro)
{
  *ischtp = &(_time_scheme.time_order);
  *istmpf = &(_time_scheme.istmpf);
  *isno2t = &(_time_scheme.isno2t);
  *isto2t = &(_time_scheme.isto2t);
  *iccvfg = &(_time_scheme.iccvfg);

  *initro = &_initro;
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global restart_auxiliary structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ileaux  --> pointer to cs_glob_restart_auxiliary->read_auxiliary
 *----------------------------------------------------------------------------*/

void
cs_f_restart_auxiliary_get_pointers(int  **ileaux)
{
  *ileaux = &(_restart_auxiliary.read_auxiliary);
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

  vcopt->theta = eqp->theta;
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
 * \brief Set init state to 1. This is necessary for fortran mapping and
 * should be changed in the future.
 *
 * \param[in] idx id of variable. 0 is viscosity, 1 density, 2 heat capacity.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_set_init_state_on(int idx)
{
  assert(idx >= 0 && idx < 3);

  if (idx == 0)
    _initvi = 1;
  else if (idx == 1)
    _initro = 1;
  else if (idx == 2)
    _initcp = 1;
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
  cs_field_define_key_int("coupled_with_vel_p", 0, CS_FIELD_VARIABLE);

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
     (cs_field_clear_key_struct_t *)cs_equation_param_clear,
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

    if (f->dim == 3)
      cs_field_set_key_int(f, cs_field_key_id("coupled"), 1);

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
     the mass fraction of water vapor */

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

      cs_field_create("algo:pressure_gradient",
                      0,
                      CS_MESH_LOCATION_CELLS,
                      3,
                      false);

      cs_field_create("algo:pressure_increment_gradient",
                      0,
                      CS_MESH_LOCATION_CELLS,
                      3,
                      false);

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

  const int key_t_ext_id = cs_field_key_id("time_extrapolated");
  const int kthetss = cs_field_key_id("st_exp_extrapolated");
  const int kthetvs = cs_field_key_id("diffusivity_extrapolated");
  const int kisso2t = cs_field_key_id("scalar_time_scheme");
  const int kscal = cs_field_key_id("scalar_id");
  const int kivisl = cs_field_key_id("diffusivity_id");

  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_rho_b = CS_F_(rho_b);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_cp = CS_F_(cp);

  /* Logging and postprocessing output */
  if (cs_glob_fluid_properties->irovar < 1) {
    cs_field_set_key_int(f_rho, cs_field_key_id("post_vis"), 0);
    cs_field_set_key_int(f_rho, cs_field_key_id("log"), 0);
    cs_field_set_key_int(f_rho_b, cs_field_key_id("post_vis"), 0);
    cs_field_set_key_int(f_rho_b, cs_field_key_id("log"), 0);
  }
  if (cs_glob_fluid_properties->ivivar < 1) {
    cs_field_set_key_int(f_mu, cs_field_key_id("post_vis"), 0);
    cs_field_set_key_int(f_mu, cs_field_key_id("log"), 0);
  }

  /* Time scheme and time stepping */

  cs_time_step_t *ts = cs_get_glob_time_step();
  cs_time_scheme_t *time_scheme = cs_get_glob_time_scheme();

  if (ts->nt_max == -1 && ts->t_max < -0.5) {
    cs_time_step_define_nt_max(10);
  }

  /* Physical properties */
  int iviext = cs_field_get_key_int(f_mu, key_t_ext_id);
  if (fabs(time_scheme->thetvi+999.) > cs_math_epzero) {
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _("in the data specification"),
       _("iviext = %d\n"
         "thetvi will be initialized automatically.\n"
         "Do not modify it.\n"),
       iviext);
  }
  else if (iviext == 0)
    time_scheme->thetvi = 0.;
  else if (iviext == 1)
    time_scheme->thetvi = 0.5;
  else if (iviext == 2)
    time_scheme->thetvi = 1.;

  if (cs_glob_fluid_properties->icp >= 0) {
    int icpext = cs_field_get_key_int(f_cp, key_t_ext_id);
    if (fabs(time_scheme->thetcp+999.) > cs_math_epzero) {
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("in the data specification"),
         _("icpext = %d\n"
           "thetcp will be initialized automatically.\n"
           "Do not modify it.\n"),
         icpext);
    }
    else if (icpext == 0)
      time_scheme->thetcp = 0.;
    else if (icpext == 1)
      time_scheme->thetcp = 0.5;
    else if (icpext == 2)
      time_scheme->thetcp = 1.;
  }

  /* NS source terms */
  int isno2t = cs_glob_time_scheme->isno2t;
  if (fabs(time_scheme->thetsn+999.) > cs_math_epzero) {
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _("in the data specification"),
       _("isno2t = %d\n"
         "thetsn will be initialized automatically.\n"
         "Do not modify it.\n"),
       isno2t);
  }
  else if (isno2t == 1)
    time_scheme->thetsn = 0.5;
  else if (isno2t == 2)
    time_scheme->thetsn = 1.;
  else if (isno2t == 0)
    time_scheme->thetsn = 0.;

  /* Turbulent variables source terms */
  int isto2t = cs_glob_time_scheme->isto2t;
  if (fabs(time_scheme->thetst+999.) > cs_math_epzero) {
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _("in the data specification"),
       _("isto2t = %d\n"
         "thetst will be initialized automatically.\n"
         "Do not modify it.\n"),
       isto2t);
  }
  else if (isto2t == 1)
    time_scheme->thetst = 0.5;
  else if (isto2t == 2)
    time_scheme->thetst = 1.;
  else if (isto2t == 0)
    time_scheme->thetst = 0.;

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    int scalar_id = cs_field_get_key_int(f, kscal);
    if (scalar_id > -1) {
      /* Scalar source terms */
      int isso2t = cs_field_get_key_int(f, kisso2t);
      cs_real_t thetss = cs_field_get_key_double(f, kthetss);

      if (fabs(thetss+1.) > cs_math_epzero) {
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("in the data specification"),
           _("Scalar (%s) isso2t = %d\n"
             "thetss will be initialized automatically.\n"
             "Do not modify it.\n"),
           f->name, isso2t);
      }
      else if (isso2t == 1) {
        thetss = 0.5;
        cs_field_set_key_double(f, kthetss, thetss);
      }
      else if (isso2t == 2) {
        thetss = 1.;
        cs_field_set_key_double(f, kthetss, thetss);
      }
      else if (isso2t == 0) {
        thetss = 0.;
        cs_field_set_key_double(f, kthetss, thetss);
      }
      /* Scalar diffusivity */
      int f_id_d = cs_field_get_key_int(f, kivisl);
      if (f_id_d >= 0.)
        iviext = cs_field_get_key_int(cs_field_by_id(f_id_d), key_t_ext_id);
      else
        iviext = 0;

      cs_real_t thetvs = cs_field_get_key_double(f, kthetvs);
      if (fabs(thetvs+1.) > cs_math_epzero) {
        cs_parameters_error
          (CS_ABORT_DELAYED,
           _("in the data specification"),
           _("Scalar (%s) iviext = %d\n"
             "thetvs will be initialized automatically.\n"
             "Do not modify it.\n"),
           f->name, iviext);
      }
      else if (iviext == 0) {
        thetvs = 0.;
        cs_field_set_key_double(f, kthetvs, thetvs);
      }
      else if (iviext == 1) {
        thetvs = 0.5;
        cs_field_set_key_double(f, kthetvs, thetvs);
      }
      else if (iviext == 2) {
        thetvs = 1.;
        cs_field_set_key_double(f, kthetvs, thetvs);
      }

    }
  }
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
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kturt = cs_field_key_id_try("turbulent_flux_model");
  const int kctheta = cs_field_key_id_try("turbulent_flux_ctheta");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int kscacp  = cs_field_key_id("is_temperature");

  const int kcpsyr = cs_field_key_id("syrthes_coupling");
  const int kclvfl = cs_field_key_id("variance_clipping");
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_p = CS_F_(p);
  cs_field_t *f_t = cs_thermal_model_field();
  cs_equation_param_t *eqp_vel
    = (f_vel != NULL) ? cs_field_get_equation_param(f_vel) : NULL;
  cs_equation_param_t *eqp_p
    = (f_p != NULL) ? cs_field_get_equation_param(f_p) : NULL;
  cs_time_step_t *ts = cs_get_glob_time_step();
  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  cs_velocity_pressure_param_t *vp_param
    = cs_get_glob_velocity_pressure_param();
  cs_velocity_pressure_model_t *vp_model
    = cs_get_glob_velocity_pressure_model();
  cs_turb_rans_model_t *rans_mdl = cs_get_glob_turb_rans_model();
  cs_vof_parameters_t *vof_param = cs_get_glob_vof_parameters();
  cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();
  cs_turbomachinery_model_t iturbo = cs_turbomachinery_get_model();

  cs_real_t *xyzp0 = fluid_props->xyzp0;

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

      if (fabs(eqp->theta+1.) > cs_math_epzero) {
        cs_parameters_error
          (CS_WARNING,
           _("advanced modification for"),
           _("(%s) of the variable theta.\n"),
           f->name);
      }
      else {
        int time_order = cs_glob_time_scheme->time_order;
        /* For the pressure, no theta-scheme */
        if (f->id == CS_F_(p)->id)
          eqp->theta = 1.;
        else if (eqp->istat == 0)
          eqp->theta = 1.;
        else if (time_order == 1)
          eqp->theta = 1.;
        else if (time_order == 2)
          eqp->theta = 0.5;
      }

    }  /* end if (f->type & CS_FIELD_VARIABLE) */
  }

  /* Diffusivity model */
  if (cs_glob_turb_model->itytur == 3) {
    cs_field_t *f_rij = CS_F_(rij);
    cs_field_t *f_eps = CS_F_(eps);
    cs_equation_param_t *eqp_rij = cs_field_get_equation_param(f_rij);
    cs_equation_param_t *eqp_eps = cs_field_get_equation_param(f_eps);
    /* Daly harlow (GGDH) on Rij and epsilon by default */
    if (cs_glob_turb_rans_model->idirsm != 0) {
      eqp_rij->idften = CS_ANISOTROPIC_RIGHT_DIFFUSION;
      eqp_eps->idften = CS_ANISOTROPIC_RIGHT_DIFFUSION;
      /* Scalar diffusivity (Shir model) elsewhere (idirsm = 0) */
    }
    else {
      eqp_rij->idften = CS_ISOTROPIC_DIFFUSION;
      eqp_eps->idften = CS_ISOTROPIC_DIFFUSION;
    }
  }

  /* ISSTPC
   * If the user has not specified anything for the slope test (=-1),
   * We impose 1 (i.e. without) for the velocity for LES
   *           0 (i.e. with) otherwise */

  if (cs_glob_turb_model->itytur == 4) {
    if (eqp_vel != NULL) {
      if (eqp_vel->isstpc == -999)
        eqp_vel->isstpc = 1;
    }
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int scalar_id = cs_field_get_key_int(f, ks);
      if (scalar_id > -1) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        if (eqp->isstpc == -999) eqp->isstpc = 0;
      }
    }
  }

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (eqp->isstpc == -999) eqp->isstpc = 0;
    }
  }

  /* BLENCV
   * If the user has not specified anything for the convective scheme
   * We impose 1 (i.e. centered) for the velocities
   *                                 the user scalars
   *                                 the thermal scalar
   *           0 (i.e. pure upwind) for the rest
   * (especially, in LES, all the variables are centered)
   *
   * For the cavitation model, we force in all cases the void fraction in upwind
   * and e display a message if the user has specified something else. */

  if (eqp_vel != NULL) {
    if (fabs(eqp_vel->blencv + 1.) < cs_math_epzero)
      eqp_vel->blencv = 1.;
  }

  if (cs_glob_vof_parameters->vof_model & CS_VOF_MERKLE_MASS_TRANSFER) {
    cs_field_t *f_void_f = CS_F_(void_f);
    cs_equation_param_t *eqp_void_f = cs_field_get_equation_param(f_void_f);
    if (fabs(eqp_void_f->blencv + 1.) < cs_math_epzero) {
      if (fabs(eqp_void_f->blencv + 1.) > cs_math_epzero) {
        cs_parameters_error
          (CS_WARNING,
           _("in the data specification"),
           _("The cavitation model requires an upwind convection scheme"
             "for the void fraction (blencv(void_f) = %g\n"
             "The user has set blencv(void_f) = %g.\n"
             "The upwind scheme for the void fraction is forced.\n"),
           0., eqp_void_f->blencv);
      }
      eqp_void_f->blencv = 0.;
    }
  }
  else if (cs_glob_vof_parameters->vof_model > 0) {
    cs_field_t *f_void_f = CS_F_(void_f);
    cs_equation_param_t *eqp_void_f = cs_field_get_equation_param(f_void_f);
    if (fabs(eqp_void_f->blencv + 1.) < cs_math_epzero) {
      eqp_void_f->blencv = 1.;
    }
  }

  if (f_t != NULL) {
    cs_equation_param_t *eqp_t = cs_field_get_equation_param(f_t);
    if (fabs(eqp_t->blencv + 1.) < cs_math_epzero) eqp_t->blencv = 1.;
  }

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int scalar_id = (f->type & CS_FIELD_USER) ?
      cs_field_get_key_int(f, ks) : -1;
    if (scalar_id > -1) {
      cs_equation_param_t *eqp_sca = cs_field_get_equation_param(f);
      if (fabs(eqp_sca->blencv + 1.) < cs_math_epzero) {
        eqp_sca->blencv = 1.;
      }
    }
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (fabs(eqp->blencv + 1.) < cs_math_epzero) {
        eqp->blencv = 0.;
      }
    }
  }

  /* nswrsm, epsrsm and epsilo:
   *
   * If the user has not specified anything (nswrsm=-1),
   * We impose
   *  at order 1:
   *        2 for the pressure
   *        1 for the other variables
   *        We initialize epsilo to 1.e-8 for the pressure
   *                                1.e-5 for the other variables
   *                      epsrsm to 10*epsilo
   *  at order 2:
   *        5 for the pressure
   *        10 for the other variables
   *        We initialize epsilo to 1.e-5
   *                      epsrsm to 10*epsilo
   */

  if (cs_glob_time_scheme->time_order == 2) {
    if (eqp_p != NULL) {
      if (eqp_p->nswrsm == -1) eqp_p->nswrsm = 5;
      if (fabs(eqp_p->epsilo + 1.) < cs_math_epzero)
      eqp_p->epsilo = 1.e-5;
      if (fabs(eqp_p->epsrsm + 1.) < cs_math_epzero)
        eqp_p->epsrsm = 10.*eqp_p->epsilo;
    }
    if (eqp_vel != NULL) {
      if (eqp_vel->nswrsm == -1) eqp_vel->nswrsm = 10;
      if (fabs(eqp_vel->epsilo + 1.) < cs_math_epzero)
        eqp_vel->epsilo = 1.e-5;
      if (fabs(eqp_vel->epsrsm + 1.) < cs_math_epzero)
        eqp_vel->epsrsm = 10.*eqp_vel->epsilo;
    }
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      int scalar_id = cs_field_get_key_int(f, ks);
      if (scalar_id > -1) {
        cs_equation_param_t *eqp_sca = cs_field_get_equation_param(f);
        if (eqp_sca->nswrsm == -1) eqp_sca->nswrsm = 10;
        if (fabs(eqp_sca->epsilo + 1.) < cs_math_epzero)
          eqp_sca->epsilo = 1.e-5;
        if (fabs(eqp_sca->epsrsm + 1.) < cs_math_epzero)
          eqp_sca->epsrsm = 10.*eqp_sca->epsilo;
      }
    }
  }

  /* For the pressure, default solver precision 1e-8
   * because the mass conservation is up to this precision. */
  if (eqp_p != NULL) {
    if (eqp_p->nswrsm == -1)
      eqp_p->nswrsm = 2;
    if (fabs(eqp_p->epsilo + 1.) < cs_math_epzero)
      eqp_p->epsilo = 1.e-8;
  }
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (eqp->nswrsm == -1) eqp->nswrsm = 1;
      if (fabs(eqp->epsilo + 1.) < cs_math_epzero) eqp->epsilo = 1.e-5;
      if (fabs(eqp->epsrsm + 1.) < cs_math_epzero)
        eqp->epsrsm = 10.*eqp->epsilo;
    }
  }

  /* dtmin dtmax cdtvar */
  if (time_opt->dtmin <= -cs_math_big_r) time_opt->dtmin = 0.1*ts->dt_ref;
  if (time_opt->dtmax <= -cs_math_big_r) time_opt->dtmax = 1000.*ts->dt_ref;

  /* For laminar cases or when using low Reynolds model: no wall function.
   * When using mixing length, Spalart-Allmaras or LES: one scale log law.
   * When using EB-RSM : all y+ wall functions
   * In all other cases: 2 scales log law.
   * Here iwallf is set automatically only if it was not set in the gui or
   * in a user subroutine. */

  cs_wall_functions_t *wall_fns = cs_get_glob_wall_functions();
  if (wall_fns->iwallf == CS_WALL_F_UNSET) {
    if (   cs_glob_turb_model->iturb == CS_TURB_MIXING_LENGTH
        || cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS
        || cs_glob_turb_model->itytur == 4) {
      wall_fns->iwallf = CS_WALL_F_1SCALE_LOG;
    }
    else if (   cs_glob_turb_model->iturb == CS_TURB_NONE
             || cs_glob_turb_model->itytur == 5) {
      wall_fns->iwallf = CS_WALL_F_DISABLED;
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      wall_fns->iwallf = CS_WALL_F_2SCALES_CONTINUOUS;
    }
    else {
      wall_fns->iwallf = CS_WALL_F_2SCALES_LOG;
    }
  }

  /* If the wall function for the velocity is the two scales wall function using
   * Van Driest mixing length (iwallf = CS_WALL_F_2SCALES_VDRIEST), then the
   * corresponding wall function for scalar should be used
   * (iwalfs = CS_WALL_F_S_VDRIEST).
   * For atmospheric flows, it is by default Louis, or Monin-Obukhov
   * Here iwalfs is set automatically only if it was not set in a user
   * subroutine. */

  if (wall_fns->iwalfs == CS_WALL_F_S_UNSET) {
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0)
      wall_fns->iwalfs = CS_WALL_F_S_LOUIS;
    else if (wall_fns->iwallf == CS_WALL_F_2SCALES_VDRIEST)
      wall_fns->iwalfs = CS_WALL_F_S_VDRIEST;
    else
      wall_fns->iwalfs = CS_WALL_F_S_ARPACI_LARSEN;
  }

  /* ypluli
   * 1/xkappa is the value that guarantees the continuity of the derivative
   * between the linear and logarithmic zones.
   *
   * In the case of invariant wall functions, we use the value of continuity
   * of the velocity profile, 10.88.
   *
   * For LES, we put back 10.88 so as to avoid checkerboard effect when we are
   * at the limit (for one scale mode indeed, ypluli=1/xkappa does not allow
   * one to compute u* in a completely satisfying manner.
   * Idem with Spalart-Allmaras.*/

  if (wall_fns->ypluli < -cs_math_big_r) {
    if (   wall_fns->iwallf == CS_WALL_F_SCALABLE_2SCALES_LOG
        || cs_glob_turb_model->itytur == 4
        || cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS
        || wall_fns->iwallf == CS_WALL_F_2SCALES_SMOOTH_ROUGH
        || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA
        || cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
      wall_fns->ypluli = 10.88;
    }
    else {
      wall_fns->ypluli = 1./cs_turb_xkappa;
    }
  }

  /* If the user did not modify icpsyr, we take by default:
   *  if there is not coupling
   *    0 for all the scalars
   *  otherwise
   *    1 for the thermal scalar if it exists
   *    0 for the others
   * The convenient modifications shall be added for the particular physics.
   * The consistency tests will be done in verini. */

  /* Count the number of scalars nscal */
  int nscal = cs_field_n_scalar_fields();

  if (nscal > 0) {

    /* We check if there is coupling */
    int nbccou = cs_syr_coupling_n_couplings();

    /* If there is coupling */
    if (nbccou != 0) {

      /* We count the number of coupled scalars */
      int nscacp = 0;
      for (int f_id = 0; f_id < n_fields; f_id ++) {
        cs_field_t *f_sca = cs_field_by_id(f_id);
        int isca = cs_field_get_key_int(f_sca, ks);
        if (isca > 0) {
          int icpsyr = cs_field_get_key_int(f_sca, kcpsyr);
          if (icpsyr == 1) nscacp++;
        }
      }

      /* If the user has not coupled any scalar */
      if (nscacp == 0) {
        /* We couple the temperature scalar of the phase */
        if (f_t != NULL)
          cs_field_set_key_int(f_t, kcpsyr, 1);
      }
    }
  }

  /* is_temperature
   * If the user has not modified "is_temperature", we take by default:
   *  passive scalar of scalars other than the temperature scalar
   *      = 0 : passive, enthalpy, or energy
   *      = 1 : temperature */

  if (nscal > 0) {
    for (int f_id = 0; f_id < n_fields; f_id ++) {
      cs_field_t *f_sca = cs_field_by_id(f_id);
      int isca = cs_field_get_key_int(f_sca, ks);
      if (isca > 0) {
        int iscacp = cs_field_get_key_int(f_sca, kscacp);
        if (iscacp == -1) {
          if (f_sca == f_t && (cs_glob_thermal_model->thermal_variable
                               == CS_THERMAL_MODEL_TEMPERATURE)) {
            iscacp = 1;
          }
          else {
            iscacp = 0;
          }
          cs_field_set_key_int(f_sca, kscacp, iscacp);
        }
      }
    }
  }

  /* Compute the hydrostatic pressure at outlet for the Dirichlet conditions
   * on pressure. It is deduced from iphydr and the gravity value (arbitrary
   * test on the norm).
   * icalhy is initialized at -1 (the user may have forced 0 or 1 and in this
   * case, we do not change it. */

  if (vp_param->icalhy != -1 && vp_param->icalhy != 0
                             && vp_param->icalhy != 1 ) {
    cs_parameters_is_in_range_int(CS_ABORT_DELAYED,
                                  _("while reading input data"),
                                  "cs_glob_velocity_pressure_param->icalhy",
                                  vp_param->icalhy,
                                  0, 1);
  }

  /* If the fluid_solid option is enabled, we force ikecou to 0. */
  if (vp_model->fluid_solid) {
    if (rans_mdl->ikecou == 1) {
      rans_mdl->ikecou = 0;
      cs_parameters_error
        (CS_WARNING,
         _("in the data specification"),
         _("The pseudo coupling of turbulent dissipation and turbulent kinetic"
           "energy (ikecou = 1) is not compatible with the use of fluid/solid"
           "option to disable the dynamic in the solid cells"
           "(fluid_solid = 1).\n"
           "The parameter ikecou is forced to 0 (no coupling).\n"
           "The calculation will be run.\n"));
    }
  }

  /* relaxv */

  if (time_opt->idtvar < 0) {
    cs_real_t relxsp = 1.-time_opt->relxst;
    if (relxsp <= cs_math_epzero) relxsp = time_opt->relxst;
    if (eqp_p != NULL) {
      if (fabs(eqp_p->relaxv+1.) <= cs_math_epzero)
        eqp_p->relaxv = relxsp;
    }
    for (int f_id = 0; f_id < n_fields; f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        if (fabs(eqp->relaxv+1.) <= cs_math_epzero)
          eqp->relaxv = time_opt->relxst;
      }
    }
  }
  else {
    for (int f_id = 0; f_id < n_fields; f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        if (fabs(eqp->relaxv+1.) <= cs_math_epzero) eqp->relaxv = 1.;
      }
    }
  }

  /* Options specific to steady case */

  if (time_opt->idtvar < 0) {
    vp_param->ipucou = 0;
    ts->dt_ref = 1.;
    time_opt->dtmin = 1.;
    time_opt->dtmax = 1.;
    for (int f_id = 0; f_id < n_fields; f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        cs_equation_param_t *eqp = cs_field_get_equation_param(f);
        eqp->istat = 0;
      }
    }
    if (eqp_vel != NULL)
      vp_param->arak = vp_param->arak/CS_MAX(eqp_vel->relaxv, cs_math_epzero);
  }

  /* With a staggered approach, no Rhie and Chow correction is needed. */
  if (vp_param->staggered == 1) vp_param->arak = 0.;

  /* Physical constant tables */

  /* iclvfl
   * If the user has not modified iclvfl, we take by default:
   *  0 for the variances
   * The convenient modifications shall be added for the particular physics.
   * If the user gives a value, we put iclvfl to 2. */

  for (int f_id = 0; f_id < n_fields; f_id ++) {
    cs_field_t *f_sca = cs_field_by_id(f_id);
    int isca = cs_field_get_key_int(f_sca, ks);
    if (isca > 0) {
      int iscavr = cs_field_get_key_int(f_sca, kscavr);
      if (iscavr > 0) {
        int iclvfl = cs_field_get_key_int(f_sca, kclvfl);
        /* Get the min clipping */
        cs_real_t scminp = cs_field_get_key_double(f_sca, kscmin);
        /* If modified put 2 */
        if (   iclvfl == -1
            && fabs(scminp+cs_math_big_r) >= cs_math_epzero) {
          cs_field_set_key_int(f_sca, kclvfl, 2);
        }
        else if (iclvfl == -1) {
          cs_field_set_key_int(f_sca, kclvfl, 0);
        }
        /* Minimum for variances is 0 or greater */
        /* Set min clipping to 0 */
        scminp = CS_MAX(0., scminp);
        cs_field_set_key_double(f_sca, kscmin, scminp);
      }
    }
  }

  for (int f_id = 0; f_id < n_fields; f_id ++) {
    cs_field_t *f_sca = cs_field_by_id(f_id);
    int isca = cs_field_get_key_int(f_sca, ks);
    if (isca > 0) {
      cs_real_t visls_0 = cs_field_get_key_double(f_sca, kvisl0);
      /* For scalars which are not variances, define the reference
       * diffusivity */
      int iscavr = cs_field_get_key_int(f_sca, kscavr);
      if (iscavr <= 0 && visls_0 < -cs_math_big_r) {
        int iscacp = cs_field_get_key_int(f_sca, kscacp);
        if (iscacp > 0) {
          /* For temperature, the diffusivity factor is directly the thermal
           * conductivity lambda = Cp * mu / Pr
           * where Pr is the (molecular) Prandtl number */
          visls_0 = fluid_props->viscl0 * fluid_props->cp0;
        }
        else {
          visls_0 = fluid_props->viscl0;
        }
        cs_field_set_key_double(f_sca, kvisl0, visls_0);
      }

      /* For fluctuation variances, the diffusivity is that of the associated
       * scalar. */
      if (iscavr >= 0) {
        cs_field_t *f_ref = cs_field_by_id(iscavr);
        visls_0 = cs_field_get_key_double(f_ref, kvisl0);
        cs_real_t visls_cmp = cs_field_get_key_double(f_sca, kvisl0);
        cs_field_set_key_double(f_sca, kvisl0, visls_0);
        if (visls_cmp > -cs_math_big_r) {
          cs_parameters_error
            (CS_WARNING,
             _("in the data specification"),
             _("The scalar %s is the fluctuation variance of the scalar  %d\n"
               "The diffusivity_ref value of the scalar %s must not be set:\n"
               "it is automatically set equal to the scalar diffusivity %d"
               "i.e. %g\n"),
             f_sca->name, iscavr, f_sca->name, iscavr, visls_0);
        }
      }
    }
  }

  /* xyzp0 : reference point for hydrostatic pressure
   * The user should specify the 3 coordinates, otherwise
   * it is set to (0.,0.,0.). */

  if (   xyzp0[0] > -0.5*cs_math_infinite_r
      && xyzp0[1] > -0.5*cs_math_infinite_r
      && xyzp0[2] > -0.5*cs_math_infinite_r) {
    fluid_props->ixyzp0 = 1;
  }
  else {
    for (int ii = 0; ii < 3; ii++) {
      xyzp0[ii] = 0.;
    }
  }

  /* VoF model enabled */
  if (cs_glob_vof_parameters->vof_model > 0) {
    if (vof_param->rho2 > vof_param->rho1) {
      fluid_props->ro0    = vof_param->rho2;
      fluid_props->viscl0 = vof_param->mu2;
    }
    else {
      fluid_props->ro0    = vof_param->rho1;
      fluid_props->viscl0 = vof_param->mu1;
    }

    /* VOF algorithm: continuity of the flux across internal faces */
    if (eqp_p != NULL)
      eqp_p->imvisf = 1;
  }

  /* Elements of albase */
  if (cs_glob_ale > CS_ALE_NONE) {
    if (!cs_restart_present() && cs_glob_ale_need_init == -999)
      cs_glob_ale_need_init = 1;
  }
  else {
    cs_glob_ale_need_init = 0;
  }

  /* Global Parameters */

  const int nbrcpl = cs_sat_coupling_n_couplings();
  if (nbrcpl >= 1 && iturbo != 0) {
    cs_glob_sat_coupling_face_interpolation_type = 1;
  }

  /* Define min/max clipping values of void fraction of VoF model */
  if (cs_glob_vof_parameters->vof_model > 0) {
    cs_real_t clvfmn = cs_field_get_key_double(CS_F_(void_f), kscmin);
    cs_real_t clvfmx = cs_field_get_key_double(CS_F_(void_f), kscmax);

    if (clvfmn < -0.5*cs_math_big_r) {
      clvfmn = 0.;
      if (cs_glob_vof_parameters->vof_model & CS_VOF_MERKLE_MASS_TRANSFER) {
        clvfmn = cs_math_epzero;
      }
    }
    if (clvfmx > 0.5*cs_math_big_r) {
      clvfmx = 1.;
      if (cs_glob_vof_parameters->vof_model & CS_VOF_MERKLE_MASS_TRANSFER) {
        clvfmx = 1.-cs_math_epzero;
      }
    }

    cs_field_set_key_double(CS_F_(void_f), kscmin, clvfmn);
    cs_field_set_key_double(CS_F_(void_f), kscmax, clvfmx);
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

    /* Build clipping field for post-processing */
    if (   f->type & CS_FIELD_VARIABLE
        && !(f->type & CS_FIELD_CDO)) {

      int k_clipping_id = cs_field_key_id("clipping_id");
      int clip_id = cs_field_get_key_int(f, k_clipping_id);

      if (clip_id >= 0) {

        /* Now create matching field
           Build name and label */

        int field_type = CS_FIELD_INTENSIVE | CS_FIELD_POSTPROCESS;
        cs_field_t *f_c
          = cs_field_create_by_composite_name(f->name,
                                              "clipped",
                                              field_type,
                                              CS_MESH_LOCATION_CELLS,
                                              f->dim,
                                              false);


        cs_field_set_key_int(f_c, k_vis, CS_POST_ON_LOCATION);
        cs_field_set_key_int(f, k_clipping_id, f_c->id);

      }

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a local equation param structure,
 *        with default options.
 *
 * \return  equation param structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t
cs_parameters_equation_param_default(void)
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
         _("    thetcp:      %g (theta for specific heat)\n"),
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
                _("      iflxmw:    %s\n"),
                _(iflxmw_value_str[cs_glob_space_disc->iflxmw]));

  const char *itbrrb_value_str[] = {N_("yes (with reconstruction)"),
                                    N_("no (without reconstruction)")};
  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "    Accurate BCs for temperature/enthalpy or"
                 " scalars with boundary field\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("      itbrrb:    %s\n"),
                _(itbrrb_value_str[cs_glob_space_disc->itbrrb]));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

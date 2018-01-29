/*============================================================================
 * General parameters management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_convection_diffusion.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_post.h"
#include "cs_parall.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_mesh_location.h"
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
  \enum cs_parameter_error_behavior_t

  \brief File acces modes

  \var CS_WARNING
       Warn only
  \var CS_ABORT_DELAYED
       Abort when \ref cs_parameters_error_barrier is called.
  \var CS_FILE_MODE_APPEND
       Abort immediately

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
        - 0: iterative process
        - 1: standard least squares method
        - 2: least squares method with extended neighborhood
        - 3: least squares method with reduced extended neighborhood
        - 4: iterative process initialized by the least squares method
  \var  cs_space_disc_t::anomax
        <a name="anomax"></a>
        non orthogonality angle of the faces, in radians.
        For larger angle values, cells with one node on the wall are kept in the
        extended support of the neighboring cells.
  \var  cs_space_disc_t::iflxmw
        method to compute interior mass flux due to ALE mesh velocity
        - 1: based on cell center mesh velocity
        - 0: based on nodes displacement

*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_var_cal_opt_t

  \brief structure containing the variable calculation options.

  \var  cs_var_cal_opt_t::iwarni
        \anchor iwarni
        \ref iwarni characterises the level of detail of the outputs for a
        variable. The quantity of information increases with its value.
        Impose the value 0 or 1 for a reasonable listing size. Impose the
        value 2 to get a maximum quantity of information, in case of problem
        during the execution.

 \var  cs_var_cal_opt_t::iconv
        \anchor iconv
        For each unknown variable to calculate, indicates if the convection is
        taken into account (1) or not (0). By default,
        \ref cs_var_cal_opt_t::iconv "iconv" is set to 0 for the pressure
        (variable \ref ipr) or f in v2f modelling (variable \ref ifb)
        and set to 1 for the other unknowns.

  \var  cs_var_cal_opt_t::istat
        \anchor istat
        For each unknown variable to calculate, indicates whether unsteady
        terms are present (1) or not (0) in the matrices. By default,
        \ref cs_var_cal_opt_t::istat "istat" is set to 0 for the pressure
        (variable \ref ipr) or f in v2f modelling (variable \ref ifb) and
        set to 1 for the other unknowns.

  \var  cs_var_cal_opt_t::idiff
        \anchor idiff
        For each unknown variable to calculate, indicates whether the
        diffusion is taken into account (1) or not (0).

  \var  cs_var_cal_opt_t::idifft
        \anchor idifft
        For each unknown variable to calculate, when diffusion is taken into
        account (\ref idiff = 1), \ref idifft indicates if the turbulent
        diffusion is taken into account (\ref idifft = 1) or not (0).

  \var  cs_var_cal_opt_t::idften
        \anchor idften
        Type of diffusivity:
        - 1: scalar diffusivity
        - 3: orthotropic diffusivity
        - 6: symmetric tensor diffusivity

  \var  cs_var_cal_opt_t::iswdyn
        \anchor iswdyn
        Dynamic relaxation type:
        - 0 no dynamic relaxation
        - 1 dynamic relaxation depending on \f$ \delta \varia^k \f$
        - 2 dynamic relaxation depending on \f$ \delta \varia^k \f$ and
        \f$ \delta \varia^{k-1} \f$.

  \var  cs_var_cal_opt_t::ischcv
        \anchor ischcv
        For each unknown variable to calculate, \ref ischcv indicates the type of
        second-order convective scheme
        - 0: Second Order Linear Upwind
        - 1: Centered \n
        - 2: Second Order with upwind-gradient reconstruction (SOLU) \n
        Useful for all the unknowns variables which are convected
        (\ref iconv = 1) and for which a second-order scheme is used
        (\ref blencv > 0).

  \var  cs_var_cal_opt_t::ibdtso
        \anchor ibdtso
        Backward differential scheme in time order.

  \var  cs_var_cal_opt_t::isstpc
        \anchor isstpc
        For each unknown variable to calculate, isstpc indicates whether a slope
        test should be used to switch from a second-order to an upwind convective
        scheme under certain conditions, to ensure stability.
        - -1: deprecated slope test activated for the considered unknown
              (for vector variable only)
        - 0: slope test activated for the considered unknown
        - 1: slope test deactivated for the considered unknown
        - 2: continuous limiter ensuring boundedness (beta limiter)
        - 3: NVD/TVD Scheme
             Then "limiter_choice" keyword must be set:
             * 0: Gamma
             * 1: SMART
             * 2: CUBISTA
             * 3: SUPERBEE
             * 4: MUSCL
             * 5: MINMOD
             * 6: CLAM
             * 7: STOIC
             * 8: OSHER
             * 9: WASEB
             * --- VOF scheme ---
             * 10: M-HRIC
             * 11: M-CICSAM
        Useful for all the unknowns variable which are convected (\ref iconv = 1)
        and for which a second-order scheme is used (\ref blencv > 0).
        The use of the slope test stabilises the calculation but may bring
        the order in space to decrease quickly.

  \var  cs_var_cal_opt_t::nswrgr
        \anchor nswrgr
        For each unknown variable, \ref nswrgr <= 1 indicates that the gradients
        are not reconstructed
         - if \ref imrgra = 0 or 4, \ref nswrgr is the number of iterations for
         the gradient reconstruction
         - if \ref imrgra = 1, 2 or 3, \ref nswrgr > 1 indicates that the
         gradients are reconstructed (but the method is not iterative, so any
         value larger than 1 for \ref nswrgr yields the same result).\n

  \var  cs_var_cal_opt_t::nswrsm
        \anchor nswrsm
        For each unknown variable, nswrsm indicates the number of iterations for
        the reconstruction of the right-hand sides of the equations with a
        first-order scheme in time (standard case), the default values are 2 for
        pressure and 1 for the other variables. With a second-order scheme in
        time (\ref optcal::ischtp "ischtp" = 2) or LES, the default values are
        5 for pressure and 10 for the other variables.

  \var  cs_var_cal_opt_t::imrgra
        \anchor imrgra
        Indicates the type of gradient reconstruction (one method for all the
        variables)
           - 0: iterative reconstruction of the non-orthogonalities
           - 1: least squares method based on the first neighbour cells (cells
        which share a face with the treated cell)
           - 2: least squares method based on the extended neighbourhood (cells
        which share a node with the treated cell)
           - 3: least squares method based on a partial extended neighbourhood
        (all first neighbours plus the extended neighbourhood cells that are
        connected to a face where the non-orthogonality angle is larger than
        parameter anomax)
           - 4: iterative reconstruction with initialisation using the least
        squares method (first neighbours)
           - 5: iterative reconstruction with initialisation using the least
        squares method based on an extended neighbourhood
           - 6: iterative reconstruction with initialisation using the least
        squares method based on a partial extended neighbourhood
        if \ref imrgra fails due to probable mesh quality problems, it is usually
        effective to use \ref imrgra = 3. Moreover, \ref imrgra = 3 is usually
        faster than \ref imrgra = 0 (but with less feedback on its use).
        It should be noted that \ref imrgra = 1, 2 or 3 automatically triggers
        a gradient limitation procedure. See \ref imligr.\n
        Useful if and only if there is \ref nswrgr > 1 for at least one variable.
        Also, pressure gradients (or other gradients deriving from a potential)
        always use an iterative reconstruction. To force a non-iterative
        reconstruction for those gradients, a negative value of this keyword
        may be used, in which case the method matching the absolute value
        of the keyword will be used.

  \var  cs_var_cal_opt_t::imligr
        \anchor imligr
        For each unknown variable, indicates the type of gradient limitation
           - -1: no limitation
           - 0: based on the neighbours
           - 1: superior order\n
        For all the unknowns, \ref imligr is initialized to -1 if \ref imrgra
        = 0 or 4 and to 1 if \ref imrgra = 1, 2 or 3.

  \var  cs_var_cal_opt_t::ircflu
        \anchor ircflu
        For each unknown variable, \ref ircflu indicates whether the convective
        and diffusive fluxes at the faces should be reconstructed:
           - 0: no reconstruction
           - 1: reconstruction \n
        Deactivating the reconstruction of the fluxes can have a stabilising
        effect on the calculation. It is sometimes useful with the
        \f$ k-\epsilon \f$ model, if the mesh is strongly non-orthogonal
        in the near-wall region, where the gradients of k and \f$ \epsilon \f$
        are strong. In such a case, setting \ref ircflu = 0 will probably help
        (switching to a first order convective scheme, \ref blencv = 0, for k
        and \f$ \epsilon \f$ might also help in that case).

  \var  cs_var_cal_opt_t::iwgrec
        \anchor iwgrec
        Gradient calculation
          - 0: standard
          - 1: weighted

  \var  cs_var_cal_opt_t::thetav
        \anchor thetav
        For each variable variable, thetav is the value of \f$ \theta \f$ used to
        express at the second-order the terms of convection, diffusion and the
        source terms which are linear functions of the solved variable (according
        to the formula \f$ \phi^{n+\theta}
        = (1-\theta) \phi^n + \theta \phi^{n+1}\f$.
        Generally, only the values 1 and 0.5 are used. The user is not allowed to
        modify this variable.
           - 1: first-order
           - 0.5: second-order \n
        Concerning the pressure, the value of \ref thetav is always 1. Concerning
        the other variables, the value \ref thetav = 0.5 is used when the
        second-order time scheme is activated by \ref ischtp = 2 (standard value
        for LES calculations), otherwise \ref thetav is set to 1.

  \var  cs_var_cal_opt_t::blencv
        \anchor blencv
        For each unknown variable to calculate, blencv indicates the proportion
        of second-order convective scheme (0 corresponds to an upwind first-order
        scheme); in case of LES calculation, a second-order scheme is
        recommended and activated by default (\ref blencv = 1).\n
        Useful for all the unknowns variable for which \ref iconv = 1.

  \var  cs_var_cal_opt_t::blend_st
        \anchor blend_st
        For each unknown variable to calculate, blend_st indicates the proportion
        of second-order convective scheme (0 corresponds to an upwind first-order
        scheme) after the slope test is activated;
        in case of LES calculation, a second-order scheme is
        recommended and activated by default (\ref blend_st = 1).\n
        Useful for all the unknowns variable for which \ref iconv = 1.


  \var  cs_var_cal_opt_t::epsilo
        \anchor epsilo
        <a name="epsilo"></a>
        For each unknown variable, relative precision for the solution of the
        linear system. The default value is \ref epsilo = \f$ 10^-8 \f$ . This
        value is set low on purpose. When there are enough iterations on the
        reconstruction of the right-hand side of the equation, the value may be
        increased (by default, in case of second-order in time, with
        \ref nswrsm = 5 or 10, \ref epsilo is increased  to \f$ 10^-5 \f$.

  \var  cs_var_cal_opt_t::epsrsm
        \anchor epsrsm
        For each unknown variable, relative precision on the reconstruction of
        the right hand-side. The default value is \ref epsrsm = \f$ 10^-8 \f$.
        This value is set low on purpose. When there are not enough iterations on
        the reconstruction of the right-hand side of the equation, the value may
        be increased (by default, in case of second-order in time, with
        \ref nswrsm = 5 or 10, \ref epsrsm is increased to
        \f$ 10^-5 \f$ ).

  \var  cs_var_cal_opt_t::epsrgr
        \anchor epsrgr
        For each unknown variable, relative precision for the iterative gradient
        reconstruction.\n Useful for all the unknowns when \ref imrgra = 0 or 4.

  \var  cs_var_cal_opt_t::climgr
        \anchor climgr
        For each unknown variable, factor of gradient limitation (high value means
        little limitation). \n
        Useful for all the unknowns variables for which \ref imligr = -1.

  \var  cs_var_cal_opt_t::extrag
        \anchor extrag
        For the variable pressure \ref ipr, extrapolation coefficient of the
        gradients at the boundaries. It affects only the Neumann conditions.
        The only possible values of \ref extrag are:
             - 0: homogeneous Neumann calculated at first-order
             - 0.5: improved homogeneous Neumann, calculated at second-order in
        the case of an orthogonal mesh and at first-order otherwise
             - 1: gradient extrapolation (gradient at the boundary face equal to
        the gradient in the neighbour cell), calculated at second-order in the
        case of an orthogonal mesh and at first-order otherwise extrag often
        allows to correct the non-physical velocities that appear on horizontal
        walls when density is variable and there is gravity. It is strongly
        advised to keep \ref extrag = 0 for the variables apart from pressure.
        See also \ref cs_stokes_model_t::iphydr "iphydr". In practice, only the
        values 0 and 1 are allowed. The value 0.5 is not allowed by default (but
        the lock can be overridden if necessary, contact the development team).

  \var  cs_var_cal_opt_t::relaxv
        \anchor relaxv
        For each variable ivar, relaxation coefficient of the variable. This
        relaxation parameter is only useful for the pressure with the unsteady
        algorithm (so as to improve the convergence in case of meshes of
        insufficient quality or and for some of the turbulent models
        (\ref iturb = 20, 21, 50 or 60 and \ref optcal::ikecou "ikecou" = 0;
        if \ref optcal::ikecou "ikecou" = 1, \ref relaxv is not used, whatever
        its value may be).
        Default values are 0.7 for turbulent variables and 1. for pressure.
        \ref relaxv also stores the value of the relaxation coefficient when
        using the steady algorithm, deduced from the value of
        \ref optcal::relxst "relxst" (defaulting to \ref relaxv = 1.
        - \ref optcal::relxst "relxst"). Useful only for the pressure and
        for turbulent variables if and only if (\f$ k-\epsilon \f$, v2f
        or \f$ k-\omega \f$ models without coupling) with the unsteady
        algorithm. Always useful with the steady algorithm.
*/

/*----------------------------------------------------------------------------*/

/*!
  \struct cs_piso_t

  \brief PISO options descriptor.

  Members of the PISO structure are publicly accessible, to allow for
  concise  syntax, as they are expected to be used in many places.

  \var  cs_piso_t::nterup
        number of interations on the pressure-velocity coupling on Navier-Stokes
  \var  cs_piso_t::epsup
        relative precision for the convergence test of the iterative process on
        pressure-velocity coupling
  \var  cs_piso_t::xnrmu
        norm  of the increment \f$ \vect{u}^{k+1} - \vect{u}^k \f$ of the
        iterative process on pressure-velocity coupling
  \var  cs_piso_t::xnrmu0
        norm of \f$ \vect{u}^0 \f$
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

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
  int       location_id;        /* Propert location id */

} cs_user_property_def_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default variable compute options */

static cs_var_cal_opt_t _var_cal_opt =
{
  .iwarni = 0,
  .iconv  = 1,
  .istat  = 1,
  .idiff  = 1,
  .idifft = 1,
  .idften = CS_ISOTROPIC_DIFFUSION,
  .iswdyn = 0,
  .ischcv = 1,
  .ibdtso = 1,
  .isstpc = 1,
  .nswrgr = 100,
  .nswrsm = 1,
  .imrgra = 0,
  .imligr = -1,
  .ircflu = 1,
  .iwgrec = 0,
  .icoupl = -1,
  .thetav = 1.,
  .blencv = 1.,
  .blend_st = 0.,
  .epsilo = 1.e-8,
  .epsrsm = 1.e-7,
  .epsrgr = 1.e-5,
  .climgr = 1.5,
  .extrag = 0.,
  .relaxv = 1.
};

/* Space discretisation options structure and associated pointer */

static cs_space_disc_t  _space_disc =
{
  .imvisf = 0,
  .imrgra = 0,
  .anomax = -1e12*10.,
  .iflxmw = 1
};

const cs_space_disc_t  *cs_glob_space_disc = &_space_disc;

/* PISO structure and associated pointer */

static cs_piso_t  _piso =
{
  .nterup = 1,
  .epsup = 1e-5,
  .xnrmu = 0.,
  .xnrmu0 = 0.,
  .n_buoyant_scal = 0
};

const cs_piso_t  *cs_glob_piso = &_piso;

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

static cs_gas_mix_species_prop_t _gas_mix_species_prop =
{
  -1.,   /* molar mass                              */
  -1.,   /* specific heat                           */
  -1.,   /* volume diffusion                        */
  -1.,   /* dynamic viscosity a                     */
  -1.,   /* dynamic viscosity b                     */
  -1.,   /* thermal conductivity a                  */
  -1.,   /* thermal conductivity b                  */
  -1.,   /* reference viscosity (Sutherland)        */
  -1.,   /* reference conductivity (Sutherland)     */
  -1.,   /* reference temperature for viscosity     */
  -1.,   /* reference temperature for conductivity  */
  -1.,   /* Sutherland temperature for viscosity    */
  -1.,   /* Sutherland temperature for conductivity */
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
                             double  **anomax,
                             int     **iflxmw);

void
cs_f_piso_get_pointers(int     **nterup,
                       double  **epsup,
                       double  **xnrmu,
                       double  **xnrmu0,
                       int     **n_buoyant_scal);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Log values of the structure */

static void
_log_func_var_cal_opt(const void *t)
{
  const char fmt_i[] = N_("      %-19s  %-4d\n");
  const char fmt_r[] = N_("      %-19s  %-12.3g\n");
  const cs_var_cal_opt_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iwarni", _t->iwarni);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iconv ", _t->iconv );
  cs_log_printf(CS_LOG_SETUP, fmt_i, "istat ", _t->istat );
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idiff ", _t->idiff );
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idifft", _t->idifft);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "idften", _t->idften);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iswdyn", _t->iswdyn);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ischcv", _t->ischcv);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ibdtso", _t->ibdtso);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "isstpc", _t->isstpc);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "nswrgr", _t->nswrgr);
  cs_log_printf(CS_LOG_SETUP, fmt_i, "nswrsm", _t->nswrsm);
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
  cs_log_printf(CS_LOG_SETUP, fmt_r, "extrag", _t->extrag);
  cs_log_printf(CS_LOG_SETUP, fmt_r, "relaxv", _t->relaxv);
}

/* Log default values of the structure */

static void
_log_func_default_var_cal_opt(const void *t)
{
  const char fmt_i[] = "      %-19s  %-12d %s\n";
  const char fmt_r[] = "      %-19s  %-12.3g %s\n";
  const cs_var_cal_opt_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP,"  var_cal_opt\n");

  cs_log_printf(CS_LOG_SETUP,_("    Printing\n"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "iwarni", _t->iwarni,
                _("Verbosity level: 0, 1 or 2"));

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
                _("Type of diffusivity: scalar (1), orthotropic (3) or symmetric "
                  "tensor (6)"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "ischcv", _t->ischcv,
                _("Type of convective scheme: 2nd order with centered-gradient "
                  "upwind reconstruction (0), centered (1), "
                  "2nd order with upwind-gradient upwind-reconstruction (SOLU) (2)"));
  cs_log_printf(CS_LOG_SETUP, fmt_i, "isstpc", _t->isstpc,
                _("0 for slope test, 1 for no slope test, 2 for min/max limiter "
                  "and 3 for NVD/TVD scheme"));
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
  cs_log_printf(CS_LOG_SETUP, fmt_r, "extrag", _t->extrag,
                _("[0.;1.] (gradients extrapolation)"));
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

/* Log values of the structure */

static void
_log_func_gas_mix_species_prop(const void *t)
{
  const char fmt[] = N_("      %-19s  %-12.3g\n");
  const cs_gas_mix_species_prop_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, fmt, "mol_mas ", _t->mol_mas);
  cs_log_printf(CS_LOG_SETUP, fmt, "cp      ", _t->cp);
  cs_log_printf(CS_LOG_SETUP, fmt, "vol_diff", _t->vol_dif);
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_a    ", _t->mu_a);
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_b    ", _t->mu_b);
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_a", _t->lambda_a);
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_b", _t->lambda_b);
  cs_log_printf(CS_LOG_SETUP, fmt, "muref   ", _t->muref);
  cs_log_printf(CS_LOG_SETUP, fmt, "lamref  ", _t->lamref);
  cs_log_printf(CS_LOG_SETUP, fmt, "trefmu  ", _t->trefmu);
  cs_log_printf(CS_LOG_SETUP, fmt, "treflam ", _t->treflam);
  cs_log_printf(CS_LOG_SETUP, fmt, "smu     ", _t->smu);
  cs_log_printf(CS_LOG_SETUP, fmt, "slam    ", _t->slam);
}

/* Log default values of the structure */

static void
_log_func_default_gas_mix_species_prop(const void *t)
{
  const char fmt[] = "      %-19s  %-12.3g %s\n";
  const cs_gas_mix_species_prop_t *_t = (const void *)t;
  cs_log_printf(CS_LOG_SETUP, fmt, "mol_mas ", _t->mol_mas,
                _("Molar mass"));
  cs_log_printf(CS_LOG_SETUP, fmt, "cp      ", _t->cp,
                _("Specific heat"));
  cs_log_printf(CS_LOG_SETUP, fmt, "vol_diff", _t->vol_dif,
                _("Volume diffusion"));
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_a    ", _t->mu_a,
                _("Dynamic viscosity a"));
  cs_log_printf(CS_LOG_SETUP, fmt, "mu_b    ", _t->mu_b,
                _("Dynamic viscosity b"));
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_a", _t->lambda_a,
                _("Thermal conductivity a"));
  cs_log_printf(CS_LOG_SETUP, fmt, "lambda_b", _t->lambda_b,
                _("Thermal conductivity b"));
  cs_log_printf(CS_LOG_SETUP, fmt, "muref   ", _t->muref,
                _("Reference thermal viscosity (Sutherland)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "lamref  ", _t->lamref,
                _("Reference thermal conductivity (Sutherland)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "trefmu  ", _t->trefmu,
                _("Reference temperature (Sutherland for viscosity)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "treflam ", _t->treflam,
                _("Reference temperature (Sutherland conductivity)"));
  cs_log_printf(CS_LOG_SETUP, fmt, "smu     ", _t->smu,
                _("Sutherland temperature for viscosity"));
  cs_log_printf(CS_LOG_SETUP, fmt, "slam    ", _t->slam,
                _("Sutherland temperature for conductivity"));
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global space disc structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   imvisf  --> pointer to cs_glob_space_disc->imvisf
 *   imrgra  --> pointer to cs_glob_space_disc->imrgra
 *   anomax  --> pointer to cs_glob_space_disc->anomax
 *   iflxmw  --> pointer to cs_glob_space_disc->iflxmw
 *----------------------------------------------------------------------------*/

void
cs_f_space_disc_get_pointers(int     **imvisf,
                             int     **imrgra,
                             double  **anomax,
                             int     **iflxmw)
{
  *imvisf = &(_space_disc.imvisf);
  *imrgra = &(_space_disc.imrgra);
  *anomax = &(_space_disc.anomax);
  *iflxmw = &(_space_disc.iflxmw);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global piso structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   nterup          --> pointer to cs_glob_piso->nterup
 *   epsup           --> pointer to cs_glob_piso->epsup
 *   xnrmu           --> pointer to cs_glob_piso->xnrmu
 *   xnrmu0          --> pointer to cs_glob_piso->xnrmu0
 *   n_buoyant_scal  --> pointer to cs_glob_piso->n_buoyant_scal
 *----------------------------------------------------------------------------*/

void
cs_f_piso_get_pointers(int     **nterup,
                       double  **epsup,
                       double  **xnrmu,
                       double  **xnrmu0,
                       int     **n_buoyant_scal)
{
  *nterup = &(_piso.nterup);
  *epsup  = &(_piso.epsup);
  *xnrmu  = &(_piso.xnrmu);
  *xnrmu0 = &(_piso.xnrmu0);
  *n_buoyant_scal = &(_piso.n_buoyant_scal);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_space_disc.
 *
 * Needed to initialize structure with GUI and user C functions.
 *
 * \return  piso information structure
 */
/*----------------------------------------------------------------------------*/

cs_space_disc_t *
cs_get_glob_space_disc(void)
{
  return &_space_disc;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_piso.
 *
 * Needed to initialize structure with GUI.
 *
 * \return  piso information structure
 */
/*----------------------------------------------------------------------------*/

cs_piso_t *
cs_get_glob_piso(void)
{
  return &_piso;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Count and set number of buoyant scalars.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_set_n_buoyant_scalars(void)
{
  const int n_fields = cs_field_n_fields()
  const int key_sca = cs_field_key_id("scalar_id");
  const int key_buo = cs_field_key_id("is_buoyant");

  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (   f->type & CS_FIELD_VARIABLE
        && cs_field_get_key_int(f, key_sca) > -1) {
      if (cs_field_get_key_int(f, key_buo)) {
        _piso.n_buoyant_scal += 1;
      }
    }
  }
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
  cs_field_define_key_int("inner_mass_flux_id", -1, 0);
  cs_field_define_key_int("boundary_mass_flux_id", -1, 0);

  cs_field_define_key_int("variable_id", -1, 0); /* inverse of ivarfl(ivar) */
  cs_field_define_key_int("scalar_id", -1, 0);   /* inverse of isca(iscal) */

  cs_field_define_key_int("scalar_diffusivity_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_double("scalar_diffusivity_ref",
                             -1.e12*10., CS_FIELD_VARIABLE); /* visls0(iscal) */

  cs_field_define_key_int("scalar_density_id", -1, CS_FIELD_VARIABLE);

  /* is the field buoyant? -1 if not, 1 if yes */
  cs_field_define_key_int("is_buoyant", 0, CS_FIELD_VARIABLE);

  cs_field_define_key_int("turbulent_flux_model", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("turbulent_flux_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_double("turbulent_schmidt",
                             -1.e12*10., CS_FIELD_VARIABLE); /* sigmas(iscal) */

  cs_field_define_key_int("gradient_weighting_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("diffusivity_tensor", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("drift_scalar_model", 0, 0);

  cs_field_define_key_int("scalar_class", 0, 0);
  cs_field_define_key_int("first_moment_id", -1, 0); /* iscavr(iscal) */

  cs_field_define_key_int("syrthes_coupling", 0, 0); /* icpsyr(iscal) */

  cs_field_define_key_int("source_term_prev_id", -1, CS_FIELD_VARIABLE);
  /* TODO merge with previous key word */
  cs_field_define_key_int("source_term_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("slope_test_upwind_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("clipping_id", -1, CS_FIELD_VARIABLE);

  cs_field_define_key_int("boundary_value_id", -1, 0);

  cs_field_define_key_int("convection_limiter_id", -1, 0);

  cs_field_define_key_int("coupling_entity", -1, 0);

  cs_field_define_key_double("min_scalar_clipping", -1.e12, 0);
  cs_field_define_key_double("max_scalar_clipping", 1.e12, 0);

  cs_field_define_key_int("measures_set_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("opt_interp_id", -1, CS_FIELD_VARIABLE);
  cs_field_define_key_int("opt_interp_analysis_id", -1, CS_FIELD_VARIABLE);

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
  cs_field_define_key_struct("var_cal_opt",
                             &_var_cal_opt,
                             _log_func_var_cal_opt,
                             _log_func_default_var_cal_opt,
                             sizeof(cs_var_cal_opt_t),
                             CS_FIELD_VARIABLE);

  /* Structure containing the solving info of the field variables
     (used for listing, not setup, so set NULL setup logging function) */
  cs_field_define_key_struct("solving_info",
                             &_solving_info,
                             NULL,
                             NULL,
                             sizeof(cs_solving_info_t),
                             CS_FIELD_VARIABLE);
  cs_field_key_disable_setup_log(cs_field_key_id("solving_info"));

  /* Restart options */
  cs_field_define_key_int("restart_file", CS_RESTART_DISABLED, 0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define field key for condensation.
 *
 * Note: this should be moved in the future to a condensation-specific file.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_define_field_key_gas_mix(void)
{
  /* Structure containing physical properties relative to
     species scalars used by the gas mixture modelling */
  cs_field_define_key_struct("gas_mix_species_prop",
                             &_gas_mix_species_prop,
                             _log_func_gas_mix_species_prop,
                             _log_func_default_gas_mix_species_prop,
                             sizeof(cs_gas_mix_species_prop_t),
                             0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read general restart info.
 *
 * This updates the previous time step info.
 */
/*----------------------------------------------------------------------------*/

void
cs_parameters_read_restart_info(void)
{
  if (cs_restart_present()) {
    cs_restart_t *r
      = cs_restart_create("main", "restart", CS_RESTART_MODE_READ);
    cs_restart_read_time_step_info(r);
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

  if (dim != 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Only user variables of dimension 1 are currently handled,\n"
                "but %s is defined with dimension %d."),
              name, dim);

  _n_user_variables++;
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
          bft_error(__FILE__, __LINE__, 0,
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
 * \brief Return a local variable calculation options structure,
 *        with default options.
 *
 * \return  variable calculations options structure
 */
/*----------------------------------------------------------------------------*/

cs_var_cal_opt_t
cs_parameters_var_cal_opt_default(void)
{
  return _var_cal_opt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the space discretization structure to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_space_disc_log_setup(void)
{
   cs_log_printf
     (CS_LOG_SETUP,
      _("\n"
        "Space discretization options\n"
        "----------------------------\n\n"
        "    imvisf:      %d (face interpolation\n"
        "                    0: arithmetic\n"
        "                    1: harmonic)\n"
        "\n"
        "    imrgra:      %d (type of gradient reconstruction\n"
        "                    0: iterative process\n"
        "                    1: standard least squares method\n"
        "                    2: least squares method with extended "
        "neighborhood\n"
        "                    3: standard least squares method with reduced "
        "extended neighborhood\n"
        "                    4: iterative process initialized by the least "
        "squares method)\n"
        "\n"
        "    anomax       %-12.3g (non-orthogonality angle (rad) above which "
        "cells are\n"
        "                    selected for the extended neighborhood)\n"
        "    iflxmw:      %d (method to compute inner mass flux due to mesh "
        "velocity in ALE\n"
        "                    0: based on mesh velocity at cell centers\n"
        "                    1: based on nodes displacement)\n"),
        cs_glob_space_disc->imvisf,
        cs_glob_space_disc->imrgra,
        cs_glob_space_disc->anomax,
        cs_glob_space_disc->iflxmw);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

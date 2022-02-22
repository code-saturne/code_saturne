/*============================================================================
 * Base time step data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_base.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_time_step.c
        base time step data.

  \enum cs_time_step_type_t

  \brief Time step type

  \var CS_TIME_STEP_STEADY
       Steady (SIMPLE) algorithm
  \var CS_TIME_STEP_CONSTANT
       Unsteady (SIMPLEC) algorithm with constant time step
  \var CS_TIME_STEP_ADAPTIVE
       Unsteady (SIMPLEC) algorithm with time-adaptive (varying) time step
  \var CS_TIME_STEP_LOCAL
       Pseudo-steady with time-and-space varying time step (SIMPLEC)

  \struct cs_time_step_t

  \brief time step descriptor

  Members of this time step are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_time_step_t::is_variable
        0 if time step is fixed in time, 1 otherwise
  \var  cs_time_step_t::is_local
        0 if time step is uniform in space, 0 if is is local
        (in which case the time value is only a reference)
  \var  cs_time_step_t::nt_prev
        Absolute time step number for previous calculation.
        In the case of a restart calculation, \ref nt_prev
        is read from the restart file. Otherwise, it is
        initialised to 0 \ref nt_prev is initialised
        automatically by the code, its value is not to be
        modified by the user.
  \var  cs_time_step_t::nt_cur
        current absolute time step number
        In case of restart, this is equal to \ref nt_prev + number
        of new iterations.
  \var  cs_time_step_t::nt_max
        maximum absolute time step number
        For the restart calculations, \ref nt_max takes into
        account the number of time steps of the previous calculations.
        For instance, after a first calculation of 3 time steps, a
        restart file of 2 time steps is realised by setting
        \ref nt_max = 3+2 = 5
  \var  cs_time_step_t::nt_ini
        number of time step for initialization
  \var  cs_time_step_t::t_prev
        Absolute time value for previous calculation.
        In the case of a restart calculation, \ref t_prev is read from
        the restart file. Otherwise it is initialised to 0.\n
        \ref t_prev is initialised automatically by the code,
        its value is not to be modified by the user.
  \var  cs_time_step_t::t_cur
        Current absolute time.
        For the restart calculations, \ref t_cur takes
        into account the physical time of the previous calculations.\n
        If the time step is uniform (\ref cs_time_step_options_t::idtvar "idtvar"
        = CS_TIME_STEP_CONSTANT or CS_ADAPTATIVE_TIME_STEP),
        \ref t_cur increases of \ref dt (value of the time step)
        at each iteration.
        If the time step is non-uniform (\ref cs_time_step_options_t::idtvar
        "idtvar"=CS_TIME_STEP_LOCAL),
        \ref t_cur increases of \ref cs_time_step_t::dt_ref "dt_ref" at each
        time step.\n
        \ref t_cur is initialised and updated automatically by the code,
        its value is not to be modified by the user.
  \var  cs_time_step_t::t_max
        maximum absolute time

  \struct cs_time_step_options_t

  \brief time step options descriptor

  Members of this time step options descriptor are publicly accessible, to
  allow for concise syntax.

  \var  cs_time_step_options_t::iptlro
        Clip the time step with respect to the buoyant effects\n

        When density gradients and gravity are present, a local thermal time
        step can be calculated, based on the Brunt-Vaisala frequency. In
        numerical simulations, it is usually wise for the time step to be
        lower than this limit, otherwise numerical instabilities may appear.\n
        \ref iptlro indicates whether the time step should be limited to the
        local thermal time step (=1) or not (=0).\n
        When \ref iptlro=1, the log shows the number of cells where the
        time step has been clipped due to the thermal criterion, as well as
        the maximum ratio between the time step and the maximum thermal time
        step. If \ref idtvar=CS_TIME_STEP_CONSTANT, since the time step is fixed
        and cannot be clipped, this ratio can be greater than 1.
        When \ref idtvar > CS_TIME_STEP_CONSTANT, this
        ratio will be less than 1, except if the constraint \ref dtmin has
        prevented the code from reaching a sufficiently low value for \ref dt.
        Useful when density gradients and gravity are present.

  \var  cs_time_step_options_t::idtvar
        Time step type: constant, adaptive, or steady algorithm variants.
        If the numerical scheme is a second-order in time, only the
        option CS_TIME_STEP_CONSTANT is allowed.
        \anchor idtvar

  \var  cs_time_step_t::dt_ref
        Reference time step.\n

        This is the time step value used in the case of a calculation run with a
        uniform and constant time step, i.e. \ref idtvar =CS_TIME_STEP_CONSTANT
        (restart calculation or not). It is the value used to initialize the
        time step in the case of an initial calculation run with a non-constant
        time step (\ref idtvar=CS_TIME_STEP_ADAPTATIVE or CS_TIME_STEP_LOCAL).
        It is also the value used to initialise the time step in the case of
        a restart calculation in which the type of time step has been changed
        (for instance, \ref idtvar=CS_TIME_STEP_ADAPTATIVE in the new calculation
        and \ref idtvar = CS_TIME_STEP_CONSTANT or CS_STEP_LOCAL_TIME in the
        previous calculation).\n
        See \ref user_initialization_time_step for examples.

  \var  cs_time_step_options_t::coumax
        Maximum Courant number (when \ref idtvar is different
        from CS_TIME_STEP_CONSTANT).

  \var  cs_time_step_options_t::cflmmx
        Max. Courant number for the continuity equation in compressible model.

  \var  cs_time_step_options_t::foumax
        Maximum Fourier number (when \ref idtvar is different from
        CS_TIME_STEP_CONSTANT).

  \var  cs_time_step_options_t::varrdt
        Maximum allowed relative increase in the calculated time step value
        between two successive time steps (to ensure stability, any decrease
        in the time step is immediate and without limit).\n
        Useful when idtvar is different from CS_TIME_STEP_CONSTANT.

  \var  cs_time_step_options_t::dtmin
        Lower limit for the calculated time step when idtvar is different from
        CS_TIME_STEP_CONSTANT.
        Take \ref dtmin = min (ld/ud, sqrt(lt/(gdelta rho/rho)), ...).

  \var  cs_time_step_options_t::dtmax
        Upper limit for the calculated time step when idtvar is different from
        CS_TIME_STEP_CONSTANT.
        Take \ref dtmax = max (ld/ud, sqrt(lt/(gdelta rho/rho)), ...).

  \var  cs_time_step_options_t::relxst
        Relaxation coefficient for the steady algorithm.
        \ref relxst = 1 : no relaxation.

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

static const char *cs_time_step_type_enum_name[]
  = {"CS_TIME_STEP_STEADY",
     "CS_TIME_STEP_CONSTANT",
     "CS_TIME_STEP_ADAPTIVE",
     "CS_TIME_STEP_LOCAL"};

static const char *cs_time_step_type_name[]
  = {N_("steady algorithm"),
     N_("constant"),
     N_("time-adaptive"),
     N_("local in time and space (pseudo-steady)")};

/* main time step structure and associated pointer */

static cs_time_step_t  _time_step = {
  .is_variable = 0,
  .is_local = 0,
  .nt_prev = 0,
  .nt_cur = 0,
  .nt_max = -1,
  .nt_ini = 2,
  .t_prev = 0.,
  .t_cur = 0.,
  .t_max = -1.,
  .dt = {0.1, 0.1, 0.1},
  .dt_ref = 0.1,
  .dt_next = 0.1
};

static cs_time_step_options_t  _time_step_options = {
  .iptlro = 0,
  .idtvar = CS_TIME_STEP_CONSTANT, /* constant time step by default */
  .coumax = 1.,
  .cflmmx = 0.99,
  .foumax = 10.,
  .varrdt = 0.1,
  .dtmin  = -1.e13,
  .dtmax  = -1.e13,
  .relxst = 0.7 /* Not used in CDO schemes */
};

const cs_time_step_t  *cs_glob_time_step = &_time_step;

const cs_time_step_options_t  *cs_glob_time_step_options = &_time_step_options;

static double _c = 0; /* compensation term for Kahan sum */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_time_step_get_pointers(int     **nt_prev,
                            int     **nt_cur,
                            int     **nt_max,
                            int     **nt_ini,
                            double  **dtref,
                            double  **t_prev,
                            double  **t_cur,
                            double  **t_max);

void
cs_f_time_step_options_get_pointers(int    **iptlro,
                                    int    **idtvar,
                                    double **coumax,
                                    double **cflmmx,
                                    double **foumax,
                                    double **varrdt,
                                    double **dtmin,
                                    double **dtmax,
                                    double **relxst);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global time step structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   nt_prev --> pointer to cs_glob_time_step->nt_prev
 *   nt_cur  --> pointer to cs_glob_time_step->nt_cur
 *   nt_max  --> pointer to cs_glob_time_step->nt_max
 *   nt_ini  --> pointer to cs_glob_time_step->nt_ini
 *   dt_ref  --> pointer to cs_glob_time_step->dt_ref
 *   t_prev  --> pointer to cs_glob_time_step->t_prev
 *   t_cur   --> pointer to cs_glob_time_step->t_cur
 *   t_max   --> pointer to cs_glob_time_step->t_ax
 *----------------------------------------------------------------------------*/

void
cs_f_time_step_get_pointers(int      **nt_prev,
                            int      **nt_cur,
                            int      **nt_max,
                            int      **nt_ini,
                            double   **dtref,
                            double   **t_prev,
                            double   **t_cur,
                            double   **t_max)
{
  *nt_prev = &(_time_step.nt_prev);
  *nt_cur = &(_time_step.nt_cur);
  *nt_max = &(_time_step.nt_max);
  *nt_ini = &(_time_step.nt_ini);
  *dtref  = &(_time_step.dt_ref);
  *t_prev = &(_time_step.t_prev);
  *t_cur = &(_time_step.t_cur);
  *t_max = &(_time_step.t_max);
}

/*----------------------------------------------------------------------------
 * Get pointers to members of the global time step structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   iptlro --> pointer to cs_glob_time_step_options->iptlro
 *   idtvar --> pointer to cs_glob_time_step_options->idtvar
 *   coumax --> pointer to cs_glob_time_step_options->coumax
 *   cflmmx --> pointer to cs_glob_time_step_options->cflmmx
 *   foumax --> pointer to cs_glob_time_step_options->foumax
 *   varrdt --> pointer to cs_glob_time_step_options->varrdt
 *   dtmin  --> pointer to cs_glob_time_step_options->dtmin
 *   dtmax  --> pointer to cs_glob_time_step_options->dtmax
 *   relxst --> pointer to cs_glob_time_step_options->relxst
 *----------------------------------------------------------------------------*/

void
cs_f_time_step_options_get_pointers(int    **iptlro,
                                    int    **idtvar,
                                    double **coumax,
                                    double **cflmmx,
                                    double **foumax,
                                    double **varrdt,
                                    double **dtmin,
                                    double **dtmax,
                                    double **relxst)
{
  *iptlro = &(_time_step_options.iptlro);
  *idtvar = &(_time_step_options.idtvar);
  *coumax = &(_time_step_options.coumax);
  *cflmmx = &(_time_step_options.cflmmx);
  *foumax = &(_time_step_options.foumax);
  *varrdt = &(_time_step_options.varrdt);
  *dtmin  = &(_time_step_options.dtmin );
  *dtmax  = &(_time_step_options.dtmax );
  *relxst = &(_time_step_options.relxst);
}

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update Kahan compensation.
 *
 * Called as a function to help avoid compiler optimizating Kahan
 * summation out (though aggressive or LTO optimizations might).
 *
 * new copensation term is (a - b) - c
 *
 * \param[in]  a  a
 * \param[in]  b  b
 * \param[in]  c  c
 */
/*----------------------------------------------------------------------------*/

static void
_update_kahan_compensation(double  *a,
                           double  *b,
                           double  *c)
{
  _c = (*a - *b) - *c;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide read/write access to cs_glob_time_step
 *
 * \return  pointer to global time step structure
 */
/*----------------------------------------------------------------------------*/

cs_time_step_t *
cs_get_glob_time_step(void)
{
  return &_time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide read/write access to cs_glob_time_step_options
 *
 * \return  pointer to global time step options structure
 */
/*----------------------------------------------------------------------------*/

cs_time_step_options_t *
cs_get_glob_time_step_options(void)
{
  return &_time_step_options;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define whether time step is variable or not
 *
 * \param[in]  is_variable  0 if time step is variable in time, 1 if it is fixed
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_define_variable(int  is_variable)
{
  if (is_variable > 0)
    _time_step.is_variable = 1;
  else
    _time_step.is_variable = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define whether time step is local in space or not
 *
 * \param[in]  is_local  0 if time step is uniform in space, 1 if it is local
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_define_local(int  is_local)
{
  if (is_local > 0)
    _time_step.is_local = 1;
  else
    _time_step.is_local = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define maximum time step number
 *
 * \param[in]  nt_max  maximum time step number (unlimited if negative)
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_define_nt_max(int  nt_max)
{
  _time_step.nt_max = nt_max;
  _time_step.t_max = -1.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define maximum time value
 *
 * \param[in]  t_max  maximum time value (unlimited if negative)
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_define_t_max(double  t_max)
{
  _time_step.t_max = t_max;
  _time_step.nt_max = -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set time values from previous (usually restarted) calculations
 *
 * \param[in]  nt_prev  previous time step number
 * \param[in]  t_prev   previous physical time
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_define_prev(int     nt_prev,
                         double  t_prev)
{
  _time_step.nt_prev = nt_prev;
  _time_step.nt_cur = nt_prev;
  _time_step.t_prev = t_prev;
  _time_step.t_cur = t_prev;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Increment the global time step.
 *
 * \param[in]  dt  time step value to increment
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_increment(double  dt)
{
  _time_step.dt[2] = _time_step.dt[1];
  _time_step.dt[1] = _time_step.dt[0];
  _time_step.dt[0] = dt;

  double z = dt - _c;
  double t = _time_step.t_cur + z;

  _update_kahan_compensation(&t, &_time_step.t_cur, &z);

  _time_step.t_cur = t;
  _time_step.nt_cur += 1;

  if (_time_step.is_local)
    cs_base_update_status("time step: %d\n",
                          _time_step.nt_cur);
  else
    cs_base_update_status("time step: %d; t = %g\n",
                          _time_step.nt_cur, _time_step.t_cur);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Redefine the current time values.
 *
 * \remark Using \ref cs_time_step_increment is preferred, but this function
 *         may be required for reverting to a previous time step.
 *
 * \param[in]  nt_cur  current time step number
 * \param[in]  t_cur   current physical time
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_redefine_cur(int     nt_cur,
                          double  t_cur)
{
  _time_step.nt_cur = nt_cur;
  _time_step.t_cur = t_cur;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the time stepping options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_step_log_setup(void)
{
  cs_log_printf
    (CS_LOG_SETUP,
     _("\nTime stepping options\n"
       "---------------------\n\n"));

  cs_time_step_type_t ts_type = cs_glob_time_step_options->idtvar;
  int ts_id = ts_type + 1; /* (this enum starts at -1) */

  if (ts_type == CS_TIME_STEP_STEADY) {

    /* Relaxation parameters */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Steady (SIMPLE) algorithm\n\n"
         "    Global parameters\n\n"
         "      idtvar: %21s (%s)\n"
         "      relxst:     %17.5g (Reference relaxation coefficient)\n\n"),
       cs_time_step_type_enum_name[ts_id],
       cs_time_step_type_name[ts_id],
       cs_glob_time_step_options->relxst);
  }

  else {

    if (ts_type ==CS_TIME_STEP_CONSTANT) {
      cs_log_printf
        (CS_LOG_SETUP,
         _("  Unsteady algorithm\n\n"
           "    Time step parameters\n\n"
           "      idtvar: %21s (%s)\n"
           "      dtref:      %17.5g (Reference time step)\n\n"),
         cs_time_step_type_enum_name[ts_id],
         cs_time_step_type_name[ts_id],
         cs_glob_time_step->dt_ref);
    }

    else {
      if (cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) {
        cs_log_printf
          (CS_LOG_SETUP,
           _("  Unsteady algorithm\n\n"));
      } else if (cs_glob_time_step_options->idtvar == CS_TIME_STEP_LOCAL) {
        cs_log_printf
          (CS_LOG_SETUP,
           _("  Space & time varying time step algorithm (pseudo-steady)\n\n"));
      }

      cs_log_printf
        (CS_LOG_SETUP,
         _("  Time step parameters:\n\n"
           "    idtvar: %21s (%s)\n"
           "    iptlro:     %17d (1: rho-related DT clipping)\n"
           "    coumax:     %17.5g (Maximum target CFL)\n"
           "    foumax:     %17.5g (Maximum target Fourier)\n"
           "    varrdt:     %17.5g (For var. DT, max. increase)\n"
           "    dtmin:      %17.5g (Minimum time step)\n"
           "    dtmax:      %17.5g (Maximum time step)\n"
           "    dtref:      %17.5g (Reference time step)\n\n"
           "  When the value of coumax or foumax is negative\n"
           "  or zero, the associated time step limitation\n"
           "  (for CFL and Fourier respectively) is ignored.\n\n"),
         cs_time_step_type_enum_name[ts_id],
         cs_time_step_type_name[ts_id],
         cs_glob_time_step_options->iptlro,
         cs_glob_time_step_options->coumax,
         cs_glob_time_step_options->foumax,
         cs_glob_time_step_options->varrdt,
         cs_glob_time_step_options->dtmin,
         cs_glob_time_step_options->dtmax,
         cs_glob_time_step->dt_ref);

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

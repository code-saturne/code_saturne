/*============================================================================
 * Base time step data.
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
#include "cs_mesh_location.h"
#include "cs_stokes_model.h"

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
        = 0 or 1), \ref t_cur increases of \ref dt (value of the time step) at each iteration.
        If the time step is non-uniform (\ref cs_time_step_options_t::idtvar "idtvar"=2),
        \ref t_cur increases of \ref cs_time_step_t::dt_ref "dt_ref" at each time step.\n
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
        step. If \ref idtvar=0, since the time step is fixed and cannot be
        clipped, this ratio can be greater than 1. When \ref idtvar > 0, this
        ratio will be less than 1, except if the constraint \ref dtmin has
        prevented the code from reaching a sufficiently low value for \ref dt.
        Useful when density gradients and gravity are present.

  \var  cs_time_step_options_t::idtvar
        Option for a variable time step
        - -1: steady algorithm
        -  0: constant time step
        -  1: time step constant in space but variable in time
        -  2: variable time step in space and in time.
        If the numerical scheme is a second-order in time, only the
        option 0 is allowed.

  \var  cs_time_step_t::dt_ref
        Reference time step.\n

        This is the time step value used in the case of a calculation run with a
        uniform and constant time step, i.e. \ref idtvar =0 (restart calculation
        or not). It is the value used to initialize the time step in the case of
        an initial calculation run with a non-constant time step(\ref idtvar=1 or
        2). It is also the value used to initialise the time step in the case of
        a restart calculation in which the type of time step has been changed
        (for instance, \ref idtvar=1 in the new calculation and \ref idtvar = 0
        or 2 in the previous calculation).\n
        See \ref user_initialization_time_step for examples.

  \var  cs_time_step_options_t::coumax
        Maximum Courant number (when \ref idtvar is different from 0).

  \var  cs_time_step_options_t::cflmmx
        Max. Courant number for the continuity equation in compressible model.

  \var  cs_time_step_options_t::foumax
        Maximum Fourier number (when \ref idtvar is different from 0).

  \var  cs_time_step_options_t::varrdt
        Maximum allowed relative increase in the calculated time step value
        between two successive time steps (to ensure stability, any decrease
        in the time step is immediate and without limit).\n
        Useful when idtvar is different from 0.

  \var  cs_time_step_options_t::dtmin
        Lower limit for the calculated time step when idtvar is different from 0.
        Take \ref dtmin = min (ld/ud, sqrt(lt/(gdelta rho/rho)), ...).

  \var  cs_time_step_options_t::dtmax
        Upper limit for the calculated time step when idtvar is different from 0.
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
  .dt = {0., 0., 0.},
  .dt_ref = 0.1,
  .dt_next = 0.1
};

static cs_time_step_options_t  _time_step_options = {
  .iptlro = 0,
  .idtvar = 0,
  .coumax = 1.,
  .cflmmx = 0.99,
  .foumax = 10.,
  .varrdt = 0.1,
  .dtmin  = -1.e12*10.,
  .dtmax  = -1.e12*10.,
  .relxst = 0.7
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
 * \brief Provide acces to cs_glob_time_step
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_time_step_t *
cs_get_glob_time_step(void)
{
  return &_time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide acces to cs_glob_time_step_options
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

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

  /* Steady */
  if (cs_glob_time_step_options->idtvar < 0) {
    /* Time step parameters */
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Steady algorithm\n\n"
         "   Global parameters\n\n"
         "    idtvar:     %14d (-1: steady algorithm)\n"
         "    relxst:     %14.5e (Reference relaxation coefficient)\n\n"),
         cs_glob_time_step_options->idtvar,
         cs_glob_time_step_options->relxst);

    /* Frozen velocity field */
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Frozen velocity field\n\n"
         "    iccvfg:      %14d (1: Frozen velocity field)\n"),
         cs_glob_stokes_model->iccvfg);
  }

  /* Unsteady */
  else {
    cs_log_printf
      (CS_LOG_SETUP,
       _("  Unsteady algorithm\n\n"
         "   Time step parameters\n\n"
         "    idtvar:      %14d (0 cst; 1,2 var (t, t-space)\n"
         "    iptlro:      %14d (1: rho-related DT clipping)\n"
         "    coumax:      %14.5e (Maximum target CFL)\n"
         "    foumax:      %14.5e (Maximum target Fourier)\n"
         "    varrdt:      %14.5e (For var. DT, max. increase)\n"
         "    dtmin:       %14.5e (Minimum time step)\n"
         "    dtmax:       %14.5e (Maximum time step)\n"
         "    dtref:       %14.5e (Reference time step)\n\n"
         "    With a non-constant time step (idtvar = 1 or 2)\n"
         "    when the value of coumax or foumax is negative\n"
         "    or zero, the associated time step limitation\n"
         "    (for CFL and Fourier respectively) is ignored.\n\n"),
         cs_glob_time_step_options->idtvar,
         cs_glob_time_step_options->iptlro,
         cs_glob_time_step_options->coumax,
         cs_glob_time_step_options->foumax,
         cs_glob_time_step_options->varrdt,
         cs_glob_time_step_options->dtmin,
         cs_glob_time_step_options->dtmax,
         cs_glob_time_step->dt_ref);

    /* Frozen velocity field */
    cs_log_printf
      (CS_LOG_SETUP,
       _("   Frozen velocity field\n\n"
         "    iccvfg:      %14d (1: Frozen velocity field)\n"),
         cs_glob_stokes_model->iccvfg);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

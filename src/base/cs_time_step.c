/*============================================================================
 * Base time step data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

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
        absolute time step number reached by previous computation
  \var  cs_time_step_t::nt_cur
        current absolute time step number
  \var  cs_time_step_t::nt_max
        maximum absolute time step number
  \var  cs_time_step_t::nt_ini
        number of time step for initialization
  \var  cs_time_step_t::t_prev
        physical time reached by previous computation
  \var  cs_time_step_t::t_cur
        current absolute time
  \var  cs_time_step_t::t_max
        maximum absolute time

  \struct cs_time_step_options_t

  \brief time step options descriptor

  Members of this time step options descriptor are publicly accessible, to
  allow for concise syntax.

  \var  cs_time_step_options_t::inpdt0
        Indicator "zero time step"
        - 0: standard calculation
        - 1: to simulate no time step
        - for non-restarted computations: only resolution (Navier-Stokes,
          turbulence, scalars) is skipped
        - for restarted computations: resolution, computation of physical
          properties, and definition of boundary conditions is skipped
          (values are read from checkpoint file).

  \var  cs_time_step_options_t::iptlro
        Clip the time step with respect to the buoyant effects
        - 0: false
        - 1: true.

  \var  cs_time_step_options_t::idtvar
        Option for a variable time step
        - -1: steady algorithm
        -  0: constant time step
        -  1: time step constant in space but variable in time
        -  2: variable time step in space and in time.

  \var  cs_time_step_options_t::dtref
        Reference time step.

  \var  cs_time_step_options_t::coumax
        Maximum Courant number (when idtvar is different from 0).

  \var  cs_time_step_options_t::cflmmx
        Max. Courant number for the continuity equation in compressible model.

  \var  cs_time_step_options_t::foumax
        Maximum Fourier number (when idtvar is different from 0).

  \var  cs_time_step_options_t::varrdt
        Allowed relative variation of dt (when idtvar different from 0).

  \var  cs_time_step_options_t::dtmin
        Minimum value of dt (when idtvar is different from 0).
        Take dtmin = min (ld/ud, sqrt(lt/(gdelta rho/rho)), ...).

  \var  cs_time_step_options_t::dtmax
        Maximum value of dt (when idtvar is different from 0).
        Take dtmax = max (ld/ud, sqrt(lt/(gdelta rho/rho)), ...).

  \var  cs_time_step_options_t::relxst
        Relaxation coefficient for the steady algorithm.

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

static cs_time_step_t  _time_step = {0, 0, 0, 0, -1, 2, 0., 0., -1.};

static cs_time_step_options_t  _time_step_options =
{
  0,          /* inpdt0 */
  0,          /* iptlro */
  0,          /* idtvar */
  -1.e12*10., /* dtref  */
  1.,         /* coumax */
  0.99,       /* cflmmx */
  10.,        /* foumax */
  0.1,        /* varrdt */
  -1.e12*10., /* dtmin  */
  -1.e12*10., /* dtmax  */
  0.7         /* relxst */
};

const cs_time_step_t  *cs_glob_time_step = &_time_step;

const cs_time_step_options_t  *cs_glob_time_step_options = &_time_step_options;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_time_step_get_pointers(int     **nt_prev,
                            int     **nt_cur,
                            int     **nt_max,
                            int     **nt_ini,
                            double  **t_prev,
                            double  **t_cur,
                            double  **t_max);

void
cs_f_time_step_options_get_pointers(int    **inpdt0,
                                    int    **iptlro,
                                    int    **idtvar,
                                    double **dtref,
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
 *   t_prev  --> pointer to cs_glob_time_step->t_prev
 *   t_cur   --> pointer to cs_glob_time_step->t_cur
 *   t_max   --> pointer to cs_glob_time_step->t_ax
 *----------------------------------------------------------------------------*/

void
cs_f_time_step_get_pointers(int      **nt_prev,
                            int      **nt_cur,
                            int      **nt_max,
                            int      **nt_ini,
                            double   **t_prev,
                            double   **t_cur,
                            double   **t_max)
{
  *nt_prev = &(_time_step.nt_prev);
  *nt_cur = &(_time_step.nt_cur);
  *nt_max = &(_time_step.nt_max);
  *nt_ini = &(_time_step.nt_ini);
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
 *   inpdt0 --> pointer to cs_glob_time_step_options->inpdt0
 *   iptlro --> pointer to cs_glob_time_step_options->iptlro
 *   idtvar --> pointer to cs_glob_time_step_options->idtvar
 *   dtref  --> pointer to cs_glob_time_step_options->dtref
 *   coumax --> pointer to cs_glob_time_step_options->coumax
 *   cflmmx --> pointer to cs_glob_time_step_options->cflmmx
 *   foumax --> pointer to cs_glob_time_step_options->foumax
 *   varrdt --> pointer to cs_glob_time_step_options->varrdt
 *   dtmin  --> pointer to cs_glob_time_step_options->dtmin
 *   dtmax  --> pointer to cs_glob_time_step_options->dtmax
 *   relxst --> pointer to cs_glob_time_step_options->relxst
 *----------------------------------------------------------------------------*/

void
cs_f_time_step_options_get_pointers(int    **inpdt0,
                                    int    **iptlro,
                                    int    **idtvar,
                                    double **dtref,
                                    double **coumax,
                                    double **cflmmx,
                                    double **foumax,
                                    double **varrdt,
                                    double **dtmin,
                                    double **dtmax,
                                    double **relxst)
{
  *inpdt0 = &(_time_step_options.inpdt0);
  *iptlro = &(_time_step_options.iptlro);
  *idtvar = &(_time_step_options.idtvar);
  *dtref  = &(_time_step_options.dtref );
  *coumax = &(_time_step_options.coumax);
  *cflmmx = &(_time_step_options.cflmmx);
  *foumax = &(_time_step_options.foumax);
  *varrdt = &(_time_step_options.varrdt);
  *dtmin  = &(_time_step_options.dtmin );
  *dtmax  = &(_time_step_options.dtmax );
  *relxst = &(_time_step_options.relxst);
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
  _time_step.nt_cur += 1;
  _time_step.t_cur += dt;
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

END_C_DECLS

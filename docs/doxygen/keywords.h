/*============================================================================
 * Doxygen documentation for specific keywords
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_equation_param.h"
#include "cs_parameters.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file keywords.h
        General keywords
*/

/*----------------------------------------------------------------------------*/

/*!
 * \defgroup field_keywords Field keywords
 */
/*!@{*/

/*!
  \var label

  Optional label associated to the field; if NULL, name will be used instead.
*/
char *label:

/*!
  \var log

  Logging flag for this field (0: off, 1: on)
*/
int log:

/*!
  \var post_vis

  Postprocessing and visualization flag for this field;
  The value may be a combination (sum) of:
  -  0 (not active)
  - \ref CS_POST_ON_LOCATION (output on main location)
  - \ref CS_POST_BOUNDARY_NR (non-reconstructed values postprocessed on boundary
    if main location is cells)
  - \ref CS_POST_MONITOR (output on main probe set)
*/
int post_vis;

/*!
  \var coupled

  If > 0, this variable is coupled using the internal coupling mechanism.

  Restricted to fields with \ref CS_FIELD_VARIABLE type.
*/
int coupled;

/*!
  \var moment_id

  If > -1, refers to the field if of which the current field is a time moment
  (see \ref cs_time_moment.h). If < 0, the current field is not a time moment.

  Restricted to fields with both \ref CS_FIELD_VARIABLE and
  \ref CS_FIELD_POSTPROCESS type.
*/
int moment_id;

/*!
  \var time_extrapolated
  Is the field time-extrapolated?
  - -1: default automatic value
  -  0: "standard" first-order: the value calculated at
        the beginning of the current time step (from the
        variables known at the end of the previous time step) is used
  -  1: second-order: the physical property \f$\phi\f$ is
        extrapolated according to the formula
        \f$\phi^{n+\theta}=[(1+\theta)\phi^n-\theta \phi^{n-1}]\f$, \f$\theta\f$ being
        given by the value of 0.5
  -  2: first-order: the physical property \f$\phi\f$ is
        extrapolated at $n+1$ according to the same formula
        as when = 1 but with \f$\theta\f$ = 1
 */
int time_extrapolated;

/*!
  \var limiter_choice
  Integer corresponding to the type of Roe-Sweby Limiter:
  - 1: minmod
  - 2: Van-Leer
  - 3: Van-Albada
  - 4: superbee

  Restricted to fields with \ref CS_FIELD_VARIABLE type.
*/
int limiter_choice;

/*!
  \var var_cal_opt
  Structure containing the equation parameters of a solved variable.

  Restricted to fields with \ref CS_FIELD_VARIABLE type.
*/
cs_equation_param_t *var_cal_opt;

/*!
  \var solving_info
  Structure containing the solving info of the field variables
  (used for log, not setup, so set NULL setup logging function)

  Restricted to fields with \ref CS_FIELD_VARIABLE type.
 */
cs_solving_info_t *solving_info;

/*!
  \var restart_file
  Indicates in which restart file the associated info may be found.

  If NULL, default rules apply.
*/
char *restart_file;

/*!@}*/

/*----------------------------------------------------------------------------*/


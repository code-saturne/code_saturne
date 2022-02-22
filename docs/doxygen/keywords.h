/*============================================================================
 * Doxygen documentation for specific keywords
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
  \var scalar_id
  Matching scalar id (inverse Fortran <tt> isca(iscal) </tt>),
  or -1 if the field does not represent a solved scalar type variable.
 */
int scalar_id;

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

/*!
  \var diffusivity_id
  Field if of the matching molecular diffusivity for a scalar.
  Negative value if the field has constant diffusivity.

  If set to 0, a matching field will be created and its value
  reset automatically to that field's id. If set directly to
  a value > 0, it is assumed that the matching diffusivity field
  has already been defined and is associated with this scalar.
  This allows both creating an associated field automatically or
  in more advanced cases, sharing a diffusivity field between
  several scalars.
*/
int diffusivity_id;

/*!
  \var diffusivity_ref
  Reference molecular diffusivity for a scalar (\f$kg.m^{-1}.s^{-1}\f$).
  Negative value if not initialized or used.

  \warning: for a temperature, the diffusivity is defined as
  \f$\lambda/C_p\f$ where \f$\lambda\f$ and \f$C_p\f$ are the conductivity
  and specific heat. When using the Graphical Interface, \f$\lambda\f$ and
  \f$C_p\f$ are specified separately, and the matching molecular diffusivity
  is computed automatically.\n
  With the compressible module, \ref diffusivity_ref (given in \ref uscfx2)
  is directly the thermal conductivity \f$W.m^{-1}.K^{-1}\f$.\n
  With gas or coal combustion, the molecular diffusivity of the enthalpy
  (\f$kg.m^{-1}.s^{-1}\f$) must be specified by the user in the variable
  ppthch::diftl0 "diftl0"(\ref cs_user_combustion).\n
  With the electric module, for the Joule effect, the diffusivity is
  specified by the user in \ref cs_user_physical_properties.c (even if
  it is constant). For the electric arcs, it is calculated from the
  thermochemical data file.
*/
double diffusivity_ref;

/*!@}*/

/*----------------------------------------------------------------------------*/


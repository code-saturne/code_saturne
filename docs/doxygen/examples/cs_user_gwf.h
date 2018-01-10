/*============================================================================
 * Code_Saturne documentation page
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

/*!
  \page cs_user_gwf_cdo Settings soils and tracers for the groundwater flow module with CDO schemes (cs_user_gwf.c)

  The add of a new soil or a new tracer takes place in \ref
  cs_user_model (see \ref cs_user_parameters_h_cdo_gwf for more
  details).

  \section cs_user_gwf_h_soil Soils

  Example for a saturated soil defined by an isotropic permeability.

  \snippet cs_user_gwf-example.c param_cdo_gwf_set_soil

  \ref cs_gwf_set_aniso_saturated_soil sets a saturated soil defined
  by an anisotropic permeability. Soils which behave according to a
  Van Genuchten model can be specified using \ref
  cs_gwf_set_aniso_genuchten_soil or \ref
  cs_gwf_set_iso_genuchten_soil. More advanced definition using a
  user-defined model is also possible using \ref cs_gwf_set_user_soil.

  \section cs_user_gwf_h_tracer Tracers

  Here is an example for a standard tracer.

  \snippet cs_user_gwf-example.c param_cdo_gwf_set_tracer

*/

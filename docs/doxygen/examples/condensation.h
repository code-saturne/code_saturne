/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*-----------------------------------------------------------------------------*/

/*!
  \page condens Settings of condensation mass source terms

  \section condens_h_intro Introduction

  Source terms modelling condensation inside the fluid domain
  on internal metal structures and at the
  boundaries can be set respectively through the subroutines
  \ref cs_user_metal_structures_source_terms and
  \ref cs_user_boundary_mass_source_terms.

  \section condens_h_metal_structures Source terms for condensation on internal metal structures

  This model can be enabled in the subroutine \ref usppmo in the file
  \ref cs_user_parameters.f90 as follows:
  \code{.f90}
  if ( ippmod(igmix).ge.0 ) then
    ! Specific condensation modelling
    !      if = -1 module not activated
    !      if =  0 condensation source terms with metal
    !                               structures activate
    icondv = -1
  endif
  \endcode

  The setting of the condensation source terms is then done in the subroutine
  \ref cs_user_metal_structures_source_terms as follows below.

  The following variables need to be declared:

  \snippet cs_user_metal_structures_source_terms-condensation.f90 loc_var_dec

  Necessary species physical properties can be retrieved as follows:

  \snippet cs_user_metal_structures_source_terms-condensation.f90 init

  The zones on which the condensation mass source term will be imposed can be
  defined as follows:

  \snippet cs_user_metal_structures_source_terms-condensation.f90 cells_selection

  Modelling of the metal side (wall thermal model and metal properties can then
  be specified as follows:

  \snippet cs_user_metal_structures_source_terms-condensation.f90 model_settings

  Finally the source term type and values have to be set as follows:

  \snippet cs_user_metal_structures_source_terms-condensation.f90 source_terms_values

  \section condens_h_boundary Boundary source terms for condensation

  The following variables need to be declared:

  \snippet cs_user_boundary_mass_source_terms-nzones_condensation.f90 loc_var_dec

  Necessary species physical properties can be retrieved as follows:

  \snippet cs_user_boundary_mass_source_terms-nzones_condensation.f90 init

  The subroutine \ref cs_user_boundary_mass_source_terms is called three times.

  At the first call the number of boundary faces and number of zones on which a
  boundary mass source term is imposed is computed according to the selection
  criteria prescribed by the user.

  \snippet cs_user_boundary_mass_source_terms-nzones_condensation.f90 zones_definition

  The above part of the subroutine is also executed at the second call. In addition, at
  the second call, condensation models are chosen.

  \snippet cs_user_boundary_mass_source_terms-nzones_condensation.f90 model_settings

  Finally at the third call, the source term types and values have to be set.

  \snippet cs_user_boundary_mass_source_terms-nzones_condensation.f90 source_terms_values

  \section boundary_mass_source-condens Boundary source terms for condensation

  The following variables need to be declared:

  \snippet cs_user_boundary_mass_source_terms-condensation.f90 loc_var_dec

  Necessary species physical properties can be retrieved as follows:

  \snippet cs_user_boundary_mass_source_terms-condensation.f90 init

  The subroutine \ref cs_user_boundary_mass_source_terms is called three times.

  At the first call the number of boundary faces and number of zones on which a
  boundary mass source term is imposed is computed according to the selection
  criteria prescribed by the user.

  \snippet cs_user_boundary_mass_source_terms-condensation.f90 zones_definition

  The above part of the subroutine is also executed at the second call. In addition, at
  the second call, condensation models are chosen.

  \snippet cs_user_boundary_mass_source_terms-condensation.f90 model_settings

  Finally at the third call, the source term types and values have to be set.

  \snippet cs_user_boundary_mass_source_terms-condensation.f90 source_terms_values

*/

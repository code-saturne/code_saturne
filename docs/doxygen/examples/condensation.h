/*============================================================================
 * code_saturne documentation page
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

/*-----------------------------------------------------------------------------*/

/*!
  \page condens Settings of condensation mass source terms

  \section condens_h_intro Introduction

  Source terms modelling condensation inside the fluid domain
  on internal metal structures and at the
  boundaries can be set respectively through the subroutines
  \ref cs_user_metal_structures_source_terms and
  \ref cs_user_wall_condensation.

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

  The condensation of steam on surfaces can be activated by adding the following lines
  in function \ref cs_user_model of file \ref cs_user_parameters.c :

  \snippet cs_user_parameters-base.c wall_condensation

  The subroutine \ref cs_user_wall_condensation is called three times.

  The first call computes the number of boundary faces and the number of zones on which
  a boundary mass source term is imposed, based on the selection criteria prescribed
  by the user.

  In this example, all faces with tag "60" in the mesh are gathered in a single
  condensation zone.
  \snippet cs_user_wall_condensation.c zones_definition

  At the second call, connectivity tables are built between the global mesh numbering
  and the one dedicated to wall condensation (see snippet above).
  In addition, parameters related to the condensation model are set.

  \snippet cs_user_wall_condensation.c model_settings

  At the third call, properties related to the solid wall are set.

  \snippet cs_user_wall_condensation.c solid_wall

  Finally, the source terms associated with the condensation phenomenon are defined.

  \snippet cs_user_wall_condensation.c source_term_values

*/

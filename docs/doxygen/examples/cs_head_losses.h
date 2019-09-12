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

  \page cs_head_losses Examples of data settings for head losses
  (cs_user_zones.c and cs_user_head_losses.c)

  \brief the \ref cs_user_head_losses function is used to compute
  the values of the head loss term, and is called at each time step
  for each previously defined head loss volume zone.

  Volume zones may be defined using the GUI, or through
  the \ref cs_user_zones (in cs_user_zones.c).

  cku is the local head loss term.

 It appears on the momentum as follows:
     \f[ \rho \der{\vect{u}}{t} = - \grad p + \vect{headloss} \: (+\: \text{other terms})\f]
                      with  \f[ \vect{headloss} = - \rho \tens{cku}\cdot \vect{u} \,\:  (\text{in } kg\cdot m^{-2} \cdot s^{-1})\f]

 For a distributed head loss, let \f${ \tens{\xi_l} = \dfrac{\tens{dh_l}}{(0.5 \rho  u^2)}}\f$ given by the litterature
    (\f$ \tens{dh_l} \f$ is the head loss per unit length)

    the source term \c tspdc is equal to \f$\tens{dh_l} = - \tens{\xi_l}(0.5\rho\vect{u}^2)\f$

    we have \f$ \tens{cku} = 0.5\tens{\xi_l}|\vect{U}| \f$


 For a singular head loss, let \f$\tens{\xi_l} = \dfrac{\tens{dh_s}}{0.5\rho\vect{u}^2}\f$ given by the litterature
    (\f$\tens{dh_s} \f$ is the singular head loss)

    the source term \c tspdc is equal to \f[\frac{\tens{dh_s}}{L} = - \frac{\tens{\xi_l}}{L} (0.5 \rho\vect{u}^2)\f]. We have \f[\tens{cku} = 0.5\frac{\tens{\xi_s}}{L}|\vect{u}|\f]

    where \f$ L \f$ is the length over which we have chosen to represent the
    singular head loss.


  \section cs_user_head_losses_examples Head loss setting examples

  Here is the list of examples:

  - \subpage base_head_losses_examples

*/
// _____________________________________________________________________________
/*!


  \page base_head_losses_examples Basic examples

  \section init_and_final Initialization and finalization

  It is useful to map a field array to a local pointer for a clear and concise
  access, such as done here for the velocity:

  \snippet cs_user_head_losses.c map_field_arrays

  Otherwise, the zone entries (see \ref cs_volume_zone_t) should contain
  the necessary information with no additional preparation.

  \section body Body

  \subsection beginning Defining a volume zone

  A volume zone may be defined using the GUI, or in the \ref cs_user_zones
  user function (in cs_user_zones.c), such as the following zone determined
  by a geometric criterion:

  \snippet cs_user_zones.c user_zones_head_loss_1

  Note that if the \ref CS_VOLUME_ZONE_HEAD_LOSS flag is not set
  (or the matching type set through the GUI), calls to \ref cs_user_head_losses
  will ignore this zone.

  \subsection head_loss_examples Head loss examples

  Note that in the following examples, we check the zone name, so we
  know which zone we are dealing with using in case of multiple zones.

  head loss tensor coefficients for each cell are organized as follows:
  \c cku11, \c cku22, \c cku33, \c cku12, \c cku13, \c cku23.

  Coefficients are set to zero (then computed based on definitions provided
  through the GUI if this is the case) before calling this function, so
  setting values to zero is usually not necessary, unless we want to fully
  overwrite a GUI-based definition.

  Note that diagonal coefficients must be positive; the calculation may
  crash if this is not the case.

  \subsection diagonal_tensor Example 1: head losses alined with an axis of the computation frame

  Using the previously defined zone, we define head losses in direction \c x

  \snippet cs_user_head_losses.c head_loss_1

  \subsection alpha_tensor Example 2: head losses along a direction at 45 degrees

  Necessary, we shall use here a 3x3 tensor to impose head losses at an angle \f$ alpha = 45^{o} \f$ with respect to x and y
  direction of the computation frame. Namely, resistance is set along components \c x by \c cku1 and \c y by \c cku2 \n.

  \image html orthogonal_reference_frame_sketch.gif "Orthogonal reference frame sketch"

  In the present example, it is chosen to set a head loss representing friction along \c X and to model a vane in Y direction
  by setting <tt> ck1 = 0 </tt>.

  \snippet cs_user_head_losses.c head_loss_2

*/
// _____________________________________________________________________________

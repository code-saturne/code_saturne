/*============================================================================
 * Code_Saturne documentation page
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

/*-----------------------------------------------------------------------------*/



/*!

  \page cs_head_losses Examples of data settings for head losses (cs_user_head_losses.f90)

  \brief User \ref cs_user_head_losses subroutine which is called at three
         different stages in the code <em>(iappel=1, 2 or 3)</em>.

   -  \c iappel = 1:
        Calculation of the number of cells where a head loss term is
        imposed: \c ncepdp.
        Called once at the beginning of the calculation.

   -  \c iappel = 2:
        Identification of the cells where a head loss term is imposed:
        array \c icepdc(ncepdc).
        Called once at the beginning of the calculation.

   -  \c iappel = 3:
        Calculation of the values of the head loss term.
        Called at each time step.

 Note that calling this subroutine completely overwrites head losses
 defined using the GUI.

  ckupdc is the local head loss term.

 It appears on the momentum as follows:
     \f[ \rho \der{\vect{u}}{t} = - \grad p + \vect{headloss} \: (+\: \text{other terms})\f]
                      with  \f[ \vect{headloss} = - \rho \tens{ckupdc}\cdot \vect{u} \,\:  (\text{in } kg\cdot m^{-2} \cdot s^{-1})\f]

 For a distributed head loss, let \f${ \tens{\xi_l} = \dfrac{\tens{dh_l}}{(0.5 \rho  u^2)}}\f$ given by the litterature
    (\f$ \tens{dh_l} \f$ is the head loss per unit length)

    the source term \c tspdc is equal to \f$\tens{dh_l} = - \tens{\xi_l}(0.5\rho\vect{u}^2)\f$

    we have \f$ \tens{ckupdc} = 0.5\tens{\xi_l}|\vect{U}| \f$


 For a singular head loss, let \f$\tens{\xi_l} = \dfrac{\tens{dh_s}}{0.5\rho\vect{u}^2}\f$ given by the litterature
    (\f$\tens{dh_s} \f$ is the singular head loss)

    the source term \c tspdc is equal to \f[\frac{\tens{dh_s}}{L} = - \frac{\tens{\xi_l}}{L} (0.5 \rho\vect{u}^2)\f]. We have \f[\tens{ckupdc} = 0.5\frac{\tens{\xi_s}}{L}|\vect{u}|\f]

    where \f$ L \f$ is the length over which we have chosen to represent the
    singular head loss.


  \section cs_user_head_losses_examples Head loss setting examples

  Here is the list of examples:

  - \subpage base_head_losses_examples

*/
// _____________________________________________________________________________
/*!


  \page base_head_losses_examples Basic examples

  \section base_loc_var Local variables to be added

  The following local variables need to be defined for the examples
  in this section:

  \snippet cs_user_head_losses.f90 loc_var_dec

  \section init_and_final Initialization and finalization

  The following initialization block needs to be added for the following examples:

  \snippet cs_user_head_losses.f90 allocate

  At the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_head_losses.f90 deallocate

  In theory Fortran 95 deallocates locally-allocated arrays automatically, but
	deallocating arrays in a symetric manner to their allocation is good pratice,
	and avoids using a different logic for C and Fortran.


  Map field array

  \snippet cs_user_head_losses.f90 map_field_array


  \section body Body

  \subsection beginning Calculation and identification of the number of cells with imposed head loss term

  \snippet cs_user_head_losses.f90 start_1

  2 calls:
  - \c iappel = 1: Calculation of the number of cells where a head loss term is imposed: \c ncepdp. Called once at the beginning of the calculation.

  - \c iappel = 2: Identification of the cells where a head loss term is imposed: array \c icepdc(ncepdc). Called once at the beginning of the calculation.

\note
  - Do not use ckupdc in this section (it is defined with <tt> iappel = 3) </tt>
  - Use icepdc in this section only with <tt> (iappel = 2) </tt>

To be completed by the user: cell selection


  \subsection head_loss_examples Head loss examples
  \subsubsection example_1 Example 1: No head loss (default)

  \snippet cs_user_head_losses.f90 example_1


  \subsubsection example_2 Example 2: Head losses defined by coordinates for zone

  (4 <= x; 2 <= y <= 8) \n No head losses else.

  \snippet cs_user_head_losses.f90 example_2


  \subsection generic Generic subsection

  For <tt> iappel = 1 </tt>, define \c ncepdp, the number of cells with head
	losses. This is valid for both examples above.

  \snippet cs_user_head_losses.f90 start_2

  Defining the number of cells with head losses

  \snippet cs_user_head_losses.f90 generic_subsection_1

  Computing the head loss coefficient values

  Third call, at each time step

  - iappel = 3:

      \c ckupdc: compute head loss coefficients in the calculation coordinates,
              organized in order <tt> k11, k22, k33, k12, k13, k23 </tt>

  Note:

   - make sure diagonal coefficients are positive. The calculation
        may crash if this is not the case, and no further check will
        be done


   \subsection diagonal_tensor Example 1: head losses in direction x

 Diagonal tensor : Example of head losses in direction \c x

  \snippet cs_user_head_losses.f90 example_3


  \subsection alpha_tensor Example 2: alpha = 45 degres

  3x3 tensor: Example of head losses at alpha = 45 degres x,y
 direction \c x resists by \c ck1 and \c y by \c ck2 \n
 <tt> ck2 = 0 </tt> represents vanes as follows:
in coordinate system \c x, \c y

  \image html orthogonal_reference_frame_sketch.gif "Orthogonal reference frame sketch"

  \snippet cs_user_head_losses.f90 example_4


  \subsection end Filling the cells of the head loss matrix

  \snippet cs_user_head_losses.f90 filling

*/
// _____________________________________________________________________________

/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

  \page cs_user_radiative_transfer Examples of data settings for radiative transfers

\image html radiative_tr_sketch.gif "Sketch of thermal flux in boundary walls"


The radiative boundary condition is based on the calculation of a new wall
temperature. This temperature is  computed with a thermal flux balance:

\f[{ Q_{conduction} = Q_{convection} + (Q_{rayt_{absorption}} - Q_{rayt_{emission}}}) \f]

Therefore :

\f[ \dfrac{xlamp}{epap}  (T_{fluid} - T_{wall})
= h_{fluid}  (T_{fluid} - T_{wall}) + epsp  (Q_{incid} - \sigma * T_{wall}) \f]


 \note In \c Code_Saturne the flux is positive when it is oriented from inside to outside.


  |  Corps                       |     Emissivity    |
  |------------------------------|------------------:|
  |  polished steel              |       0.06        |
  |  oxidized steel              |       0.80        |
  |  steel rough                 |       0.94        |
  |  polished aluminium          |       0.04        |
  |  oxidiezd aluminium (inside) |       0.09        |
  |  oxidized aluminium (wet air)|       0.90        |
  |  brick                       |       0.93        |
  |  concrete                    |       0.93        |
  |  paper                       |       0.8 to 0.9  |
  |  water                       |       0.96        |


  \section bound_faces Boundary faces identification

   Boundary faces may be identified using the \ref getfbr subroutine. The syntax of this subroutine is described in the \ref cs_user_boundary_conditions subroutine, but a more thorough description can be found in the user guide.

\note These usefull constant are definded \n
\f$ TKELVI = 273.16D0 \f$ \n
\f$ SIG = 5.6703D-8 \f$


  \section init_fin Initialization and finalization


   The following initialization block needs to be added for the following examples:

  \snippet cs_user_radiative_transfer_bcs.f90 allocate

  At the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_radiative_transfer_bcs.f90 deallocate

  In theory Fortran 95 deallocates locally-allocated arrays automatically, but deallocating arrays in a symetric manner to their allocation is good pratice, and avoids using a different logic for C and Fortran.

 <b> Remaining initialisation</b>


  ivar: number of the thermal variable

   \snippet cs_user_radiative_transfer_bcs.f90 ivar

  Min and Max values for the wall temperatures (clipping otherwise)

 \f$ T_{min} \f$ and \f$T_{max} \f$ are given in Kelvin.

  \snippet  cs_user_radiative_transfer_bcs.f90 temp

   \section assign2 Assign boundary conditions to boundary wall

    \subsection zone_def Zones definition

 We define zones of wall boundary, and we assign a type.
   This allows to apply the boundary conditions and realize
   balance sheets by treating them separately for each zone.

 For each boundary face ifac (not just the faces of wall)
   the user defines his own choice by a number of zone
   \c izfrdp(ifac) from color of the boundary face
     or more generally, their properties (color, groups ...),
     or boundary conditions specified in \ref cs_user_boundary_conditions,
     or even of their coordinates.

\warning It is essential that ALL boundary faces
   have been assigned to a zone.
   The number of zones (the value of \c izfrdp(ifac)) is
   arbitrarily chosen by the user, but must be a
   positive integer and less than or equal to \c nbzrdm
   (value set in parameter \ref radiat.h).


\subsection wall_carac Wall caracteristics

\warning The unit of the temperature is the Kelvin

\subsubsection manda Mandatory data

   - \c isothp(ifac) boundary face type
               -  \c itpimp -> Gray wall with fixed inside temperature
               -  \c ipgrno -> Gray wall with fixed outside temperature
               -  \c iprefl -> Reflecting wall with fixed outside temperature
               -  \c ifgrno -> Gray wall with fixed conduction flux
               -  \c ifrefl -> Reflecting wall with fixed conduction flux

   - \c tintp(ifac) inside wall temperature (Kelvin)
                  initialize thwall at the first time step.
                  If \c isothp = \c itpimp, the value of thwall is fixed to \c tintp
                  In the other case, \c tintp is only for initialization.
\subsubsection data Other data (depending of the isothp)

  - \c rcodcl = conduction flux
  - \c epsp   = emissivity
  - \c xlamp  = conductivity (\f$W.m^{-1}.K^{-1}\f$)
  - \c epap   = thickness (\f$m\f$)
  - \c textp  = outside temperature (\f$K\f$)


\section ex Examples of boundary conditions

Here is a list of examples:

  \subsection ex1 Gray or black wall with profil of fixed inside temperature

  For wall boundary faces, selection criteria: color 1  \n

  \snippet   cs_user_radiative_transfer_bcs.f90  example_1


  \subsection ex2 Gray or black wall with fixed outside temperature \f$ T_{ext} \f$

  For wall boundary faces, selection criteria: color 2  \n

  \snippet  cs_user_radiative_transfer_bcs.f90 example_2


  \subsection ex3 Reflecting wall (\f$ epsp = 0 \f$) with fixed outside temperature \f$ T_{ext} \f$

  For wall boundary faces, selection criteria: color 3 \n

  \snippet  cs_user_radiative_transfer_bcs.f90 example_3


  \subsection ex4 Gray or black wall and fixed conduction flux through the wall

For wall boundary faces which have the color 4: \n

\f[
\begin{array}{rcl}
\frac{\texttt{xlamp}}{\texttt{epap}} \cdot (T_{wall} - T_{ext})
&=& \text{fixed conduction flux in } W.m^{-2} \\
&=& \texttt{rodcl(ifac,ivar,3)}
\end{array}
\f]

If the conduction flux is zero then the wall is adiabatic. The array \f$ \texttt{rcodcl(ifac,ivar,3)}\f$ has the value of the flux. \n
Flux density (< 0 if gain for the fluid)
 - For temperature \f$T\f$, in \f$ W.m^{-2}\f$:

 \f[ rcodcl(ifac,ivar,3)=C_p (viscls+\frac{visct}{\sigma})\cdot \grad{T}\cdot \vect{n} \f]

 - For enthalpy \f$h\f$, in \f$ W.m^{-2} \f$:
 \f[ RCODC(IFAC,IVAR,3)=(viscls+\frac{visct}{\sigma})\cdot \grad{H} \cdot \vect{n}\f]

\snippet  cs_user_radiative_transfer_bcs.f90 example_4


\subsection ex5 Reflecting wall and fixed conduction flux through the wall

For wall boundary faces which have the color 5:\n

\f[
\frac{xlamp}{epap} \cdot (T_{wall} - T_{ext}) = \text{fixed conduction flux}
\f]
and \f$ epsp = 0 \f$

If the conduction flux is zero then the wall is adiabatic.
 Flux density (< 0 if gain for the fluid)
  - For temperatures \f$T\f$,    in \f$ W.m^{-2} \f$:
    \f[  rcodcl(ifac,ivar,3) = C_p  (viscls+\frac{visct}{\sigma}) \cdot \grad{T}\cdot \vect{n} \f]
  - For enthalpies \f$h\f$,      in \f$ W.m^{-2} \f$:
    \f[  rcodcl(ifac,ivar,3) =    (viscls+\frac{visct}{\sigma})  \cdot \grad{H} \cdot \vect{n} \f]

\snippet  cs_user_radiative_transfer_bcs.f90 example_5


\section w Warning

For all boundary faces that are not wall it is MANDATORY to impose a number of
zone in the array \c izfrdp. For each zone, informations will be displayed in the listing.

\snippet cs_user_radiative_transfer_bcs.f90 w

Verification that all boundary faces have been treated

\snippet  cs_user_radiative_transfer_bcs.f90 check


\section end_loop End of the loop on the boundary faces

\snippet  cs_user_radiative_transfer_bcs.f90 end_radiative
\section format_radiative_trans Format
\snippet cs_user_radiative_transfer_bcs.f90 format_radiative

*/

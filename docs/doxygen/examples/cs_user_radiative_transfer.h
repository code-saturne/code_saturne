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

/*----------------------------------------------------------------------------*/

/*!

  \page cs_user_radiative_transfer Examples of data settings for radiative transfers

  \section radiat_activ Activation of the module

  The module can be activated in the \ref usppmo routine in
  \ref cs_user_parameters.f90. The corresponding keyword is \c iirayo in the
  \ref cs_glob_rad_transfer_options structure.

  This member can take the values:
   - \ref radiat::iirayo "iirayo" = 0: module desactivated.
   - \ref radiat::iirayo "iirayo" = 1: the module is activated and the Discrete
                                       Ordinates Method is used.
   - \ref radiat::iirayo "iirayo" = 2: the module is activated and the P1 model
                                       is used.

  \section radiat_param Radiation module specific parameters.

  When the module is activated, its specific input parameters should be set in
  the \ref cs_user_radiative_transfer_parameters function of the
  \ref cs_user_radiative_transfer.c file.

  \section cs_user_radiative_transfer_h_cs_user_radiative_transfer_parameters Calculation options for the radiative transfer module.

  Radiative transfer parameters may be defined using the
  \ref cs_user_radiative_transfer_parameters function.

  \snippet cs_user_radiative_transfer.c cs_user_radiative_transfer_parameters

  \section cs_user_radiative_transfer_h_boundary_conditions Radiative transfer boundary conditions

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


  \subsection bound_faces Boundary faces identification

   Boundary faces may be identified using the \ref getfbr subroutine. The syntax of this subroutine is described in the \ref cs_user_boundary_conditions subroutine, but a more thorough description can be found in the user guide.

\note These usefull constant are definded \n
\f$ TKELVI = 273.16D0 \f$ \n
\f$ SIG = 5.6703D-8 \f$


  \subsection init_fin Initialization and finalization


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

   \subsection assign2 Assign boundary conditions to boundary wall

    \subsubsection zone_def Zones definition

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


\paragraph wall_carac Wall caracteristics

\warning The unit of the temperature is the Kelvin

\paragraph manda Mandatory data

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
\paragraph data Other data (depending of the isothp)

  - \c rcodcl = conduction flux
  - \c epsp   = emissivity
  - \c xlamp  = conductivity (\f$W.m^{-1}.K^{-1}\f$)
  - \c epap   = thickness (\f$m\f$)
  - \c textp  = outside temperature (\f$K\f$)


\subsection ex Examples of boundary conditions

Here is a list of examples:

  \subsubsection ex1 Gray or black wall with profil of fixed inside temperature

  For wall boundary faces, selection criteria: color 1  \n

  \snippet   cs_user_radiative_transfer_bcs.f90  example_1


  \subsubsection ex2 Gray or black wall with fixed outside temperature \f$ T_{ext} \f$

  For wall boundary faces, selection criteria: color 2  \n

  \snippet  cs_user_radiative_transfer_bcs.f90 example_2


  \subsubsection ex3 Reflecting wall (\f$ epsp = 0 \f$) with fixed outside temperature \f$ T_{ext} \f$

  For wall boundary faces, selection criteria: color 3 \n

  \snippet  cs_user_radiative_transfer_bcs.f90 example_3


  \subsubsection ex4 Gray or black wall and fixed conduction flux through the wall

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


\subsubsection ex5 Reflecting wall and fixed conduction flux through the wall

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

\subsubsection w Warning

For all boundary faces that are not wall it is MANDATORY to impose a number of
zone in the array \c izfrdp. For each zone, informations will be displayed in the listing.

\snippet cs_user_radiative_transfer_bcs.f90 w

Verification that all boundary faces have been treated

\snippet  cs_user_radiative_transfer_bcs.f90 check

\subsection end_loop End of the loop on the boundary faces

\snippet  cs_user_radiative_transfer_bcs.f90 end_radiative
\subsection format_radiative_trans Format
\snippet cs_user_radiative_transfer_bcs.f90 format_radiative

\section abso_flux Absorption coefficient and net radiation flux

The absorption coefficient and the net radiation flux for the radiative module
can be defined in \ref cs_user_radiative_transfer.f90 through the \ref usray3
and \ref usray5 subroutines.

\subsection abso Absorption coefficient

The absorption coefficient is defined in \ref usray3 .

\subsubsection arg Arguments of usray3

\snippet cs_user_radiative_transfer.f90 arg_1

\subsubsection var Local variables to be added

\snippet cs_user_radiative_transfer.f90 loc_var_dec_1

\subsubsection abso_coeff_computation Computation of the absorption coefficient

\snippet cs_user_radiative_transfer.f90 abso_coeff

\subsubsection format_1 Format

\snippet cs_user_radiative_transfer.f90 format_1

\subsection usray5_flux Net radiation flux

The net radiation flux is computed in \ref usray5 .

\subsubsection arg2 Arguments of usray5

\snippet cs_user_radiative_transfer.f90 arg_2

\subsubsection var2 Local variables to be added

\snippet cs_user_radiative_transfer.f90 loc_var_dec_2

\subsubsection init Initialization

At the end of the subroutine, if \c iok is different from zero, some faces
have been forgotten and the calculation stops.

\snippet cs_user_radiative_transfer.f90 init

\subsubsection net_flux_computation Computation of the net radiation flux

\snippet cs_user_radiative_transfer.f90 net_flux

\subsubsection format_2 Format

\snippet cs_user_radiative_transfer.f90 format_2

*/

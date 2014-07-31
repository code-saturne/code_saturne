/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

  \page cs_user_radiative_transfert_parameters Examples of data settings for radiative tranferts (usray2.f90)
  
\image html radiative_tr_sketch.gif "Sketch of thermal flux in boundary walls"


The radiative boundary condition is based on the calculation of a new wall temperature. This temperature is  computed with a thermal flux balance:

\f[{ Q_{conduction} = Q_{convection} + (Q_{rayt_{absorption}} - Q_{rayt_{emission}}}) \f]

Therefore :

\f[ \dfrac{XLAMP}{EPAP}  (T_{fluide} - T_{parop}) = H_{fluide}  (T_{fluide} - T_{parop}) + EPSP  (Q_{incid} - \sigma * T_{parop}) \f]


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


  \section bound_faces Bondary faces identification

   Boundary faces may be identified using the \ref getfbr subroutine. The syntax of this subroutine is described in the \ref cs_user_boundary_conditions subroutine, but a more thorough description can be found in the user guide.

\note These usefull constant are definded \n
\f$ TKELVI = 273.16D0 \f$ \n
\f$ SIG = 5.6703D-8 \f$ 


  \section init_fin Initialization and finalization


   The following initialization block needs to be added for the following examples:

  \snippet cs_user_radiative_transfert_parameters.f90 allocate

  At the end of the subroutine, it is recommended to deallocate the work array:

  \snippet cs_user_radiative_transfert_parameters.f90 deallocate
 
  In theory Fortran 95 deallocates locally-allocated arrays automatically, but deallocating arrays in a symetric manner to their allocation is good pratice, and avoids using a different logic for C and Fortran.

 <b> Remaining initialisation</b>


  IVAR: number of the thermal variable

   \snippet cs_user_radiative_transfert_parameters.f90 ivar

  Min and Max values for the wall temperatures (clipping otherwise)
 
 \f$ T_{min} \f$ and \f$^T_{max} \f$ are given in Kelvin.

  \snippet  cs_user_radiative_transfert_parameters.f90 temp

   \section assign2 Assign boundary conditions to boundary wall

    \subsection zone_def Zones definition

 We define zones of wall boundary, and we assign a type.
   This allows to apply the boundary conditions and realize
   balance sheets by treating them separately for each zone.

 For each boundary face ifac (not just the faces of wall)
   the user defines his own choice by a number of zone
   \c IZFRDP(ifac) from color of the boundary face
     or more generally, their properties (color, groups ...),
     or boundary conditions specified in \ref cs_user_boundary_conditions,
     or even of their coordinates.
\warning It is essential that ALL boundary faces
   have been assigned to a zone.
 The number of zones (the value of \c IZFRDP(ifac)) is
   arbitrarily chosen by the user, but must be a
   positive integer and less than or equal to \c NBZRDM
   (value set in parameter \ref radiat.h).
 


\subsection wall_carac Wall caracteristics

\warning The unit of the temperature is the Kelvin

\subsubsection manda Mandatory data

   - \ref isothp(ifac) boundary face type
               -   = \ref itpimp -> Gray wall with fixed inside temperature
               -   = \ref ipgrno -> Gray wall with fixed outside temperature
               -   = \ref iprefl -> Reflecting wall with fixed outside temperature
               -   = \ref ifgrno -> Gray wall with fixed conduction flux
               -   = \ref ifrefl -> Reflecting wall with fixed conduction flux

   - \ref tintp(ifac) inside wall temperature (Kelvin)
                  initialize thwall at the first time step.
                  If \ref isothp = \ref itpimp, the value of thwall is fixed to \c tintp
                  In the other case, \c tintp is only for initialization.
\subsubsection data Other data (depending of the isothp)

  - \c rcodcl = conduction flux
  - \c epsp   = emissivity
  - \c xlamp  = conductivity (\f$W \cdot m^{-1} \cdot K^{-1}\f$)
  - \c epap   = thickness (\f$m\f$)
  - \c textp  = outside temperature (\f$K\f$)
\section ex Examples

Here is a list of examples:

   - \subpage radiative_transfert_parameters_examples



*/
// __________________________________________________________________________________
/*!
 \page radiative_transfert_parameters_examples Radiative transfert parameters examples

  \section ex1 Example 1

  For wall boundary faces, selection criteria: color 1  \n
  Gray or black wall with profil of fixed inside temperature

  \snippet   cs_user_radiative_transfert_parameters.f90  example_1
  
   \section ex2 Example 2

  For wall boundary faces, selection criteria: color 2  \n
  Gray or black wall with fixed outside temperature \f$ T_{EXTP} \f$
   
  \snippet  cs_user_radiative_transfert_parameters.f90 example_2

  \section ex3 Example 3 
  
  For wall boundary faces, selection criteria: color 3 \n
  Reflecting wall (EPSP = 0) with fixed outside temperature \f$ T_{EXTP} \f$

  \snippet  cs_user_radiative_transfert_parameters.f90 example_3

  \section ex4 Example 4

For wall boundary faces which have the color 4: \n
Gray or black wall and fixed conduction flux through the wall

\f[ \frac{\vect{XLAMP}}{EPAP} \cdot (T_{parop} - T_{extp}) = \text{fixed conduction flux in } W\cdot m^{-2} \f]
\f[ \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: \: = \vect{RODCL}(IFAC,IVAR,3) \f]

If the conduction flux is zero then the wall is adiabatic. The array \f$ \vect{RCODCL}(IFAC,IVAR,3)\f$ has the value of the flux. \n
Flux density (< 0 if gain for the fluid)
 - For temperature T, in \f$ W\cdot m^{-2}\f$:
 
 \f[ \vect{RCODCL}(IFAC,IVAR,3)=C_p (VISCLS+\frac{VISCT}{\sigma})\cdot \grad{T} \f]

 - For enthalpy H, in \f$ W \cdot m^{-2} \f$:
 \f[ \vect{RCODC}(IFAC,IVAR,3)=(VISCLS+\frac{VISCT}{\sigma})\cdot \grad{H} \f]


\snippet  cs_user_radiative_transfert_parameters.f90 example_4

\section ex5 Example 5

For wall boundary faces which have the color 5:\n
reflecting wall and fixed conduction flux through the wall

\f[ \frac{\vect{XLAMP}}{EPAP} \cdot (T_{parop} - T_{extp}) = \text{fixed conduction flux}\f]
and \f$ EPSP = 0 \f$

If the conduction flux is zero then the wall is adiabatic.
 Flux density (< 0 if gain for the fluid)
  - For temperatures T,    in \f$ W\cdot m^{-2} \f$:
    \f[  \vect{RCODCL}(IFAC,IVAR,3) = C_p  (VISCLS+\frac{VISCT}{\sigma}) \cdot \grad{T} \f]
  - For enthalpies H,      in \f$ W \cdot m^{-2} \f$:
    \f[  \vect{RCODCL}(IFAC,IVAR,3) =    (VISCLS+\frac{VISCT}{\sigma})  \cdot \grad{H} \f]

\snippet  cs_user_radiative_transfert_parameters.f90 example_5

\section w Warning 

For all boundary faces that are not wall it is MANDATORY to impose a number of zone in the array \ref izfrdp. For each zone, informations will be displayed in the listing.

\snippet cs_user_radiative_transfert_parameters.f90 w

\section ex6 Example 6 

Verification that all boundary faces have been treated

\snippet  cs_user_radiative_transfert_parameters.f90 example_6 

\section end_loop End of the loop on the boundary faces

\snippet  cs_user_radiative_transfert_parameters.f90 end_radiative 
\section format_radiative_trans Format
\snippet cs_user_radiative_transfert_parameters.f90 format_radiative

*/  

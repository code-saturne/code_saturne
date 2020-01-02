!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Function:
! ---------

!> \file cs_user_boundary_conditions.f90
!>
!> \brief User subroutine which fills boundary conditions arrays
!> (\c icodcl, \c rcodcl) for unknown variables.
!>
!> See \subpage cs_user_boundary_conditions_examples for examples.
!>
!> \section cs_user_boundary_conditions_intro Introduction
!>
!> Here one defines boundary conditions on a per-face basis.
!>
!> Boundary faces may be selected using the \ref getfbr subroutine.
!>
!> \code getfbr(string, nelts, lstelt) \endcode
!>  - string is a user-supplied character string containing selection criteria;
!>  - nelts is set by the subroutine. It is an integer value corresponding to
!>    the number of boundary faces verifying the selection criteria;
!>  - lstelt is set by the subroutine. It is an integer array of size nelts
!>    containing the list of boundary faces verifying the selection criteria.
!>
!> string may contain:
!>  - references to colors (ex.: 1, 8, 26, ...)
!>  - references to groups (ex.: inlet, group1, ...)
!>  - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!>
!> These criteria may be combined using logical operators (\c and,\c or) and
!> parentheses.
!>
!> \par Example
!> \code 1 and (group2 or group3) and y < 1 \endcode
!> will select boundary faces
!> of color 1, belonging to groups 'group2' or 'group3' and with face center
!> coordinate y less than 1.
!>
!> Operators priority, from highest to lowest:
!>  '( )' > 'not' > 'and' > 'or' > 'xor'
!>
!> Similarly, interior faces and cells can be identified using the \ref getfac
!> and \ref getcel subroutines (respectively). Their syntax are identical to
!> \ref getfbr syntax.
!>
!> For a more thorough description of the criteria syntax, see the user guide.
!>
!>
!> \section bc_types Boundary condition types
!>
!> Boundary conditions may be assigned in two ways.
!>
!>
!> \subsection std_bcs For "standard" boundary conditions:
!>
!> One defines a code in the \c itypfb
!> array (of dimensions number of boundary faces).
!> This code will then be used by a non-user subroutine to assign the
!> following conditions.
!> The available codes are:
!>  - \c ientre: Inlet
!>  - \c isolib: Free outlet
!>  - \c isymet: Symmetry
!>  - \c iparoi: Wall (smooth)
!>  - \c iparug: Rough wall
!>
!> These integers are defined elsewhere (in paramx.f90 module).
!> Their value is greater than or equal to 1 and less than or  equal to
!> ntypmx (value fixed in paramx.h)
!>
!> In addition, some values must be defined:
!>  - Inlet (more precisely, inlet/outlet with prescribed flow, as the flow
!>           may be prescribed as an outflow):
!>    - Dirichlet conditions on variables other than pressure are mandatory
!>      if the flow is incoming, optional if the flow is outgoing (the code
!>      assigns zero flux if no Dirichlet is specified); thus,
!>      at face \c ifac, for the variable \c ivar: \c rcodcl(ifac, ivar, 1)
!>
!>
!>  - Smooth wall: (= impermeable solid, with smooth friction)
!>    - Velocity value for sliding wall if applicable:
!>                  - \c rcodcl(ifac, iu, 1) = fluid velocity in the x direction
!>                  - \c rcodcl(ifac, iv, 1) = fluid velocity in the y direction
!>                  - \c rcodcl(ifac, iw, 1) = fluid velocity in the z direction
!>    - Specific code and prescribed temperature value at wall if applicable:
!>                  - \c icodcl(ifac, ivar)    = 5
!>                  - \c rcodcl(ifac, ivar, 1) = prescribed temperature
!>    - Specific code and prescribed flux value at wall if applicable:
!>                  - \c icodcl(ifac, ivar)    = 3
!>                  - \c rcodcl(ifac, ivar, 3) = prescribed flux
!>    .
!>    Note that the default condition for scalars (other than k and epsilon)
!>    is homogeneous Neumann.
!>
!>
!>  - Rough wall: (= impermeable solid, with rough friction)
!>    - Velocity value for sliding wall if applicable:
!>                  - \c rcodcl(ifac, iu, 1) = fluid velocity in the x direction
!>                  - \c rcodcl(ifac, iv, 1) = fluid velocity in the y direction
!>                  - \c rcodcl(ifac, iw, 1) = fluid velocity in the z direction
!>    - Value of the dynamic roughness height to specify in
!>                  - \c rcodcl(ifac, iu, 3)
!>    - Value of the scalar roughness height (if required) to specify in
!>                  - \c rcodcl(ifac, iv, 3) (values for iw are not used)
!>    - Specific code and prescribed temperature value at wall if applicable:
!>                  - \c icodcl(ifac, ivar)    = 6
!>                  - \c rcodcl(ifac, ivar, 1) = prescribed temperature
!>    - Specific code and prescribed flux value at rough wall, if applicable:
!>                  - \c icodcl(ifac, ivar)    = 3
!>                  - \c rcodcl(ifac, ivar, 3) = prescribed flux
!>    .
!>    Note that the default condition for scalars (other than k and epsilon)
!>    is homogeneous Neumann.
!>
!>  - Symmetry (= slip wall):
!>    - Nothing to specify
!>
!>  - Free outlet (more precisely free inlet/outlet with prescribed pressure)
!>    - Nothing to prescribe for pressure and velocity. For scalars and
!>      turbulent values, a Dirichlet value may optionally be specified.
!>      The behavior is as follows:
!>          - pressure is always handled as a Dirichlet condition
!>          - if the mass flow is inflowing:
!>              one retains the velocity at infinity
!>              Dirichlet condition for scalars and turbulent values
!>               (or zero flux if the user has not specified a
!>                Dirichlet value)
!>          - if the mass flow is outflowing:
!>              one prescribes zero flux on the velocity, the scalars,
!>              and turbulent values
!>    .
!>    Note that the pressure will be reset to p0 on the first free outlet
!>    face found.
!>
!>
!> \subsection nonstd_bcs For "non-standard" conditions:
!>
!> Other than (inlet, free outlet, wall, symmetry), one defines
!>  - on one hand, for each face:
!>    - an admissible \c itypfb value (i.e. greater than or equal to 1 and
!>      less than or equal to \c ntypmx; see its value in paramx.h).
!>      The values predefined in paramx.h:
!>      \c ientre, \c isolib, \c isymet, \c iparoi, \c iparug are in this range,
!>      and it is preferable not to assign one of these integers to \c itypfb
!>      randomly or in an inconsiderate manner. To avoid this, one may use
!>      \c iindef if one wish to avoid checking values in paramx.h. \c iindef
!>      is an admissible value to which no predefined boundary condition
!>      is attached.
!>      Note that the \c itypfb array is reinitialized at each time step to
!>      the non-admissible value of 0. If one forgets to modify \c itypfb for
!>      a given face, the code will stop.
!>
!>  - and on the other hand, for each face and each variable:
!>    - a code
!>                      - \c icodcl(ifac, ivar)
!>    - three real values
!>                      - \c rcodcl(ifac, ivar, 1)
!>                      - \c rcodcl(ifac, ivar, 2)
!>                      - \c rcodcl(ifac, ivar, 3)
!>
!> \anchor icodcl \anchor rcodcl
!> The value of \c icodcl is taken from the following:
!>  - 1: Dirichlet      (usable for any variable)
!>  - 3: Neumann        (usable for any variable)
!>  - 4: Symmetry       (usable only for the velocity and components of
!>                       the Rij tensor)
!>  - 5: Smooth wall    (usable for any variable except for pressure)
!>  - 6: Rough wall     (usable for any variable except for pressure)
!>  - 9: Free outlet    (usable only for velocity)
!>  - 13: Dirichlet for the advection operator and
!>        Neumann for the diffusion operator
!>
!> The values of the 3 \c rcodcl components are:
!>  - \c rcodcl(ifac, ivar, 1):
!>     - Dirichlet for the variable          if \c icodcl(ifac, ivar) = 1 or 13
!>     - Wall value (sliding velocity, temp) if \c icodcl(ifac, ivar) = 5
!>     .
!>     The dimension of \c rcodcl(ifac, ivar, 1) is that of the
!>     resolved variable, for instance:
!>        - U (velocity in m/s),
!>        - T (temperature in degrees)
!>        - H (enthalpy in J/kg)
!>        - F (passive scalar in -)
!>  - \c rcodcl(ifac, ivar, 2):
!>       "exterior" exchange coefficient (between the prescribed value
!>                        and the value at the domain boundary)
!>                        rinfin = infinite by default
!>     - For velocities U,                in kg/(m2 s):
!>        \c rcodcl(ifac, ivar, 2) =          (viscl+visct) / d
!>     - For the pressure P,              in  s/m:
!>        \c rcodcl(ifac, ivar, 2) =                     dt / d
!>     - For temperatures T,              in Watt/(m2 degres):
!>        \c rcodcl(ifac, ivar, 2) = Cp*(viscls+visct/turb_schmidt) / d
!>     - For enthalpies H,                in kg /(m2 s):
!>        \c rcodcl(ifac, ivar, 2) =    (viscls+visct/turb_schmidt) / d
!>     - For other scalars F              in:
!>        \c rcodcl(ifac, ivar, 2) =    (viscls+visct/turb_schmidt) / d
!>            (d has the dimension of a distance in m)
!>
!>  - \c rcodcl(ifac, ivar, 3) if \c icodcl(ifac, ivar) = 3 or 13:
!>      Flux density (< 0 if gain, n outwards-facing normal)
!>     - For velocities U,                in kg/(m s2) = J:
!>        \c rcodcl(ifac, ivar, 3) =         -(viscl+visct) * (grad U).n
!>     - For pressure P,                  in kg/(m2 s):
!>        \c rcodcl(ifac, ivar, 3) =                    -dt * (grad P).n
!>     - For temperatures T,              in Watt/m2:
!>        \c rcodcl(ifac, ivar, 3) = -Cp*(viscls+visct/turb_schmidt) * (grad T).n
!>     - For enthalpies H,                in Watt/m2:
!>        \c rcodcl(ifac, ivar, 3) = -(viscls+visct/turb_schmidt) * (grad H).n
!>     - For other scalars F              in:
!>        \c rcodcl(ifac, ivar, 3) = -(viscls+visct/turb_schmidt) * (grad F).n
!>
!>  - \c rcodcl(ifac, ivar, 3) if \c icodcl(ifac, ivar) = 6:
!>      Roughness for the rough wall law
!>     - For velocities U, dynamic roughness
!>         \c rcodcl(ifac, iu, 3) = roughd
!>     - For other scalars, thermal roughness
!>         \c rcodcl(ifac, iv, 3) = rough
!>
!>
!> Note that if the user assigns a value to \c itypfb equal to \c ientre, \c isolib,
!> \c isymet, \c iparoi, or \c iparug and does not modify \c icodcl (zero value by
!>  default), \c itypfb will define the boundary condition type.
!>
!> To the contrary, if the user prescribes \c icodcl(ifac, ivar) (nonzero),
!> the values assigned to \c rcodcl will be used for the considered face
!> and variable (if \c rcodcl values are not set, the default values will
!> be used for the face and variable, so:
!>                          - \c rcodcl(ifac, ivar, 1) = 0.d0
!>                          - \c rcodcl(ifac, ivar, 2) = rinfin
!>                          - \c rcodcl(ifac, ivar, 3) = 0.d0)
!>
!> Especially, one may have for example:
!>  - set \c itypfb(ifac) = \c iparoi which prescribes default wall
!>    conditions for all variables at face ifac,
!>  - and define IN ADDITION for variable ivar on this face specific
!>    conditions by specifying \c icodcl(ifac, ivar) and the 3 \c rcodcl values.
!>
!> The user may also assign to \c itypfb a value not equal to \c ientre, \c isolib,
!> \c isymet, \c iparoi, \c iparug, \c iindef but greater than or equal to 1 and less
!> than or equal to ntypmx (see values in param.h) to distinguish groups
!> or colors in other subroutines which are specific to the case and in
!> which itypfb is accessible. In this case though it will be necessary
!> to prescribe boundary conditions by assigning values to icodcl and to
!> the 3 \c rcodcl fields (as the value of \c itypfb will not be predefined in
!> the code).
!>
!>
!> \subsection comp_bcs Boundary condition types for compressible flows
!>
!> For compressible flows, only predefined boundary conditions may
!> be assigned among: \c iparoi, \c isymet, \c iesicf, \c isspcf, \c isopcf, \c iephcf, \c ieqhcf
!>
!>  - \c iparoi : standard wall
!>  - \c isymet : standard symmetry
!>
!>  - \c iesicf, \c isspcf, \c isopcf, \c iephcf, \c ieqhcf : inlet/outlet
!>
!> For inlets/outlets, we can prescribe
!> a value for turbulence and passive scalars in \c rcodcl(.,.,1)
!> for the case in which the mass flux is incoming. If this is not
!> done, a zero flux condition is applied.
!>
!> - \c iesicf: prescribed inlet/outlet (for example supersonic inlet)
!>           the user prescribes the velocity and all thermodynamic variables
!> - \c isspcf: supersonic outlet
!>           the user does not prescribe anything
!> - \c isopcf: subsonic outlet with prescribed pressure
!>           the user presribes the pressure
!> - \c iephcf: mixed inlet with prescribed total pressure and enthalpy
!>           the user prescribes the total pressure and total enthalpy
!> - \c ieqhcf: subsonic inlet with prescribed mass and enthalpy flow
!>           to be implemented
!>
!>
!> \subsection cons_rul Consistency rules
!>
!> A few consistency rules between \c icodcl codes for variables with
!> non-standard boundary conditions:
!>
!>  - Codes for velocity components must be identical
!>  - Codes for Rij components must be identical
!>  - If code (velocity or Rij) = 4
!>    one must have code (velocity and Rij) = 4
!>  - If code (velocity or turbulence) = 5
!>    one must have code (velocity and turbulence) = 5
!>  - If code (velocity or turbulence) = 6
!>    one must have code (velocity and turbulence) = 6
!>  - If scalar code (except pressure or fluctuations) = 5
!>    one must have velocity code = 5
!>  - If scalar code (except pressure or fluctuations) = 6
!>    one must have velocity code = 6
!>
!>
!> \remarks
!>   - Caution: to prescribe a flux (nonzero) to Rij, the viscosity to take
!>              into account is viscl even if visct exists
!>              (visct=rho cmu k2/epsilon)
!>   - One have the ordering array for boundary faces from the previous time
!>       step (except for the fist one, where \c itrifb has not been set yet).
!>   - The array of boundary face types \c itypfb has been reset before
!>       entering the subroutine.
!>
!>
!> \subsubsection cs_user_bc_cell_id Cell values of some variables
!>
!> Cell value field ids
!>
!> - Density:                        \c irom
!> - Dynamic molecular viscosity:    \c iviscl
!> - Turbulent viscosity:            \c ivisct
!> - Specific heat:                  \c icp
!> - Diffusivity(lambda):            \c field_get_key_int(ivarfl(isca(iscal)), &
!>                                      kivisl, ...)
!>
!>
!> \subsubsection fac_id Faces identification
!>
!> - Density:                               \c field id \c ibrom
!> - Boundary mass flux (for convecting \c ivar):
!>     field id \c iflmab
!>     using \c field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
!> - For other values: take as an approximation the value in the adjacent cell
!>                     i.e. as above with \c iel = ifabor(ifac).
!>
!> Please refer to the
!> <a href="../../theory.pdf#boundary"><b>boundary conditions</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atchem
use atincl
use atsoil
use ctincl
use cs_fuel_incl
use mesh
use field
use turbomachinery
use iso_c_binding
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer, allocatable, dimension(:) :: lstelt

integer          ielt, nlelt, ifac

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

call getfbr('1', nlelt, lstelt)

if (ttcabs.lt.3.8d0) then
 do ielt = 1, nlelt
    ifac = lstelt(ielt)
    rcodcl(ifac, isca(2),1) = 20.d0 + 100.d0*ttcabs
 enddo
else
 do ielt = 1, nlelt
    ifac = lstelt(ielt)
    rcodcl(ifac, isca(2),1) = 400.d0
 enddo
endif

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions

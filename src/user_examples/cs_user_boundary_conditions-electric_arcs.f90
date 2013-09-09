!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine cs_user_boundary_conditions &
!=====================================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   rcodcl )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Fill boundary conditions arrays (icodcl, rcodcl) for unknown variables.


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be selected using the 'getfbr' subroutine.

!  getfbr(string, nelts, eltlst):
!  - string is a user-supplied character string containing selection criteria;
!  - nelts is set by the subroutine. It is an integer value corresponding to
!    the number of boundary faces verifying the selection criteria;
!  - lstelt is set by the subroutine. It is an integer array of size nelts
!    containing the list of boundary faces verifying the selection criteria.

!  string may contain:
!  - references to colors (ex.: 1, 8, 26, ...)
!  - references to groups (ex.: inlet, group1, ...)
!  - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!  These criteria may be combined using logical operators ('and', 'or') and
!  parentheses.
!  Example: '1 and (group2 or group3) and y < 1' will select boundary faces
!  of color 1, belonging to groups 'group2' or 'group3' and with face center
!  coordinate y less than 1.

!  Operators priority, from highest to lowest:
!    '( )' > 'not' > 'and' > 'or' > 'xor'

! Similarly, interior faces and cells can be identified using the 'getfac'
! and 'getcel' subroutines (respectively). Their syntax are identical to
! 'getfbr' syntax.

! For a more thorough description of the criteria syntax, see the user guide.


! Boundary condition types
! ========================

! Boundary conditions may be assigned in two ways.


!    For "standard" boundary conditions:
!    -----------------------------------

!     (inlet, free outlet, wall, symmetry), one defines a code in the 'itypfb'
!     array (of dimensions number of boundary faces).
!     This code will then be used by a non-user subroutine to assign the
!     following conditions.
!     Thus:

!     Code      |  Boundary type
!     --------------------------
!      ientre   |   Inlet
!      isolib   |   Free outlet
!      isymet   |   Symmetry
!      iparoi   |   Wall (smooth)
!      iparug   |   Rough wall

!     These integers are defined elsewhere (in paramx.f90 module).
!     Their value is greater than or equal to 1 and less than or  equal to
!     ntypmx (value fixed in paramx.h)


!     In addition, some values must be defined:


!     - Inlet (more precisely, inlet/outlet with prescribed flow, as the flow
!              may be prescribed as an outflow):

!       -> Dirichlet conditions on variables other than pressure are mandatory
!         if the flow is incoming, optional if the flow is outgoing (the code
!         assigns zero flux if no Dirichlet is specified); thus,
!         at face 'ifac', for the variable 'ivar': rcodcl(ifac, ivar, 1)


!     - Smooth wall: (= impermeable solid, with smooth friction)

!       -> Velocity value for sliding wall if applicable
!         at face ifac, rcodcl(ifac, iu, 1)
!                       rcodcl(ifac, iv, 1)
!                       rcodcl(ifac, iw, 1)
!       -> Specific code and prescribed temperature value at wall if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 5
!                       rcodcl(ifac, ivar, 1) = prescribed temperature
!       -> Specific code and prescribed flux value at wall if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 3
!                       rcodcl(ifac, ivar, 3) = prescribed flux

!        Note that the default condition for scalars (other than k and epsilon)
!        is homogeneous Neumann.


!     - Rough wall: (= impermeable solid, with rough friction)

!       -> Velocity value for sliding wall if applicable
!         at face ifac, rcodcl(ifac, iu, 1)
!                       rcodcl(ifac, iv, 1)
!                       rcodcl(ifac, iw, 1)
!       -> Value of the dynamic roughness height to specify in
!                       rcodcl(ifac, iu, 3)
!       -> Value of the scalar roughness height (if required) to specify in
!                       rcodcl(ifac, iv, 3) (values for iw are not used)
!       -> Specific code and prescribed temperature value at wall if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 6
!                       rcodcl(ifac, ivar, 1) = prescribed temperature
!       -> Specific code and prescribed flux value at rough wall, if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 3
!                       rcodcl(ifac, ivar, 3) = prescribed flux

!        Note that the default condition for scalars (other than k and epsilon)
!        is homogeneous Neumann.

!     - Symmetry (= slip wall):

!       -> Nothing to specify


!     - Free outlet (more precisely free inlet/outlet with prescribed pressure)

!       -> Nothing to prescribe for pressure and velocity. For scalars and
!          turbulent values, a Dirichlet value may optionally be specified.
!          The behavior is as follows:
!              * pressure is always handled as a Dirichlet condition
!              * if the mass flow is inflowing:
!                  one retains the velocity at infinity
!                  Dirichlet condition for scalars and turbulent values
!                   (or zero flux if the user has not specified a
!                    Dirichlet value)
!                if the mass flow is outflowing:
!                  one prescribes zero flux on the velocity, the scalars,
!                  and turbulent values

!       Note that the pressure will be reset to p0 on the first free outlet
!       face found


!    For "non-standard" conditions:
!    ------------------------------

!     Other than (inlet, free outlet, wall, symmetry), one defines
!      - on one hand, for each face:
!        -> an admissible 'itypfb' value (i.e. greater than or equal to 1 and
!           less than or equal to ntypmx; see its value in paramx.h).
!           The values predefined in paramx.h:
!           'ientre', 'isolib', 'isymet', 'iparoi', 'iparug' are in this range,
!           and it is preferable not to assign one of these integers to 'itypfb'
!           randomly or in an inconsiderate manner. To avoid this, one may use
!           'iindef' if one wish to avoid checking values in paramx.h. 'iindef'
!           is an admissible value to which no predefined boundary condition
!           is attached.
!           Note that the 'itypfb' array is reinitialized at each time step to
!           the non-admissible value of 0. If one forgets to modify 'typfb' for
!           a given face, the code will stop.

!      - and on the other hand, for each face and each variable:
!        -> a code             icodcl(ifac, ivar)
!        -> three real values  rcodcl(ifac, ivar, 1)
!                              rcodcl(ifac, ivar, 2)
!                              rcodcl(ifac, ivar, 3)
!     The value of 'icodcl' is taken from the following:
!       1: Dirichlet      (usable for any variable)
!       3: Neumann        (usable for any variable)
!       4: Symmetry       (usable only for the velocity and components of
!                          the Rij tensor)
!       5: Smooth wall    (usable for any variable except for pressure)
!       6: Rough wall     (usable for any variable except for pressure)
!       9: Free outlet    (usable only for velocity)
!     The values of the 3 'rcodcl' components are:
!      rcodcl(ifac, ivar, 1):
!         Dirichlet for the variable          if icodcl(ifac, ivar) =  1
!         Wall value (sliding velocity, temp) if icodcl(ifac, ivar) =  5
!         The dimension of rcodcl(ifac, ivar, 1) is that of the
!           resolved variable: ex U (velocity in m/s),
!                                 T (temperature in degrees)
!                                 H (enthalpy in J/kg)
!                                 F (passive scalar in -)
!      rcodcl(ifac, ivar, 2):
!         "exterior" exchange coefficient (between the prescribed value
!                          and the value at the domain boundary)
!                          rinfin = infinite by default
!         For velocities U,                in kg/(m2 s):
!           rcodcl(ifac, ivar, 2) =          (viscl+visct) / d
!         For the pressure P,              in  s/m:
!           rcodcl(ifac, ivar, 2) =                     dt / d
!         For temperatures T,              in Watt/(m2 degres):
!           rcodcl(ifac, ivar, 2) = Cp*(viscls+visct/sigmas) / d
!         For enthalpies H,                in kg /(m2 s):
!           rcodcl(ifac, ivar, 2) =    (viscls+visct/sigmas) / d
!         For other scalars F              in:
!           rcodcl(ifac, ivar, 2) =    (viscls+visct/sigmas) / d
!              (d has the dimension of a distance in m)
!
!      rcodcl(ifac, ivar, 3) if icodcl(ifac, ivar) <> 6:
!        Flux density (< 0 if gain, n outwards-facing normal)
!                            if icodcl(ifac, ivar)= 3
!         For velocities U,                in kg/(m s2) = J:
!           rcodcl(ifac, ivar, 3) =         -(viscl+visct) * (grad U).n
!         For pressure P,                  in kg/(m2 s):
!           rcodcl(ifac, ivar, 3) =                    -dt * (grad P).n
!         For temperatures T,              in Watt/m2:
!           rcodcl(ifac, ivar, 3) = -Cp*(viscls+visct/sigmas) * (grad T).n
!         For enthalpies H,                in Watt/m2:
!           rcodcl(ifac, ivar, 3) = -(viscls+visct/sigmas) * (grad H).n
!         For other scalars F              in:
!           rcodcl(ifac, ivar, 3) = -(viscls+visct/sigmas) * (grad F).n

!      rcodcl(ifac, ivar, 3) if icodcl(ifac, ivar) = 6:
!        Roughness for the rough wall law
!         For velocities U, dynamic roughness
!           rcodcl(ifac, iu, 3) = roughd
!         For other scalars, thermal roughness
!           rcodcl(ifac, iv, 3) = rought


!      Note that if the user assigns a value to itypfb equal to ientre, isolib,
!       isymet, iparoi, or iparug and does not modify icodcl (zero value by
!       default), itypfb will define the boundary condition type.

!      To the contrary, if the user prescribes icodcl(ifac, ivar) (nonzero),
!        the values assigned to rcodcl will be used for the considered face
!        and variable (if rcodcl values are not set, the default values will
!        be used for the face and variable, so:
!                                 rcodcl(ifac, ivar, 1) = 0.d0
!                                 rcodcl(ifac, ivar, 2) = rinfin
!                                 rcodcl(ifac, ivar, 3) = 0.d0)
!        Especially, one may have for example:
!        -> set itypfb(ifac) = iparoi which prescribes default wall
!        conditions for all variables at face ifac,
!        -> and define IN ADDITION for variable ivar on this face specific
!        conditions by specifying icodcl(ifac, ivar) and the 3 rcodcl values.


!      The user may also assign to itypfb a value not equal to ientre, isolib,
!       isymet, iparoi, iparug, iindef but greater than or equal to 1 and less
!       than or equal to ntypmx (see values in param.h) to distinguish groups
!       or colors in other subroutines which are specific to the case and in
!       which itypfb is accessible.  In this case though it will be necessary
!       to prescribe boundary conditions by assigning values to icodcl and to
!       the 3 rcodcl fields (as the value of itypfb will not be predefined in
!       the code).


! Boundary condition types for compressible flows
! ===============================================

! For compressible flows, only predefined boundary conditions may
! be assigned

!    iparoi, isymet, iesicf, isspcf, isopcf, ierucf, ieqhcf

!    iparoi : standard wall
!    isymet : standard symmetry

!    iesicf, isspcf, isopcf, ierucf, ieqhcf : inlet/outlet

! For inlets/outlets, we can prescribe
!  a value for turbulence and passive scalars in rcodcl(.,.,1)
!  for the case in which the mass flux is incoming. If this is not
!  done, a zero flux condition is applied.

! iesicf : prescribed inlet/outlet (for example supersonic inlet)
!         the user prescribes the velocity and all thermodynamic variables
! isspcf : supersonic outlet
!         the user does not prescribe anything
! isopcf : subsonic outlet with prescribed pressure
!         the user presribes the pressure
! ierucf : subsonic inlet with prescribed velocity and density
!         the user prescribes the velocity and density
! ieqhcf : subsonic inlet with prescribed mass and enthalpy flow
!         to be implemented


! Consistency rules
! =================

!       A few consistency rules between 'icodcl' codes for variables with
!       non-standard boundary conditions:

!           Codes for velocity components must be identical
!           Codes for Rij components must be identical
!           If code (velocity or Rij) = 4
!             one must have code (velocity and Rij) = 4
!           If code (velocity or turbulence) = 5
!             one must have code (velocity and turbulence) = 5
!           If code (velocity or turbulence) = 6
!             one must have code (velocity and turbulence) = 6
!           If scalar code (except pressure or fluctuations) = 5
!             one must have velocity code = 5
!           If scalar code (except pressure or fluctuations) = 6
!             one must have velocity code = 6


! Remarks
! =======

!       Caution: to prescribe a flux (nonzero) to Rij, the viscosity to take
!                into account is viscl even if visct exists
!                (visct=rho cmu k2/epsilon)

!       One have the ordering array for boundary faces from the previous time
!         step (except for the fist one, where 'itrifb' has not been set yet).
!       The array of boundary face types 'itypfb' has been reset before
!         entering the subroutine.


!       Note how to access some variables (for variable 'ivar'
!                                              scalar   'iscal'):

! Cell values  (let iel = ifabor(ifac))

! * Density:                                 propce(iel, ipproc(irom))
! * Dynamic molecular viscosity:             propce(iel, ipproc(iviscl))
! * Turbulent viscosity:                     propce(iel, ipproc(ivisct))
! * Specific heat:                           propce(iel, ipproc(icp)
! * Diffusivity(lambda):                     propce(iel, ipproc(ivisls(iscal)))

! Boundary face values

! * Density:                                 propfb(ifac, ipprob(irom))
! * Mass flux (for convecting 'ivar'):
!     field id 'iflmab', using field_get_key_int(ivarfl(ivar), kimasf, iflmab)

! * For other values: take as an approximation the value in the adjacent cell
!                     i.e. as above with iel = ifabor(ifac).


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb(nfabor)   ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb(nfabor)   ! ia ! --> ! boundary face types                            !
! izfppp(nfabor)   ! ia ! --> ! boundary face zone number                      !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!  (nfabor,nvar,3) !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use ctincl
use elincl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(nfabor,*)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

!< [loc_var_dec]
integer          ifac, ii, iel
integer          idim
integer          izone,iesp
integer          ilelt, nlelt

double precision uref2, d2s3
double precision rhomoy, dhy, ustar2
double precision xkent, xeent
double precision z1   , z2

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! --- For boundary faces of color 1 assign an inlet for all phases
!     ============================================================
!        and assign a cathode for "electric" variables.
!        =============================================
!
!< [example_1]
call getfbr('1', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac) = ientre

  ! Zone Number (from 1 to n)
  izone = 1

  ! Zone localization for a given face
  izfppp(ifac) = izone

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 0.d0

  ! Turbulence

  if (itytur.eq.2 .or. itytur.eq.3                  &
       .or. iturb.eq.50 .or. iturb.eq.60            &
       .or. iturb.eq.70) then

    uref2 = rcodcl(ifac,iu,1)**2                    &
           +rcodcl(ifac,iv,1)**2                    &
           +rcodcl(ifac,iw,1)**2
    uref2 = max(uref2,1.d-12)

    ! Turbulence example computed using equations valid for a pipe.

    ! We will be careful to specify a hydraulic diameter adapted
    ! to the current inlet.

    ! We will also be careful if necessary to use a more precise
    ! formula for the dynamic viscosity use in the calculation of
    ! the Reynolds number (especially if it is variable, it may be
    ! useful to take the law from 'usphyv'. Here, we use by default
    ! the 'viscl0" value.

    ! Regarding the density, we have acess to its value at boundary
    ! faces (romb) so this value is the one used here (specifically,
    ! it is consistent with the processing in 'usphyv', in case of
    ! variable density)

    ! Hydraulic diameter
    dhy     = 0.075d0

    ! Calculation of friction velocity squared (ustar2)
    ! and of k and epsilon at the inlet (xkent and xeent) using
    ! standard laws for a circular pipe
    ! (their initialization is not needed here but is good practice).

    rhomoy = propfb(ifac,ipprob(irom))
    ustar2 = 0.d0
    xkent  = epzero
    xeent  = epzero

    call keendb                                            &
    !==========
     ( uref2, dhy, rhomoy, viscl0, cmu, xkappa,            &
       ustar2, xkent, xeent)

    if (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = xkent
      rcodcl(ifac,iep,1) = xeent

    elseif (itytur.eq.3) then

      rcodcl(ifac,ir11,1) = d2s3*xkent
      rcodcl(ifac,ir22,1) = d2s3*xkent
      rcodcl(ifac,ir33,1) = d2s3*xkent
      rcodcl(ifac,ir12,1) = 0.d0
      rcodcl(ifac,ir13,1) = 0.d0
      rcodcl(ifac,ir23,1) = 0.d0
      rcodcl(ifac,iep,1)  = xeent

    elseif (iturb.eq.50) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iep,1)  = xeent
      rcodcl(ifac,iphi,1) = d2s3
      rcodcl(ifac,ifb,1)  = 0.d0

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

    endif

  endif

  ! --- Handle Scalars

  ! Enthalpy in J/kg (ihm)
  ! On this example we impose the value of the enthalpy
  ! the arbitrary value of 1.d6 corresponds to a temperature of 2200 Kelvin
  ! for argon at atmospheric pressure (see dp_ELE)

  ii = ihm
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 1.d6

  !  Electric potential  (ipotr)
  ! (could corresponds also to the real part of the electrical potential if
  ! Joule Effect by direct conduction)

  ! In the Cathode example (electric arc applications),
  ! we impose a constant value of the electrical potential which is zero,
  ! assuming that the potential is equal to "ipotr + an arbitrary constant"
  ! (What is important for electric arc is the difference between anode and
  ! cathode potentials)

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 0.d0

  ! Mass fraction of the (n-1) gas mixture components

  if (ngazg .gt. 1) then
    do iesp=1,ngazg-1
      ii = iycoel(iesp)
      icodcl(ifac,isca(ii))   = 1
      rcodcl(ifac,isca(ii),1) = 0.d0
    enddo
  endif

  ! Specific model for Joule effect by direct conduction:
  ! Imaginary part of the potentiel (ipoti) is imposed to zero

  if (ippmod(ieljou).ge. 2) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    rcodcl(ifac,isca(ii),1) = 0.d0
  endif

  ! Specific model for Electric arc:
  ! Vector Potential: Zero flux by default beacuse we don't a lot about
  ! vector potential (what we know, is that A is equal to zero at the infinite)

  ! All the boundary conditions for A are zero flux, except on some chosen faces
  ! where we need to impose a value in order to have a stable calculation
  ! (well defined problem)
  ! These faces are chosen where we are sure that the electrical current density
  ! remains very low generally far from the center of the electric arc and from
  ! the electrodes (see above)

  if (ippmod(ielarc).ge.2) then
    do idim= 1,ndimve
      ii = ipotva(idim)
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0
    enddo
  endif

enddo
!< [example_1]

! --- For boundary faces of color 5 assign an free outlet for all phases
!     ==================================================================
!        and example of electrode for Joule Effect by direct conduction.
!        ==============================================================
!
!< [example_2]
call getfbr('5', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac)   = isolib

  ! Zone Number (from 1 to n)
  izone = 2

  ! Zone location for a given face

  izfppp(ifac) = izone

  ! --- Handle Scalars

  ! Enthalpy in J/kg  (By default zero flux with ISOLIB)
  !    Nothing to do

  ! Mass fraction of the (n-1) gas mixture components
  ! (Zero flux by defaut with ISOLIB)
  !    Nothing to do

  ! Specific model for Joule Effect by direct conduction:

  ! If you want to make a simulation with an imposed Power PUISIM
  ! (you want to get PUISIM imposed in useli1 and PUISIM = Amp x Volt)
  ! you need to impose IELCOR=1 in useli1
  ! The boundary conditions will be scaled by COEJOU coefficient
  ! for example the electrical potential will be multiplied bu COEJOU
  ! (both real and imaginary part of the electrical potential if needed)

  ! COEJOU is automatically defined in order that the calculated dissipated power
  ! by Joule effect (both real and imaginary part if needed) is equal to PUISIM

  ! At the beginning of the calculation, COEJOU ie equal to 1;
  ! COEJOU is writing and reading in the result files.

  ! If you don't want to calculate with by scaling,
  ! you can impose directly the value.

  if (ippmod(ieljou).ge. 1) then
    ii = ipotr
    icodcl(ifac,isca(ii))   = 1
    if (ielcor.eq.1) then
      rcodcl(ifac,isca(ii),1) = 500.d0*coejou
    else
      rcodcl(ifac,isca(ii),1) = 500.d0
    endif
  endif

  if (ippmod(ieljou).ge. 2) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    if (ielcor.eq.1) then
      rcodcl(ifac,isca(ii),1) = sqrt(3.d0)*500.d0*coejou
    else
      rcodcl(ifac,isca(ii),1) = sqrt(3.d0)*500.d0
    endif
  endif

enddo
!< [example_2]

! --- For boundary faces of color 2 assign a free outlet for all phases
!     =================================================================
!        and example of anode for electric arc.
!        =====================================
!
!< [example_3]
call getfbr('2', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac)   = isolib

  ! Zone number (from 1 to n)
  izone = 3

  ! Zone localization for a given face
  izfppp(ifac) = izone

  ! --- Handle scalars

  ! Enthalpy in J/kg  (Zero flux by default with ISOLIB)
  !    Nothing to do

  ! Real component of the electrical potential

  ! For electric arc model,
  ! ======================
  ! *  we generally calculate the "electric variables" assuming that the total
  !    intensity of the electrical current is imposed (COUIMP is the value of
  !    the imposed total current).

  !    In that case, you need to impose IELCOR=1 in useli1
  !    The "electrical variables" will be scaled by COEPOT coefficient :
  !    for example the electrical potential will be multiplied by COEPOT,
  !    Joule effect will be multipied by COEPOT * COEPOT and so on (see uselrc.f90)

  !    COEJOU is defined in uselrc.fr : different possibilities are described in
  !    uselrc.f90, depending on the different physics you want to simulate
  !    (scaling from current, from power, special model for restriking ...)

  !    The variable DPOT is defined: it corresponds to the electrical potential
  !    difference between the electrodes (Anode potential - cathode Potential).
  !    DPOT is calculated in uselrc.f90. DPOT is saved at each time step, and
  !    for a following calculation

  !    DPOT is the value of the boundary condition on anode assuming that
  !    the cathode potential is equel to zero.

  ! *  It is also possible to fix the value of the potential on the anode.
  !    (for example, 1000 Volts).

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1

  if (ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    rcodcl(ifac,isca(ii),1) = dpot
  else
    rcodcl(ifac,isca(ii),1) = 1000.d0
  endif

  ! Mass fraction of the (n-1) gas mixture components
  !  zero flux by default with ISOLIB
  !  nothing to do

  ! vector Potential
  !  zero flux by default with ISOLIB
  !  nothing to do

enddo
!< [example_3]

! --- For boundary faces of color 3 assign a wall for all phases
!     ==========================================================
!        and example of potential vector Dirichlet condition
!        ===================================================

!< [example_4]
call getfbr('3', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac)   = iparoi

  ! Zone number (from 1 to n)
  izone = 4

  ! Zone localization for a given face
  izfppp(ifac) = izone

  ! Wall: zero flow (zero flux for pressure)
  !       friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  ! --- Handle scalars
  !  Enthalpy in J/kg  (Zero flux by default)
  !     Nothing to do

  ! Real component of the electrical potential
  !    Zero flux by default
  !    Nothing to do

  ! Specific model for Electric arc :
  ! ================================

  ! Vector potential  A (Ax, Ay, Az)

  ! Zero flux by default because we don't a lot about vector potential
  ! (what we know, is that A is equal to zero at the infinite)

  ! All the boundary conditions for A are zero flux, except on some chosen faces
  ! where we need to impose a value in order to have a stable calculation
  ! These faces are chosen where we are sure that the electrical current density
  ! remains very low generally far from the center of the electric arc and from
  ! the electrodes:

  ! On the following example, we choose to impose a "dirichlet" value for the
  ! 3 components of A on a small zone of the boundary located near the vertical
  ! free outlet of the computation domain.

  ! In this example, the electric arc is at the center of the computational domain,
  ! located on z axis (near x = 0 and y = 0).
  ! The x (1st) and y (the 3rd) coordinates are contained between
  ! -2.5 cm nd 2.5 cm:

  !    Ax(t, x,y,z) = Ax(t-dt, x=2.5cm, y=2.5cm, z)
  !    Ay(t, x,y,z) = Ay(t-dt, x=2.5cm, y=2.5cm, z)
  !    Az(t, x,y,z) = Az(t-dt, x=2.5cm, y=2.5cm, z)

  if (ippmod(ielarc).ge.2) then
    if (cdgfbo(1,ifac) .le.  2.249d-2  .or.                      &
        cdgfbo(1,ifac) .ge.  2.249d-2  .or.                      &
        cdgfbo(3,ifac) .le. -2.249d-2  .or.                      &
        cdgfbo(3,ifac) .ge.  2.249d-2      ) then
      iel = ifabor(ifac)
      do idim = 1, ndimve
        ii = ipotva(idim)
        icodcl(ifac,isca(ii))   = 1
        rcodcl(ifac,isca(ii),1) = rtpa(iel,isca(ii))
      enddo
    endif
  endif

enddo
!< [example_4]

! --- For boundary faces of color 51 assign a wall
!     ============================================
!        and restriking model for electric arc (anode boundaray condition)
!        =================================================================

!< [example_5]
call getfbr('51', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac)   = iparoi

  ! Zone number (from 1 to n)
  izone = 5

  ! Zone localization for a given face
  izfppp(ifac) = izone

  ! ---- Enthalpy (J/kg) :
  !      imposed heat transfer coefficient

  ii=ihm
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 2.d4
  rcodcl(ifac,isca(ii),2) = 1.d5

  ! Real electrical potential: anode boundary condition;
  ! dpot calculated in uselrc.f

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1

  if (ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    rcodcl(ifac,isca(ii),1) = dpot
  else
    rcodcl(ifac,isca(ii),1) = 100.d0
  endif

  ! Restriking modeling:
  ! ===================
  !    example to fit depending on the case, the geometry etc...
  !     and also in agreement with uselrc.fr

  if (ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    if (iclaq.eq.1 .and. ntcabs.le.ntdcla+30) then

      z1 = zclaq - 2.d-4
      if (z1.le.0.d0) z1 = 0.d0
      z2 = zclaq + 2.d-4
      if (z2.ge.2.d-2) z2 = 2.d-2

      if (cdgfbo(3,ifac).ge.z1 .and. cdgfbo(3,ifac).le.z2) then
        icodcl(ifac,isca(ii))   = 1
        rcodcl(ifac,isca(ii),1) = dpot
      else
        icodcl(ifac,isca(ii))   = 3
        rcodcl(ifac,isca(ii),3) = 0.d0
      endif
    endif
  endif

  ! Vector potential : Zero flux

  if (ippmod(ielarc).ge.2) then
    do idim= 1,ndimve
      ii = ipotva(idim)
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0
    enddo
  endif

enddo
!< [example_5]

! --- For boundary faces of color 4 assign a symetry
!     ==============================================

!< [example_6]
call getfbr('4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Symmetries

  itypfb(ifac)   = isymet

  ! Zone number (from 1 to n)
  izone = 6

  ! Zone localization for a given face
  izfppp(ifac) = izone

  ! For all scalars, by default a zero flux condition is assumed
  ! (for potentials also)

  ! In Joule effect direct conduction,
  ! we can use an anti-symetry condition for the imaginary component of the
  ! electrical potential depending on the electrode configuration:

  if (ippmod(ieljou).ge. 2) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    rcodcl(ifac,isca(ii),1) = 0.d0
  endif

enddo
!< [example_6]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_user_boundary_conditions

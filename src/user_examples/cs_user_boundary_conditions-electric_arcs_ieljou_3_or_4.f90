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

integer          ifac  , ii     , iel
integer          i     , ntf    , nb    , id , itrouv
integer          izone
integer          nborne(nbtrmx)
integer          ilelt , nlelt

double precision rnbs2,capaeq
double precision sir(nelemx)   ,sii(nelemx)
double precision sirb(nbtrmx,6),siib(nbtrmx,6)
double precision ur(nbtrmx,6)  ,ui(nbtrmx,6)
double precision sirt(nbtrmx)  ,siit(nbtrmx)
character*200    chain

integer, allocatable, dimension(:) :: lstelt

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

! 2.1 - Computation of intensity (A/m2) for each electrode
!       -----------------------------------------------------

!   Pre - initialisation

do i= 1,nbelec
  sir(i) = 0.d0
  sii(i) = 0.d0
enddo

do ntf= 1,nbtrf
  sirt(ntf) = 0.d0
  siit(ntf) = 0.d0
enddo

if (ntcabs.lt.(ntpabs+2)) then
  do ntf = 1,nbtrf
    uroff(ntf) = 0.d0
    uioff(ntf) = 0.d0
  enddo
endif

!     Loop on selected boundary faces

do i=1,nbelec

  chain = ' '
  write(chain,3000) ielecc(i)

  if ( ielect(i).ne. 0 ) then
    call getfbr(chain,nlelt,lstelt)
    !==========

    do ilelt = 1, nlelt

      ifac = lstelt(ilelt)

      iel = ifabor(ifac)

      do id=1,ndimve
        sir(i) = sir(i)                                           &
             +propce(iel,ipproc(idjr(id)))*surfbo(id,ifac)
      enddo

      if ( ippmod(ieljou) .eq. 4 ) then
        do id=1,ndimve
          sii(i) = sii(i)                                         &
               +propce(iel,ipproc(idji(id)))*surfbo(id,ifac)
        enddo
      endif

    enddo

  endif

enddo


! 2.2 - Definition of Voltage on each termin of transformers
!       ----------------------------------------------------------------

!  2.2.1 Computation of Intensity on each termin of transformers

do i=1,nbelec
  sirb(ielect(i),ielecb(i)) = 0.d0
  if ( ippmod(ieljou) .eq. 4 ) then
    siib(ielect(i),ielecb(i)) = 0.d0
  endif
enddo

do i=1,nbelec
  if ( ielect(i).ne. 0 ) then
    sirb(ielect(i),ielecb(i)) = sirb(ielect(i),ielecb(i))         &
                               +sir(i)
    if ( ippmod(ieljou) .eq. 4 ) then
       siib(ielect(i),ielecb(i)) = siib(ielect(i),ielecb(i))      &
                                  +sii(i)
    endif
  endif
enddo

!  2.2.2 RVoltage on each termin

do ntf=1,nbtrf

!      Primary and Secondary in Triangle

  if (ibrpr(ntf) .eq. 0 .and. ibrsec(ntf) .eq. 0 ) then

    nborne(ntf) = 3

    rnbs2 = 3.d0*rnbs(ntf)*rnbs(ntf)
    ur(ntf,1)=  1.154675d0*tenspr(ntf)/rnbs(ntf)                  &
      + (zr(ntf)*sirb(ntf,1)-zi(ntf)*siib(ntf,1))/rnbs2

    ur(ntf,2)= -0.5773d0*tenspr(ntf)/rnbs(ntf)                    &
      + (zr(ntf)*sirb(ntf,2)-zi(ntf)*siib(ntf,2))/rnbs2
    ur(ntf,3)= -0.5773d0*tenspr(ntf)/rnbs(ntf)                    &
      + (zr(ntf)*sirb(ntf,3)-zi(ntf)*siib(ntf,3))/rnbs2

    ui(ntf,1)=  0.d0                                              &
      + (zi(ntf)*sirb(ntf,1)+zr(ntf)*siib(ntf,1))/rnbs2
    ui(ntf,2)= -1.d0*tenspr(ntf)/rnbs(ntf)                        &
      + (zi(ntf)*sirb(ntf,2)+zr(ntf)*siib(ntf,2))/rnbs2
    ui(ntf,3)=  1.d0*tenspr(ntf)/rnbs(ntf)                        &
      + (zi(ntf)*sirb(ntf,3)+zr(ntf)*siib(ntf,3))/rnbs2

  else

    write(nfecra, *) 'Matrice sur le Transfo a ecrire'
    call csexit(1)

  endif
enddo

!  2.2.3 Total intensity for a transformer
!         (zero valued WHEN Offset established)

do ntf=1,nbtrf
  sirt(ntf) = 0.d0
  if ( ippmod(ieljou) .eq. 4 ) then
    siit(ntf) = 0.d0
  endif
enddo

do i=1,nbelec
  if ( ielect(i).ne. 0 ) then
    sirt(ielect(i)) = sirt(ielect(i)) + sir(i)
    if ( ippmod(ieljou) .eq. 4 ) then
      siit(ielect(i)) = siit(ielect(i)) + sii(i)
    endif
  endif
enddo

!  2.2.4 Take in account of Offset

capaeq = 3.d0

do ntf=1,nbtrf
  uroff(ntf) = uroff(ntf) + sirt(ntf)/capaeq
  if ( ippmod(ieljou) .eq. 4 ) then
    uioff(ntf) = uioff(ntf) + siit(ntf)/capaeq
  endif
enddo

! A reference transformer is assumed to have an Offset zero valued

if ( ntfref .gt. 0 ) then
  uroff(ntfref) = 0.d0
  uioff(ntfref) = 0.d0
endif

do ntf=1,nbtrf
  do nb=1,nborne(ntf)
    ur(ntf,nb) = ur(ntf,nb) + uroff(ntf)
    if ( ippmod(ieljou) .eq. 4 ) then
      ui(ntf,nb) = ui(ntf,nb) + uioff(ntf)
    endif
  enddo
enddo

! Print of UROFF (real part of offset potential)

write(nfecra,1500)
do ntf=1,nbtrf
  write(nfecra,2000) ntf,uroff(ntf)
enddo
write(nfecra,1501)

!  2.2.5 Take in account of Boundary Conditions


!     Loop on selected Boundary Faces

do i=1,nbelec

  CHAIN = ' '
  write(chain,3000) ielecc(i)

  call getfbr(chain,nlelt,lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    iel = ifabor(ifac)

    itypfb(ifac) = iparoi

!     - Zone number
    izone = i

!      - Allocation of zone number
    izfppp(ifac) = izone

    if ( ielect(i) .ne. 0 ) then
      icodcl(ifac,isca(ipotr))   = 1
      rcodcl(ifac,isca(ipotr),1) = ur(ielect(i),ielecb(i))

      if ( ippmod(ieljou).eq.4  ) then
        icodcl(ifac,isca(ipoti))   = 1
        rcodcl(ifac,isca(ipoti),1) = ui(ielect(i),ielecb(i))
      endif

    else

      ii = ipotr
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0

      if ( ippmod(ieljou).eq. 4   ) then
        ii = ipoti
        icodcl(ifac,isca(ii))   = 3
        rcodcl(ifac,isca(ii),3) = 0.d0
      endif
    endif

  enddo

enddo

! 3 - Test, if not any reference transformer
!      a piece of wall may be at ground.

if ( ntfref .eq. 0 ) then

  itrouv = 0
  do ifac = 1, nfabor

    if ( itypfb(ifac) .eq. iparoi ) then

      if (icodcl(ifac,isca(ipotr)) .eq. 1 ) then

        if ( ippmod(ieljou).eq.3 ) then

          if ( abs(rcodcl(ifac,isca(ipotr),1)).lt.1.e-20 ) then
            itrouv = 1
          endif

        else if ( ippmod(ieljou).eq.4 ) then

          if (icodcl(ifac,isca(ipoti)) .eq. 1 ) then

            if (abs(rcodcl(ifac,isca(ipotr),1)).lt.1.e-20         &
           .and.abs(rcodcl(ifac,isca(ipoti),1)).lt.1.e-20 ) then
              itrouv = 1
            endif
          endif
        endif
      endif
    endif
  enddo

  if ( itrouv .eq. 0 ) then
    write(nfecra,1000)
    call csexit (1)
  endif

endif

!--------
! Formats
!--------

 1000 format(1X,' ERROR in JOULE : '                          ,/, &
       1X,' ====================   '                          ,/, &
      10X,' Lack of reference : choose a transformer for wich',/, &
      10X,' offset is assumed zero or a face at ground on the',/, &
          ' boundary')
 1500 format(/,2X,' ** INFORMATIONS ON TRANSFOMERS           ',/, &
         2X,'    ---------------------------------------'/    ,/, &
         1X,'      ---------------------------------'         ,/, &
         1X,'      Number of Transfo        UROFF    '        ,/, &
         1X,'      ---------------------------------')
 1501 format(1X,'      ---------------------------------')
 2000 format(10x,i6,12x,e12.5)
 3000 format(i7)

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_user_boundary_conditions

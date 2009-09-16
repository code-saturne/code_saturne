!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine usclim &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Fill boundary conditions arrays (icodcl, rcodcl) for unknown variables.


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.

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

! Similarly, interior faces and cells can be identified using the 'getfac'
! and 'getcel' subroutines (respectively). Their syntax are identical to
! 'getfbr' syntax.

! For a more thorough description of the criteria syntax, it can be referred
! to the user guide.


! Boundary condition types
! ========================

! Boundary conditions may be assigned in two ways.


!    For "standard" boundary conditions:
!    -----------------------------------

!     (inlet, free outlet, wall, symmetry), one defines a code in the 'itypfb'
!     array (of dimensions number of boundary faces, number of phases).
!     This code will then be used by a non-user subroutine to assign the
!     following conditions (scalars in particular will receive the conditions
!     of the phase to which they are assigned). Thus:

!     Code      |  Boundary type
!     --------------------------
!      ientre   |   Inlet
!      isolib   |   Free outlet
!      isymet   |   Symmetry
!      iparoi   |   Wall (smooth)
!      iparug   |   Rough wall

!     These integers are defined elsewhere (in paramx.h header).
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
!                       rcodcl(ifac, iu, 3) (values for iv and iw are not used)
!       -> Specific code and prescribed temperature value at wall if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 6
!                       rcodcl(ifac, ivar, 1) = prescribed temperature
!                       rcodcl(ifac, ivar, 3) = dynamic roughness height
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
!           rcodcl(ifac, ivar, 3) = roughd
!         For other scalars, thermal roughness
!           rcodcl(ifac, ivar, 3) = rought


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
!        -> set itypfb(ifac, iphas) = iparoi which prescribes default wall
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


!       Note how to access some variables (for phase    'iphas'
!                                              variable 'ivar'
!                                              scalar   'iscal'):

! Cell values  (let iel = ifabor(ifac))

! * Density:                                 propce(iel, ipproc(irom(iphas)))
! * Dynamic molecular viscosity:             propce(iel, ipproc(iviscl(iphas)))
! * Turbulent viscosity:                     propce(iel, ipproc(ivisct(iphas)))
! * Specific heat:                           propce(iel, ipproc(icp(iphas))
! * Diffusivity(lambda):                     propce(iel, ipproc(ivisls(iscal)))

! Boundary face values

! * Density:                                 propfb(ifac, ipprob(irom(iphas)))
! * Mass flux (for convecting 'ivar'):       propfb(ifac, ipprob(ifluma(ivar)))

! * For other values: take as an approximation the value in the adjacent cell
!                     i.e. as above with iel = ifabor(ifac).


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! maxelt           !  e ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfac(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb(nfabor    ! ia ! <-- ! indirection for boundary faces ordering)       !
!  (nfabor, nphas) !    !     !                                                !
! itypfb           ! ia ! --> ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! idevel(nideve)   ! ia ! <-- ! integer work array for temporary developpement !
! ituser(nituse    ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfavor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!                  !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! coefu            ! ra ! --- ! work array                                     !
!  (nfabor, 3)     !    !     !  (computation of pressure gradient)            !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary developpement    !
! rtuser(nituse    ! ra ! <-- ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "ihmpre.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml), maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac, iel, ii, ivar, iphas
integer          ilelt, nlelt
double precision uref2, d2s3
double precision rhomoy, dh, ustar2
double precision xintur
double precision xkent, xeent

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in definition of boundary conditions',/,   &
'@    =======',/,                                                 &
'@     The user subroutine ''usclim'' must be completed.',/,      &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initialization
!===============================================================================

idebia = idbia0
idebra = idbra0

d2s3 = 2.d0/3.d0

!===============================================================================
! 2.  Assign boundary conditions to boundary faces here

!     One may use selection criteria to filter boundary case subsets
!       Loop on faces from a subset
!         Set the boundary condition for each face
!===============================================================================

! --- For boundary faces of color 2 and x <= 0.01,
!     assign an inlet for all phases
call getfbr('2 and x < 0.01', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel = ifabor(ifac)

  do iphas = 1, nphas

    itypfb(ifac,iphas) = ientre

    rcodcl(ifac,iu(iphas),1) = 1.1d0
    rcodcl(ifac,iv(iphas),1) = 1.1d0
    rcodcl(ifac,iw(iphas),1) = 1.1d0

    uref2 = rcodcl(ifac,iu(iphas),1)**2  &
           +rcodcl(ifac,iv(iphas),1)**2  &
           +rcodcl(ifac,iw(iphas),1)**2
    uref2 = max(uref2,1.d-12)


    !   Turbulence example computed using equations valid for a pipe.

    !   We will be careful to specify a hydraulic diameter adapted
    !     to the current inlet.

    !   We will also be careful if necessary to use a more precise
    !     formula for the dynamic viscosity use in the calculation of
    !     the Reynolds number (especially if it is variable, it may be
    !     useful to take the law from 'usphyv'. Here, we use by default
    !     the 'viscl0" value given in 'usini1'.
    !   Regarding the density, we have acess to its value at boundary
    !     faces (romb) so this value is the one used here (specifically,
    !     it is consistent with the processing in 'usphyv', in case of
    !     variable density)

    !     Hydraulic diameter
    dh     = 0.075d0

    !   Calculation of friction velocity squared (ustar2)
    !     and of k and epsilon at the inlet (xkent and xeent) using
    !     standard laws for a circular pipe
    !     (their initialization is not needed here but is good practice).
    rhomoy = propfb(ifac,ipprob(irom(iphas)))
    ustar2 = 0.d0
    xkent  = epzero
    xeent  = epzero

    call keendb                                            &
    !==========
        ( uref2, dh, rhomoy, viscl0(iphas), cmu, xkappa,   &
          ustar2, xkent, xeent )

    ! itytur is a flag equal to iturb/10
    if    (itytur(iphas).eq.2) then

      rcodcl(ifac,ik(iphas),1)  = xkent
      rcodcl(ifac,iep(iphas),1) = xeent

    elseif(itytur(iphas).eq.3) then

      rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir12(iphas),1) = 0.d0
      rcodcl(ifac,ir13(iphas),1) = 0.d0
      rcodcl(ifac,ir23(iphas),1) = 0.d0
      rcodcl(ifac,iep(iphas),1)  = xeent

    elseif(iturb(iphas).eq.50) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iep(iphas),1)  = xeent
      rcodcl(ifac,iphi(iphas),1) = d2s3
      rcodcl(ifac,ifb(iphas),1)  = 0.d0

    elseif(iturb(iphas).eq.60) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

    endif

    ! --- Handle scalars attached to the current phase
    if(nscal.gt.0) then
      do ii = 1, nscal
        if(iphsca(ii).eq.iphas) then
          rcodcl(ifac,isca(ii),1) = 1.d0
        endif
      enddo
    endif

  enddo

enddo

! --- Prescribe at boundary faces of color 3 an inlet for all phases
call getfbr('3', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel  = ifabor(ifac)

  do iphas = 1, nphas

    itypfb(ifac,iphas) = ientre

    rcodcl(ifac,iu(iphas),1) = 1.1d0
    rcodcl(ifac,iv(iphas),1) = 1.1d0
    rcodcl(ifac,iw(iphas),1) = 1.1d0

    uref2 = rcodcl(ifac,iu(iphas),1)**2   &
           +rcodcl(ifac,iv(iphas),1)**2   &
           +rcodcl(ifac,iw(iphas),1)**2
    uref2 = max(uref2,1.d-12)

    ! Turbulence example computed using turbulence intensity data.

    ! We will be careful to specify a hydraulic diameter adapted
    !   to the current inlet.

    ! Hydraulic diameter

    dh     = 0.075d0
    ! Turbulence intensity
    xintur = 0.02d0

    ! Calculation of k and epsilon at the inlet (xkent and xeent) using
    !   the turbulence intensity and standard laws for a circular pipe
    !   (their initialization is not needed here but is good practice)
    xkent  = epzero
    xeent  = epzero

    call keenin                                                   &
    !==========
        ( uref2, xintur, dh, cmu, xkappa, xkent, xeent )

    ! itytur is a flag equal to iturb/10
    if    (itytur(iphas).eq.2) then

      rcodcl(ifac,ik(iphas),1)  = xkent
      rcodcl(ifac,iep(iphas),1) = xeent

    elseif(itytur(iphas).eq.3) then

      rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir12(iphas),1) = 0.d0
      rcodcl(ifac,ir13(iphas),1) = 0.d0
      rcodcl(ifac,ir23(iphas),1) = 0.d0
      rcodcl(ifac,iep(iphas),1)  = xeent

    elseif(iturb(iphas).eq.50) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iep(iphas),1)  = xeent
      rcodcl(ifac,iphi(iphas),1) = d2s3
      rcodcl(ifac,ifb(iphas),1)  = 0.d0

    elseif(iturb(iphas).eq.60) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

    endif

    ! --- Handle scalars attached to the current phase
    if(nscal.gt.0) then
      do ii = 1, nscal
        if(iphsca(ii).eq.iphas) then
          rcodcl(ifac,isca(ii),1) = 1.d0
        endif
      enddo
    endif

  enddo
enddo

! --- Prescribe at boundary faces of group 'outlet' an outlet for all phases
call getfbr('outlet', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Outlet: zero flux for velocity and temperature, prescribed pressure
  !         Note that the pressure will be set to P0 at the first
  !         free outlet face (isolib)

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = isolib
  enddo

enddo

! --- Prescribe at boundary faces of color 5 a wall for all phases
call getfbr('5', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = iparoi
  enddo

  ! If sliding wall with velocity u(1) = 1:
  ! rcodcl(ifac, iu(1), 1) = 1.d0

  ! If sliding wall with velocity u = 0: nothing to do

  if(nscal.gt.0) then

    ! If temperature prescribed to 20 with wall law (scalar ii=1):
    ! ii = 1
    ! icodcl(ifac, isca(ii))   = 5
    ! rcodcl(ifac, isca(ii), 1) = 20.d0

    ! If temperature prescribed to 50 with no wall law (simple Dirichlet)
    !   with exchange coefficient 8 (scalar ii=2):
    ! ii = 2
    ! icodcl(ifac, isca(ii))    = 1
    ! rcodcl(ifac, isca(ii),1)  = 50.d0
    ! rcodcl(ifac, isca(ii), 2) = 8.d0

    ! If flux prescribed to 4.d0 (scalar ii=3):
    ! ii = 3
    ! icodcl(ifac, isca(ii))    = 3
    ! rcodcl(ifac, isca(ii), 3) = 4.D0

  endif
enddo

! --- Prescribe at boundary faces of color 7 a rough wall for all phases
call getfbr('7', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       rough friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = iparug

  ! Roughness for velocity: 1cm
    rcodcl(ifac,iu(iphas),3) = 0.01d0
  enddo

  ! If sliding wall with velocity u(1) = 1:
  ! rcodcl(ifac, iu(1), 1) = 1.d0

  ! If sliding wall with velocity u = 0: nothing to do
  if(nscal.gt.0) then

    ! If temperature prescribed to 20 with rough wall law (scalar ii=1)
    !   with roughness of 1 cm:
    ! ii = 1
    ! icodcl(ifac, isca(ii))   = 5
    ! rcodcl(ifac, isca(ii), 1) = 20.d0
    ! rcodcl(ifac, isca(ii), 3) = 0.01.d0

    ! If flux prescribed to 4.d0 (scalar ii=3):
    ! ii = 3
    ! icodcl(ifac, isca(ii))    = 3
    ! rcodcl(ifac, isca(ii), 3) = 4.D0

  endif

enddo

! --- Prescribe at boundary faces of color 4 a symmetry
call getfbr('4', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Symmetries

  do iphas = 1, nphas
    itypfb(ifac,iphas)   = isymet
  enddo

enddo

! From this point on, the examples given are destined to an ADVANCED USE,
!   the user being led to define complex boundary condition combinations.


! --- Example of specific boundary conditions fully defined by the user,
!     on the basis of wall conditions.
!     selection (mass flow computation, specific logging, ...)
!     We prescribe for color 1234 a wall for all phases, with in addition:
!       - a Dirichlet condition on velocity of phase 1 (sliding wall
!         with no-slip condition)
!       - a Dirichlet condition on the first scalar.

call getfbr('1234', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  do iphas = 1, nphas
    itypfb(ifac,iphas) = iparoi
  enddo

  iphas = 1
  icodcl(ifac,iu(iphas) )  = 1
  rcodcl(ifac,iu(iphas),1) = 1.d0
  rcodcl(ifac,iu(iphas),2) = rinfin
  rcodcl(ifac,iu(iphas),3) = 0.d0
  icodcl(ifac,iv(iphas) )  = 1
  rcodcl(ifac,iv(iphas),1) = 0.d0
  rcodcl(ifac,iv(iphas),2) = rinfin
  rcodcl(ifac,iv(iphas),3) = 0.d0
  icodcl(ifac,iw(iphas) )  = 1
  rcodcl(ifac,iw(iphas),1) = 0.d0
  rcodcl(ifac,iw(iphas),2) = rinfin
  rcodcl(ifac,iw(iphas),3) = 0.d0

  ivar = isca(1)
  icodcl(ifac,ivar )  = 1
  rcodcl(ifac,ivar,1) = 10.d0
  rcodcl(ifac,ivar,2) = rinfin
  rcodcl(ifac,ivar,3) = 0.d0

enddo

! --- Example of specific boundary conditions fully defined by the user,
!     with no definition of a specific type.
!     We prescribe at color 5678 a homogeneous Neumann condition for
!     all variables (whatever the phase).

call getfbr('5678', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! CAUTION: the value of itypfb must be assigned to iindef

  do iphas = 1, nphas
    itypfb(ifac,iphas) = iindef
  enddo

  do ii = 1, nvar
    icodcl(ifac,ii )  = 3
    rcodcl(ifac,ii,1) = 0.d0
    rcodcl(ifac,ii,2) = rinfin
    rcodcl(ifac,ii,3) = 0.d0
  enddo

enddo

! --- Example of specific boundary conditions fully defined by the user,
!     with the definition of a specific type, for example for future
!     selection (mass flow computation, specific logging, ...)
!     We prescribe for color 6789 a homogeneous Neumann condition for
!     all variables (whatever the phase), except for the first
!     scalar, for which we select a homogeneous Dirichlet.

call getfbr('6789', nlelt, lstelt)
!==========
do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! CAUTION: the value of itypfb must be different from
!          iparoi, ientre, isymet, isolib, iindef,
!          greater than or equal to 1, and
!          less than or equal to ntypmx;
!          these integers are defined in paramx.h

  do iphas = 1, nphas
    itypfb(ifac,iphas) = 89
  enddo

  do ii = 1, nvar
    icodcl(ifac,ii )  = 3
    rcodcl(ifac,ii,1) = 0.d0
    rcodcl(ifac,ii,2) = rinfin
    rcodcl(ifac,ii,3) = 0.d0
  enddo

  icodcl(ifac,isca(1) )  = 1
  rcodcl(ifac,isca(1),1) = 0.d0
  rcodcl(ifac,isca(1),2) = rinfin
  rcodcl(ifac,isca(1),3) = 0.d0

enddo

!----
! Formats
!----

!----
! End
!----

return
end subroutine

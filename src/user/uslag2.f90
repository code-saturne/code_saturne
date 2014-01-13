!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine uslag2 &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   dt     , rtpa   , propce ,                                     &
   ettp   , tepa   )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   -------------------------------------

!    User subroutine (Mandatory intervention)

!    User subroutine for the boundary conditions associated to the particles
!    (inlet and treatment of the other boundaries)


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! itrifb(nfabor)   ! ia ! <-- ! indirection for the sorting of the             !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! ifrlag(nfabor    ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous          !
! (ncelet,*)       !    !     ! time step                                      !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use cstnum
use cstphy
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)

! Local variables

integer          ifac , izone, nbclas, iclas
integer          icoal , ilayer
integer          ilelt, nlelt

double precision pis6 , mp0 , temp

integer, dimension(ndlaim) :: iczpar
double precision, dimension(ndlagm) :: rczpar

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
  !==========
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP AT THE ENTRANCE OF THE BOUNDARY           ',/,&
'@    ========                                                ',/,&
'@     CONDITIONS OF THE LAGRANGIAN MODULE:                   ',/,&
'@     THE USER SUBROUTINE uslag2 MUST BE FILLED              ',/,&
'@                                                            ',/,&
'@  The calculation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  Memory management
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 2. Initialization
!===============================================================================

pis6 = pi / 6.d0

!===============================================================================
! 3. Construction of the boundary zones
!===============================================================================


!     Definition of the boundary zones
!     --------------------------------

!     For the Lagrangian module, the user defines nfrlag boundary zones
!     from the color of the boundary faces, or more generally from their
!     properties (colors, groups..) or from the boundary conditions prescribed
!     in cs_user_boundary_conditions, or even from their coordinates. To do
!     that, we fill the ifrlag(nfabor) array which gives for every boundary
!     face the number of the zone to which it belongs ifrlag(ifac)
!
!     Be careful, all the boundary faces must have been assigned.
!
!     The number of the zones (thus the values of ifrlag(ifac)) is arbitrarily
!     chosen by the user, but must be a positive integer and inferior or equal
!     to nflagm (parameter prescribed in lagpar.h).
!
!     Afterwards, we assign to every zone a type named itylag that will be used
!     to impose global boundary conditions.

izone = -1

! ---> First zone numbered izone=1 ( = color 10)
call getfbr('10',nlelt,lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 1
  ifrlag(ifac) = izone

enddo

! ---> Second zone numbered izone=2 ( = part of color 4)
call getfbr('4 and y < 1.0',nlelt,lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 2
  ifrlag(ifac) = izone

enddo

! ---> Third zone numbered izone=3 ( = inlet)
do ifac = 1, nfabor
  if(itypfb(ifac).eq.ientre) then
    izone        = 4
    ifrlag(ifac) = izone
  endif
enddo

! ---> Nth zone numbered izone=5 (= color 3)
call getfbr('3',nlelt,lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 5
  ifrlag(ifac) = izone

enddo


!===============================================================================
! 4. Injection per particle class into the calculation domain
!===============================================================================

!   TO PROVIDE INFORMATION ABOUT THE PARTICLE CLASSES,
!   WE FOLLOW A TWO-STEP PROCEDURE:

!   1) FIRST, THE NUMBER OF PARTICLE CLASSES IS PRESCRIBED
!      FOR EACH BOUNDARY ZONE: IUSNCL (by default, this parameter is equal to zero)

!   2) AFTERWARDS, FOR EACH ZONE AND FOR EACH CLASS, WE PRESCRIBE
!      THE PHYSICAL PROPERTIES OF THE PARTICLES
!




! --> Number of particle classes entering the domain
!   We assign here the number of classes for each zone previously identified.
!
!   This number is zero by default.
!   The maximal number of classes is nclagm (defined in lagpar.h)

! ---> First zone numbered izone = 1: 1 class injected
izone     = 1
nbclas    = 1
iusncl(izone) = nbclas

! ---> Second zone numbered izone = 2: 0 class injected
izone     = 2
nbclas    = 0
iusncl(izone) = nbclas

! ---> Third zone numbered izone = 4 : 0 class injected
izone     = 4
nbclas    = 0
iusncl(izone) = nbclas

! ---> Zone numbered izone = 5 : 0 class injected
izone     = 5
nbclas    = 0
iusncl(izone) = nbclas


! --> For every class associated with a zone,
!     we give the following information.


!     iusncl number of classes per zone
!     iusclb boundary conditions for the particles
!     = ientrl -> zone of particle inlet
!     = isortl -> particle outlet
!     = irebol -> rebound of the particles
!     = idepo1 -> definitive deposition
!     = idepo2 -> definitive deposition, but the particle remains in memory
!                 (useful only if iensi2 = 1)
!     = idepfa -> deposition of the particle with DLVO forces
!     = iencrl -> fouling (coal only iphyla = 2)
!     = isymtl -> symmetry condition for the particles (zero flux)

!     Array iczpar:
!     ============
!        ijnbp : number of particles per class and per zone
!        ijfre : injection frequency. If ijfre = 0, then the injection
!                occurs only at the first absolute iteration.
!        iclst : number of the group to which the particle belongs
!                (only if one wishes to calculate statistics per group)
!        ijuvw : type of condition on the velocity
!                  = -1 imposed flow velocity
!                  =  0 imposed velocity along the normal direction of the
!                      boundary face, with norm equal to rczpar(iuno)
!                  =  1 imposed velocity: we prescribe   rczpar(iupt)
!                                                        rczpar(ivpt)
!                                                        rczpar(iwpt)
!                  =  2 user-defined profile
!        ijprtp : type of temperature condition
!                  =  1 imposed temperature: we prescribe rczpar(itpt)
!                  =  2 user-defined profile
!        ijprdp : type of diameter condition
!                  =  1 imposed diameter: we prescribe  rczpar(idpt)
!                                                       rczpar(ivdpt)
!                  =  2 user-defined profile
!        ijprtp : type of temperature condition
!                  =  1 imposed temperature: we prescribe rczpar(itpt)
!                  =  2 user-defined profile
!        inuchl : number of the coal of the particle (only if iphyla = 2)
!        irawcl : type of coal injection composition (only if iphyla = 2)
!                  = 0 coal injected with an user-defined composition
!                  = 1 raw coal injection

!     Array rczpar:
!     ============
!        iuno  : Norm of the velocity (m/s)
!        iupt  : Velocity along the X axis, for each class and for each zone (m/s)
!        ivpt  : Velocity along the Y axis, for each class and for each zone (m/s)
!        iwpt  : Velocity along the Z axis, for each class and for each zone (m/s
!        ipoit : Statistical weight (number of samples) associated
!                to the particle (automatically computed to respect a mass
!                flow rate if it is defined)
!        idebt : Mass flow rate (kg/s)

!        Physical characteristics:
!          idpt   : diameter (m)
!          ivdpt  : standard deviation of the diameter (m)
!          iropt  : density (kg/m3)
!          itpt   : temperature in Celsius degress (no enthalpy)
!          icpt   : specific heat (J/kg/K)
!          iepsi  : emissivity (if =0 then no radiative effect is taken into account)

!         If coal (iphyla=2)
!            ihpt    : temperature in Kelvin degres (no enthalpy) (layer by layer)
!            ifrmwt  : mass fraction of moisture in the particle
!            ifrmch  : mass fraction of reactive coal in the particle (layer by layer)
!            ifrmck  : mass fraction of coke coal in the particle (layer by layer)
!            irdck   : coke diameter (m)
!            ird0p   : coal particle initial diameter (m)
!            irhock0 : density of coke at the end of the pyrolysis (kg/m³) (layer by layer)


! ---> EXAMPLE : First zone, numbered IZONE = 1 (NBCLAS classes)
!        IUSCLB : adherence of the particle to a boundary face
!        IJNBP  : 10 particles for each class,
!        IJFRE  : injection every other time step
!        IJUVW, IUPT, IVPT, IWPT : imposed velocity on 1.1D0, 0.0D0, 0.0D0
!        ICPT   : cp equal to 10000
!        ITPT   : temperature equal to 25 Celsius degress
!        IDPT   : diameter equal to 50.E-6 m
!        IEPSI  : emissivity equal to 0.7
!        IVDPT  : constant diameter ==> standard deviation null
!        IROPT  : density
!        IPOIT  : statistical weight (number of physical particles
!                 represented by one statistical particle)
!        IDEBT  : mass flow rate


izone         = 1
nbclas        = iusncl(izone)
iusclb(izone) =  ientrl

do iclas  = 1, nbclas

  ! Ensure defaults are set

  call lagr_init_zone_class_param(iczpar, rczpar)
  !==============================

  ! Now define parameters for this class and zone

  iczpar(ijnbp) = 10
  iczpar(ijfre) = 2

  if (nbclst.gt.0) then
    iczpar(iclst) = 1
  endif

  iczpar(ijuvw) = -1
  rczpar(iupt)  = 1.1d0
  rczpar(ivpt)  = 0.0d0
  rczpar(iwpt)  = 0.0d0
  iczpar(ijprpd)= 1
  rczpar(ipoit) = 1.d0
  rczpar(idebt) = 0.d0

  ! if the physics is " simple"

  if (iphyla.eq.0 .or. iphyla.eq.1) then

    ! Mean value and standard deviation of the diameter

    iczpar(ijprdp)= 1
    rczpar(idpt)  = 50.d-6
    rczpar(ivdpt) = 0.d0

    ! Density

    rczpar(iropt) = 2500.d0

    if (iphyla.eq.1) then

      ! Temperature and Cp

      if ( itpvar.eq.1 ) then
        iczpar(ijprtp) = 1
        rczpar(itpt)    = 20.d0

        rczpar(icpt)    = 1400.d0
        rczpar(iepsi)   = 0.7d0
      endif

    endif

    ! Coal

  else if ( iphyla.eq.2 ) then

    ! CAUTION :   1) To transport and burn coal particles with the Lagrangian
    !                module, a specific physics for the dispersed phase must
    !                be activated for the carrier phase.
    !
    !             2) The physical properties of the coal particles are known
    !                from the thermo-chemical file: dp_FCP
    !
    !             3) For the current phase ICLAS, and for the current boundary
    !                zone NB, we assign to the coal particles the properties of
    !                the coal ICOAL of the ICOAL class taken from the file dp_FCP.
    !
    !             4) icoal : number of the coal between 1 and ncharb defined by
    !                the user in the file dp_FCP.


    ! Mean value and standard deviation of the diameter

    rczpar(idpt)  = 1.0d-5
    rczpar(ivdpt) = 0.d0

    ! Number of the coal
    icoal = 1
    iczpar(inuchl) = icoal

    ! Temperature (in K)
    temp = 800.d0
    do ilayer = 1, nlayer
      rczpar(ihpt(ilayer)) = temp
    enddo

    ! Raw coal or user defined coal injection condition
    iczpar(irawcl) = 1

    if (iczpar(irawcl) .eq.  0) then
      ! Example of user-defined injection, coal after devolatilisation

      ! Density
      rczpar(iropt) = xashch(icoal)*rho0ch(icoal) +                          &
                      (1.0d0-xwatch(icoal)-xashch(icoal)) * rho0ch(icoal)    &
                       * (1.d0-(y1ch(icoal)+y2ch(icoal))/2.0d0)

      ! Specific heat
      rczpar(icpt) = cp2ch(icoal)

      ! water mass fraction in the particle
      rczpar(ifrmwt) = 0.0d0
      ! reactive coal mass fraction in the particle
      do ilayer = 1, nlayer
        rczpar(ifrmch(ilayer)) = 0.0d0
      enddo
      ! coke mass fraction in the particle
      do ilayer = 1, nlayer
        rczpar(ifrmck(ilayer)) = ((1.0d0-xwatch(icoal)-xashch(icoal))        &
                                   * rho0ch(icoal)                           &
                                   * (1.d0-(y1ch(icoal)+y2ch(icoal))/2.0d0)) &
                                  / rczpar(iropt)
      enddo

      ! coke diameter
      rczpar(irdck)   = rczpar(idpt)
      ! initial particle diameter
      rczpar(ird0p)   = rczpar(idpt)
      ! coke density after pyrolysis
      do ilayer = 1, nlayer
        rczpar(irhock0(ilayer)) = rczpar(iropt)*rczpar(ifrmck(ilayer))
      enddo
    endif

  endif

  ! Complete definition of parameters for this class and zone

  call lagr_define_zone_class_param(iclas, izone, iczpar, rczpar)
  !================================

enddo

! ---> Second zone, numbered izone = 2
!        IUSCLB : rebound of the particle

izone     = 2
iusclb (izone)         =  irebol


! same procedure for the other zones...

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine uslag2

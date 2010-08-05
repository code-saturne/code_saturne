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

subroutine uslag2 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   ,                                              &
   rdevel , rtuser , ra     )

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
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


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
!                  !    !     !                                                !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! nideve nrdeve    ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse nrtuse    ! i  ! <-- ! sizes of ituser and rtuser arrays              !
!                  !    !     !                                                !
! ifacel           ! ia ! <-- ! interior faces -> cells connectivity           !
! (2, nfac)        !    !     !                                                !
! ifabor           ! ia ! <-- ! boundary faces -> cells connectivity           !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! ia ! <-- ! boundary face family numbers                   !
! (nfabor)         !    !     !                                                !
! ifmcel           ! ia ! <-- ! cell family numbers                            !
! (ncelet)         !    !     !                                                !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml,nprfml    !    !     !                                                !
! maxelt           !  i ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ipnfac           ! ia ! <-- ! interior faces -> vertices index (optional)    !
!   (lndfac)       !    !     !                                                !
! nodfac           ! ia ! <-- ! interior faces -> vertices list (optional)     !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! ia ! <-- ! boundary faces -> vertices index (optional)    !
!   (lndfbr)       !    !     !                                                !
! nodfbr           ! ia ! <-- ! boundary faces -> vertices list  (optional)    !
!   (nfabor+1)     !    !     !                                                !
! itrifb(nfabor    ! ia ! <-- ! indirection for the sorting of the             !
!  nphas      )    !    !     ! boundary faces                                 !
! itypfb(nfabor    ! ia ! <-- ! type of the boundary faces                     !
!  nphas      )    !    !     !                                                !
! ifrlag(nfabor    ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! --- ! complementary user array of integers           !
! ia(*)            ! ia ! <-- ! macro array of integers                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
! (ndim,nnod)      !    !     !                                                !
! volume           ! ra ! <-- ! cell volumes                                   !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous          !
! (ncelet,*)       !    !     ! time step                                      !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! rdevel(nrdeve    ! ra ! <-- ! dev. complementary array of reals              !
! rtuser(nrtuse    ! ra ! <-- ! user complementary array of reals              !
! ra(*)            ! ra ! --- ! macro array of reals                           !
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

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifac , izone, nbclas, iclas
integer          icha
integer          ilelt, nlelt

double precision pis6 , mp0 , temp

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

idebia = idbia0
idebra = idbra0

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
!     from the color of the boundary faces, or more generally from their properties
!     (colors, groups..) or from the boundary conditions prescribed in usclim, or
!     even from their coordinates. To do that, we fill the ifrlag(nfabor) array
!     which gives for every boundary face the number of the zone to which it belongs ifrlag(ifac)
!
!     Be careful, all the boundary faces must have been affected.
!
!     The number of the zones (thus the values of ifrlag(ifac)) is arbitrarily
!     chosen by the user, but must be a positive integer and inferior or equal to
!     nflagm (parameter prescribed in lagpar.h).
!
!     Afterwards, we assign to every zone a type named itylag that will be used
!     to impose global boundary conditions.

izone = -1

! ---> First zone numbered izone=1 ( = color 10)
CALL GETFBR('10',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 1
  ifrlag(ifac) = izone

enddo

! ---> Second zone numbered izone=2 ( = part of color 4)
CALL GETFBR('4 and Y < 1.0',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 2
  ifrlag(ifac) = izone

enddo

! ---> Third zone numbered izone=3 ( = inlet phase 1)
do ifac = 1, nfabor
  if(itypfb(ifac,1).eq.ientre) then
    izone        = 4
    ifrlag(ifac) = izone
  endif
enddo

! ---> Nth zone numbered izone=5 (= color 3)
CALL GETFBR('3',NLELT,LSTELT)
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
!     we give the followong information.


!     iusncl number of classes per zone
!     iusclb boundary conditions for the particles
!     = ientrl -> zone of particle inlet
!     = isortl -> particle outlet
!     = irebol -> rebound of the particles
!     = idepo1 -> definitive deposition
!     = idepo2 -> definitive deposition, but the particle remains in memory
!                 (useful only if iensi2 = 1)
!     = idepo3 -> deposition and resuspension possible
!                 following the conditions of the flow
!     = idepfa -> deposition of the particle with attachment force,
!                 the velocity is conserved, and resuspension is possible
!                 (Possible if ladlvo = 1 )
!     = iencrl -> fouling (coal only iphyla = 2)
!     = jbord1 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord2 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord3 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord4 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord5 -> user-defined particle/boundary interaction (cf. uslabo)



!     Array iuslag :
!     ================
!        ijnbp : number of particles per class and per zone
!        ijfre : injection frequency. If ijfre = 0, then the injection
!                occurs only at the first absolute iteration.
!        iclst : number of the group to which the particle belongs
!                (only if one wishes to calculate statistics per group)
!        ijuvw : type of condition on the velocity
!                  = -1 imposed flow velocity
!                  =  0 imposed velocity along the normal direction of the
!                      boundary face, with norm equal to RUSLAG(ICLAS,IZONE,IUNO)
!                  =  1 imposed velocity: we prescribe   RUSLAG(ICLAS,IZONE,IUPT)
!                                                        RUSLAG(ICLAS,IZONE,IVPT)
!                                                        RUSLAG(ICLAS,IZONE,IWPT)
!                  =  2 user-defined profile
!        ijprtp : type of temperature condition
!                  =  1 imposed temperature: we prescribe RUSLAG(ICLAS,IZONE,ITPT)
!                  =  2 user-defined profile
!        ijprdp : type of diameter condition
!                  =  1 imposed diameter: we prescribe  RUSLAG(ICLAS,IZONE,IDPT)
!                                                       RUSLAG(ICLAS,IZONE,IVDPT)
!                  =  2 user-defined profile
!        inuchl : number of the coal of the particle (only if iphyla = 2)

!     Array ruslag :
!     ===============
!        iuno  : Norm of the velocity (m/s)
!        iupt  : Velocity along the X axis, for each class and for each zone (m/s)
!        ivpt  : Velocity along the Y axis, for each class and for each zone (m/s)
!        iwpt  : Velocity along the Z axis, for each class and for each zone (m/s)
!        idebt : Mass flow rate (kg/s)
!        ipoit : Statistical weight (number of samples) associated
!                to the particle (automatically computed to respect a mass
!                flow rate if it is defined)

!        Physical characteristics:
!          idpt   : diameter (m)
!          ivdpt  : standard deviation of the diameter (m)
!          itpt   : temperature in Celsius degress (no enthalpy)
!          icpt   : specific heat (J/kg/K)
!          iepsi  : emissivity (if =0 then no radiative effect is taken into account)
!          iropt  : density (kg/m3)

!         If coal (iphyla=2)
!            ihpt  : temperature in Celsius degress (no enthalpy)
!            imcht : mass of reactive coal (kg)
!            imckt : masse of coke (kg)


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


izone     = 1
nbclas    = iusncl(izone)
iusclb (izone)         =  ientrl
do iclas  = 1, nbclas

  iuslag (iclas,izone,ijnbp) = 10
  iuslag (iclas,izone,ijfre) = 2

  if (nbclst.gt.0) then
    iuslag(iclas,izone,iclst) = 1
  endif

  iuslag (iclas,izone,ijuvw) = -1
  ruslag (iclas,izone,iupt)  = 1.1d0
  ruslag (iclas,izone,ivpt)  = 0.0d0
  ruslag (iclas,izone,iwpt)  = 0.0d0
  iuslag (iclas,izone,ijprpd)= 1
  ruslag (iclas,izone,ipoit) = 1.d0
  ruslag (iclas,izone,idebt) = 0.d0

!    if the physics is " simple"

  if ( iphyla.eq.0 .or. iphyla.eq.1 ) then

!        Mean value and standard deviation of the diameter

    iuslag (iclas,izone,ijprdp)= 1
    ruslag (iclas,izone,idpt)  = 50.d-6
    ruslag (iclas,izone,ivdpt) = 0.d0

!        Density

    ruslag(iclas,izone,iropt) = 2500.d0

    if ( iphyla.eq.1 ) then

!        Temperature and Cp

      if ( itpvar.eq.1 ) then
        iuslag (iclas,izone,ijprtp) = 1
        ruslag(iclas,izone,itpt)    = 20.d0

        ruslag(iclas,izone,icpt)    = 1400.d0
        ruslag(iclas,izone,iepsi)   = 0.7d0
      endif

    endif

!    Coal

  else if ( iphyla.eq.2 ) then

!    CAUTION :   1) To transport and burn coal particles with the Lagrangian
!                   module, a specific physics for the dispersed phase must
!                   be activated for the carrier phase.
!
!                2) The physical properties of the coal particles are known
!                   from the thermo-chemical file: dp_FCP
!
!                3) For the current phase ICLAS, and for the current boundary zone
!                   NB, we assign to the coal particles the properties of the coal ICHA
!                   of the icha class taken from the file dp_FCP.
!
!                4) icha : number of the coal between 1 and ncharb defined by the user
!                   in the file dp_FCP.
!


    icha = ichcor(iclas)
    temp = 800.d0

!        Number of the coal

    iuslag(iclas,izone,inuchl) = icha

!        Temperature and Cp

    ruslag(iclas,izone,ihpt) = temp
    ruslag(iclas,izone,icpt) = cp2ch(icha)

!        Mean value and standard deviation of the diameter

    ruslag (iclas,izone,idpt)  = diam20(iclas)
    ruslag (iclas,izone,ivdpt) = 0.d0

!        Density

    ruslag(iclas,izone,iropt) =  rho0ch(icha)

!        Mass of reactive coal and
!        mass of coke (null if the coal has never burnt)

    mp0 = pis6 * ( ruslag(iclas,izone,idpt)**3 )                  &
               * ruslag(iclas,izone,iropt)
    ruslag(iclas,izone,imcht) = mp0 * (1.d0-xashch(icha))
    ruslag(iclas,izone,imckt) = 0.d0

  endif

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

return

end subroutine

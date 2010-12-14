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

subroutine uslaed &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  , ibord  ,                                              &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   auxl1  , auxl2  , auxl3  , w1     , w2     , w3     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose :
! ----------

!   Subroutine of the Lagrangian particle-tracking module :
!   -------------------------------------

!     User subroutine (non-mandatory intervention)

!     Integration of the sde for the user-defined variables.
!     The variables are constant by default.


!                                         d T       T - PIP
!     The sde must be of the form:       ----- = - ---------
!                                         d t         Tca


!     T : IIIIeme user-defined variable, given for the ip particle by
!            T = ETTP(IP,JVLS(IIII))
!            T = ETTPA(IP,JVLS(IIII))

!     Tca : Characteristic time for the sde
!           to be prescribed in the array auxl1

!     PIP : Coefficient of the sde (pseudo right member)
!           to be prescribed in the array auxl2
!
!           If the chosen scheme is first order (nordre=1)
!           then, at the first and only passage pip is expressed
!           as a function of the quantities of the previous time step contained in ettpa
!
!           If the chosen scheme is second order (nordre=2)
!           then, at the first passage (nor=1) pip is expressed as
!           a function of the quantities of the previous time step contained in ettpa,
!           and at the second passage (nor=2) pip is expressed as
!           a function of the quantities of the current time step contained in ettp

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! ibord            ! ia ! <-- ! number of the boundary face of part./wall      !
!   (nbpmax)       !    !     ! interaction                                    !
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! <-- ! complementary user array of integers           !
! ia(*)            ! ia ! --- ! macro array of integers                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! transported variables at the current           !
! (ncelet,*)       !    !     ! and previous time step                         !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! tempct           ! ra ! <-- ! characteristic thermal time and                !
!  (nbpmax,2)      !    !     ! implicit source term of return coupling        !
! tsvar            ! ra ! <-- ! prediction 1st substep for the ivar variable,  !
! (nbpmax,nvp1)    !    !     ! used for the correction at the 2nd substep     !
!                  !    !     !                                                !
! auxl1(nbpmax)    ! ra ! ---   work array                                     !
! auxl2(nbpmax)    ! ra ! --- ! work array                                     !
! auxl3(nbpmax)    ! ra ! --- ! work array                                     !
! w1..w3(ncelet    ! ra ! --- ! work arrays                                    !
! rdevel(nrdeve    ! ra ! <-- ! dev complementary array of reals               !
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
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse

integer          itepa(nbpmax,nivep)  , ibord(nbpmax)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax), auxl2(nbpmax), auxl3(nbpmax)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          npt , iel , iiii , ipl , iphas

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

! We enter this subroutine only if additional variables have been defined in
! uslag1; we must then necessarily define how they are solved.


if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP IN THE LAGRANGIAN MODULE                  ',/,&
'@    =========                                               ',/,&
'@     THE USER SUBROUTINE uslaed MUST BE FILLED              ',/,&
'@                                                            ',/,&
'@  The calculation will not be run                           ',/,&
'@                                                            ',/,&
'@  Additional variables have been declared in                ',/,&
'@    uslag1 (NVLS=)                                          ',/,&
'@  The subroutine uslaed must be filled to precise           ',/, &
'@    the stochastic differential equation to be solved       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! FIXME : TODO write a user example


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initializations
!===============================================================================

idebia = idbia0
idebra = idbra0

iphas = ilphas

!===============================================================================
! 2. Characteristic time of the current sde
!===============================================================================

! Loop on the additional variables

do iiii = 1,nvls

!      Number of the treated variable in ettp

  ipl = jvls(iiii)

  do npt = 1,nbpart

    if ( itepa(npt,jisor).gt.0 ) then

      iel = itepa(npt,jisor)

!     Characteristic time tca of the differential equation
!     This example must be adapted to the case

      auxl1(npt) = 1.d0

!     Prediction at the first substep
!     This example must be adapted to the case

      if (nor.eq.1) then
        auxl2(npt) = ettpa(npt,ipl)
      else

!     Correction at the second substep
!     This example must be adapted to the case

        auxl2(npt) = ettp(npt,ipl)
      endif

    endif
  enddo

!===============================================================================
! 3. Integration of the variable ipl
!===============================================================================

  call lagitg                                                     &
  !==========
   ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ipl    ,                                                     &
     itepa(1,jisor)  , ibord  ,                                   &
     ettp   , ettpa  , auxl1  , auxl2  , tsvar  )

enddo

!===============================================================================

!----
! End
!----

end subroutine

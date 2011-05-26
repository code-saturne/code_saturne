!-------------------------------------------------------------------------------

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

subroutine cfdttv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iwarnp ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   wcf    ,                                                       &
   wflmas , wflmab , viscb  , w1     , w2     , w3     ,          &
   w4     , w5     , w6     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA CONTRAINTE LIEE AU CFL POUR L'ALGO COMPRESSIBLE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,nvar)    !    !     !  source de masse                               !
!                  !    !     ! pour ivar=ipr, smacel=flux de masse            !
! wcf(ncelet)      ! tr ! --> ! contrainte compressible                        !
! wflmas(nfac)     ! tr ! --- ! tab de trav aux faces internes                 !
! wflmab(nfabor    ! tr ! --- ! tab de trav aux faces de bord                  !
! viscb(nfabor     ! tr ! --- ! tab de trav aux faces de bord                  !
! w1..6 (ncelet    ! tr ! --- ! tableaux de travail                            !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstnum
use cstphy
use optcal
use entsor
use parall
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iwarnp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision wcf(ncelet)
double precision wflmas(nfac), wflmab(nfabor), viscb(nfabor)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia, ifinra
integer          ifac  , iel   , ivar  , iscal
integer          init

double precision, allocatable, dimension(:) :: viscf
double precision, allocatable, dimension(:,:) :: coefu
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10, w11, w12

!===============================================================================
!===============================================================================
! 0.  INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(viscf(nfac))
allocate(coefu(nfabor,3))

! Allocate work arrays
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(w10(ncelet), w11(ncelet), w12(ncelet))

idebia = idbia0
idebra = idbra0

iscal  = irho
ivar   = isca(iscal)

!===============================================================================
! 1. CALCUL DE LA CONDITION CFL ASSOCIEE A LA MASSE VOLUMIQUE
!===============================================================================

! ---> Calcul du "flux de masse" associe a la masse volumique

do ifac = 1, nfac
  wflmas(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  wflmab(ifac) = 0.d0
enddo

call cfmsfl                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   wflmas , wflmab ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , w10    , w11    , w12    ,          &
   viscf  , viscb  , coefu  ,                                     &
   ra     )

! ---> Sommation sur les faces (depend de si l'on explicite ou non
!                               le terme de convection)

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,wflmas,wflmab,w1)

do ifac = 1, nfac
  wflmas(ifac) = max(0.d0,wflmas(ifac))
enddo
do ifac = 1, nfabor
  wflmab(ifac) = max(0.d0,wflmab(ifac))
enddo
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,wflmas,wflmab,w2)


! ---> Calcul du coefficient CFL/Dt

do iel = 1, ncel
  wcf(iel) = max( -dble(iconv(ivar))*w1(iel)/volume(iel),         &
       max( dble(1-iconv(ivar))*w2(iel)/volume(iel), 0.d0 ) )
enddo

! Free memory
deallocate(viscf)
deallocate(coefu)
deallocate(w7, w8, w9)
deallocate(w10, w11, w12)

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine

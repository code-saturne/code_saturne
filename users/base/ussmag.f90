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

subroutine ussmag &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smagor , mijlij , mijmij ,                                     &
   w1     , w2     , w3     , w4     ,                            &
   ra     )

!===============================================================================
! FONCTION :
! --------

! MODIFICATION UTILISATEUR DE LA CONSTANTE DE SMAGORINSKY DE LA PHASE
! IPHAS DANS LE CAS DE L'UTILISATION D'UN MODELE DYNAMIQUE

!              SMAGOR = Mij.Lij / Mij.Mij

! EN FAIT, DES MOYENNES LOCALES DU NUMERATEUR ET DU DENOMINATEUR
! SONT REALISEES AVANT L'APPEL A USSMAG, SOIT

!              SMAGOR = < Mij.Lij > / < Mij.Mij >

! DANS CET ROUTINE, Mij.Lij ET Mij.Mij SONT PASSES EN ARGUMENT
! AVANT LA MOYENNE LOCALE.
! DANS L'EXEMPLE CI-DESSOUS ON REALISE UNE MOYENNE LOCALE DU
! RAPPORT.
!              SMAGOR = < Mij.Lij / Mij.Mij >

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iphas            ! i  ! <-- ! phase number                                   !
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
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smagor(ncelet    ! tr ! <-- ! constante de smagorinsky dans le cas           !
! , nphas)         !    !     ! d'un modlele dynamique                         !
! mijlij(ncelet    ! tr ! <-- ! mij.lij avant moyenne locale                   !
! mijmij(ncelet    ! tr ! <-- ! mij.mij avant moyenne locale                   !
! w1..4(ncelet)    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use cstnum
use optcal
use cstphy
use entsor
use parall
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iphas

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smagor(ncelet), mijlij(ncelet), mijmij(ncelet)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! --- Memoire
idebia = idbia0
idebra = idbra0

!===============================================================================
! 2.  MOYENNE SPATIALE SUR LE VOISINAGE ETENDU

!     Dans le cas ou l'utilisateur souhaite utilise le voisinage
!       etendu, il est fortement conseille de passer le mode de
!       calcul des gradients en IMRGRA = 2, afin de conserver
!       la totalite du voisinage etendu. En effet, le calcul de
!       moyenne locale est generalement degradee en voisinage reduit
!       (IMRGRA = 3).

!===============================================================================

!     On calcule le rapport
do iel = 1, ncel
  if(abs(mijmij(iel)).le.epzero) then
    w1(iel) = smagmx(iphas)**2
  else
    w1(iel) = mijlij(iel)/mijmij(iel)
  endif
enddo

!     On passe dans le filtre local
call cfiltr ( w1     , smagor , w2     , w3     )
!==========


!----
! FIN
!----

return
end subroutine

!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine lagpoi &
!================

 ( idbia0 , idbra0 ,                                              &
   lndnod ,                                                       &
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   icocel , itycel , ifrlag , itepa  ,                            &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , statis ,                                     &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     RESOLUTION DE L'EQUATION DE POISSON POUR LES VITESSE MOYENNES
!                 DES PARTICULES
!       ET CORRECTION DES VITESSES INSTANTANNEES
!                 DES PARTICULES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! icocel           ! te ! --> ! connectivite cellules -> faces                 !
! (lndnod)         !    !     !    face de bord si numero negatif              !
! itycel           ! te ! --> ! connectivite cellules -> faces                 !
! (ncelet+1)       !    !     !    pointeur du tableau icocel                  !
! ifrlag           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module lagrangien                     !
! itepa            ! te ! --> ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! moyennes statistiques                          !
!(ncelet,nvlsta    !    !     !                                                !
! w1...w3(ncel)    ! tr ! --- ! tableau de travail                             !
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
use optcal
use entsor
use cstphy
use cstnum
use pointe
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          lndnod
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          icocel(lndnod) , itycel(ncelet+1)
integer          ifrlag(nfabor) ,  itepa(nbpmax,nivep)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ifinia, ifinra
integer          npt , iel , ifac
integer          iphila , iphil
integer          iw1   , iw2   , iw3   , iw4 , iw5
integer          iw6   , iw7   , iw8   , iw9
integer          idtr   , ifmala , ifmalb
integer          iviscf , iviscb , idam   , ixam
integer          idrtp  , ismbr  , irovsd
integer          icoefap , icoefbp
integer          ivar0
integer          inc, iccocg
integer          nswrgp , imligp , iwarnp
integer          iphydp
double precision epsrgp , climgp , extrap

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idtr   = idebra
iviscf = idtr   + ncelet
iviscb = iviscf + nfac
idam   = iviscb + nfabor
ixam   = idam   + ncelet
idrtp  = ixam   + nfac*2
ismbr  = idrtp  + ncelet
irovsd = ismbr  + ncelet
ifmala = irovsd + ncelet
ifmalb = ifmala + nfac

iphila  = ifmalb + nfabor
iphil   = iphila  + ncelet
iw1    = iphil   + ncelet
iw2    = iw1    + ncelet
iw3    = iw2    + ncelet
iw4    = iw3    + ncelet
iw5    = iw4    + ncelet
iw6    = iw5    + ncelet
iw7    = iw6    + ncelet
iw8    = iw7    + ncelet
iw9    = iw8    + ncelet
ifinra = iw9    + ncelet
call rasize('lagpoi',ifinra)
!==========

do iel=1,ncel
  if ( statis(iel,ilpd) .gt. seuil ) then
    statis(iel,ilvx) = statis(iel,ilvx)                           &
                      /statis(iel,ilpd)
    statis(iel,ilvy) = statis(iel,ilvy)                           &
                      /statis(iel,ilpd)
    statis(iel,ilvz) = statis(iel,ilvz)                           &
                      /statis(iel,ilpd)
    statis(iel,ilfv) = statis(iel,ilfv)                           &
                      /( dble(npst) * volume(iel) )
  else
    statis(iel,ilvx) = 0.d0
    statis(iel,ilvy) = 0.d0
    statis(iel,ilvz) = 0.d0
    statis(iel,ilfv) = 0.d0
  endif
enddo

call lageqp                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , propce , propfa , propfb ,                            &
   ra(iviscf) , ra(iviscb) ,                                      &
   ra(idam) , ra(ixam) ,                                          &
   ra(idrtp) , ra(ismbr) , ra(irovsd) ,                           &
   ra(ifmala) , ra(ifmalb) ,                                      &
   statis(1,ilvx) , statis(1,ilvy) , statis(1,ilvz) ,             &
   statis(1,ilfv) ,                                               &
   ra(iphila) , ra(iphil) ,                                       &
   w1     , w2     , w3     , ra(iw1) , ra(iw2) ,                 &
   ra(iw3) , ra(iw4) , ra(iw5) , ra(iw6) ,                        &
   ra(iw7) , ra(iw8) , ra(iw9) ,                                  &
   ra     )

! Calcul du gradient du Correcteur PHI
! ====================================


!       On alloue localement 2 tableaux de NFABOR pour le calcul
!         de COEFA et COEFB de W1,W2,W3

icoefap = ifinra
icoefbp = icoefap + nfabor
ifinra  = icoefbp + nfabor
call rasize ('lageqp',ifinra)
!==========

do ifac = 1, nfabor
  iel = ifabor(ifac)
  ra(icoefap+ifac-1) = ra(iphil+iel-1)
  ra(icoefbp+ifac-1) = zero
enddo

inc = 1
iccocg = 1
nswrgp = 100
imligp = -1
iwarnp = 2
epsrgp = 1.d-8
climgp = 1.5d0
extrap = 0.d0


! En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(ra(iphil))
  !==========
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
ivar0 = 0

!    Sans prise en compte de la pression hydrostatique

iphydp = 0

call grdcel                                                       &
!==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   ra(iphil) , ra(iphil) , ra(iphil)    ,                         &
   ra(iphil) , ra(icoefap) , ra(icoefbp) ,                        &
   w1       , w2    , w3 ,                                        &
   ra(iw1)  , ra(iw2) , ra(iw3) ,                                 &
   ra     )

! CORRECTION DES VITESSES MOYENNES ET RETOUR AU CUMUL

do iel = 1,ncel
  if ( statis(iel,ilpd) .gt. seuil ) then
    statis(iel,ilvx) = statis(iel,ilvx) - w1(iel)
    statis(iel,ilvy) = statis(iel,ilvy) - w2(iel)
    statis(iel,ilvz) = statis(iel,ilvz) - w3(iel)
  endif
enddo

do iel = 1,ncel
  if ( statis(iel,ilpd) .gt. seuil ) then
    statis(iel,ilvx) = statis(iel,ilvx)*statis(iel,ilpd)
    statis(iel,ilvy) = statis(iel,ilvy)*statis(iel,ilpd)
    statis(iel,ilvz) = statis(iel,ilvz)*statis(iel,ilpd)
    statis(iel,ilfv) = statis(iel,ilfv)                           &
                      *( dble(npst) * volume(iel) )
  endif
enddo

! CORRECTION DES VITESSES INSTANTANNES

do npt = 1,nbpart
  if ( itepa(npt,jisor).gt.0 ) then
    iel = itepa(npt,jisor)
    ettp(npt,jup) = ettp(npt,jup) - w1(iel)
    ettp(npt,jvp) = ettp(npt,jvp) - w2(iel)
    ettp(npt,jwp) = ettp(npt,jwp) - w3(iel)
  endif
enddo

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

end subroutine

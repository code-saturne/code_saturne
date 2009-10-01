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

subroutine diverv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     ,                                                       &
   div    , ux     , vy     , wz     ,                            &
   coefax , coefay , coefaz ,                                     &
   coefbx , coefby , coefbz ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser ,                                              &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!        CALCUL DE LA DIVERGENCE D'UN VECTEUR

!   (On ne s'embete pas, on appelle 3 fois le gradient)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! div(ncelet)      ! tr ! --> ! divergence du vecteur                          !
! ux,uy,uz         ! tr ! --> ! composante du vecteur                          !
! (ncelet)         !    !     !                                                !
! coefax,...       ! tr ! ->  ! conditions aux limites pour les                !
! coefbz           !    !     ! faces de bord                                  !
! (nfabor)         !    !     !                                                !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "parall.h"
include "period.h"
include "lagpar.h"
include "lagran.h"

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
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet)
double precision div(ncelet)
double precision ux(ncelet) , vy(ncelet) , wz(ncelet)
double precision coefax(nfabor) , coefay(nfabor) , coefaz(nfabor)
double precision coefbx(nfabor) , coefby(nfabor) , coefbz(nfabor)
double precision w1(ncelet) , w2(ncelet) , w3(ncelet)
double precision w4(ncelet) , w5(ncelet) , w6(ncelet)
double precision w7(ncelet) , w8(ncelet) , w9(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ivar0
integer          iphydp,idimte,itenso
integer          iel
integer          inc, iccocg
integer          nswrgp, imligp, iwarnp
double precision epsrgp, climgp, extrap

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! En periodique et parallele, echange avant calcul du gradient

!    Parallele
if(irangp.ge.0) then
  call parcom(ux)
  !==========
  call parcom(vy)
  !==========
  call parcom(wz)
  !==========
endif

!    Periodique
if(iperio.eq.1) then
  idimte = 1
  itenso = 0
  call percom                                                     &
  !==========
( idimte , itenso ,                                               &
  ux     , ux     , ux     ,                                      &
  vy     , vy     , vy     ,                                      &
  wz     , wz     , wz     )
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0

!    Sans prise en compte de la pression hydrostatique

iphydp = 0

inc = 1
iccocg = 1
nswrgp = 100
imligp = -1
iwarnp = 2
epsrgp = 1.d-8
climgp = 1.5d0
extrap = 0.d0

!===============================================================================
! 1. Calcul du gradient de UX DANS W1
!===============================================================================

call grdcel                                                       &
!==========
( idebia , idebra ,                                               &
  ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,  &
  nnod   , lndfac , lndfbr , ncelbr , nphas  ,                    &
  nideve , nrdeve , nituse , nrtuse ,                             &
  ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,  &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                    &
  ipnfac , nodfac , ipnfbr , nodfbr ,                             &
  idevel , ituser , ia     ,                                      &
  xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,  &
  ux     , ux     , ux     ,                                      &
  ux     , coefax , coefbx ,                                      &
  w1     , w4     , w5     ,                                      &
  w6     , w7     , w8     ,                                      &
  rdevel , rtuser , ra     )

!===============================================================================
! 2. Calcul du gradient de VY DANS W2
!===============================================================================

call grdcel                                                       &
!==========
( idebia , idebra ,                                               &
  ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,  &
  nnod   , lndfac , lndfbr , ncelbr , nphas  ,                    &
  nideve , nrdeve , nituse , nrtuse ,                             &
  ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,  &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                    &
  ipnfac , nodfac , ipnfbr , nodfbr ,                             &
  idevel , ituser , ia     ,                                      &
  xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,  &
  vy     , vy     , vy     ,                                      &
  vy     , coefay , coefby ,                                      &
  w4     , w2     , w5     ,                                      &
  w6     , w7     , w8     ,                                      &
  rdevel , rtuser , ra     )

!===============================================================================
! 3. Calcul du gradient de VZ DANS W3
!===============================================================================

call grdcel                                                       &
!==========
( idebia , idebra ,                                               &
  ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,  &
  nnod   , lndfac , lndfbr , ncelbr , nphas  ,                    &
  nideve , nrdeve , nituse , nrtuse ,                             &
  ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,  &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                    &
  ipnfac , nodfac , ipnfbr , nodfbr ,                             &
  idevel , ituser , ia     ,                                      &
  xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,  &
  wz     , wz     , wz     ,                                      &
  wz     , coefaz , coefbz ,                                      &
  w5     , w6     , w3     ,                                      &
  w7     , w8     , w9     ,                                      &
  rdevel , rtuser , ra     )

!===============================================================================
! 4. Calcul de la divergence du vecteur (UX,VY,WZ)
!===============================================================================

do iel = 1,ncel
  div(iel) = w1(iel) + w2(iel) + w3(iel)
enddo

!----
! FIN
!----

end subroutine

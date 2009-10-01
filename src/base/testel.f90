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

subroutine testel &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  , nvar   ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    ,                                                       &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!-------------------------------------------------------------------------------
! Arguments
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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
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
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
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
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "entsor.h"
include "albase.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr , nphas , nvar
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtp(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra, ifinia, ifinra
integer          ifac  , iel   , ivar  , iphas
integer          inc   , iccocg, iphydp
integer          iuiph , iviph , iwiph
integer          nswrgp, imligp, iwarnp
integer          ipclip
integer          iw1   , iw2   , iw3
integer          indwri, indact, ipart, idimt, ientla, ivarpr

double precision epsrgp, climgp, extrap
double precision xx, yy, zz
double precision rbid(1)

character*32     namevr

!===============================================================================

!===============================================================================
! 0.  INITIALISATIONS
!===============================================================================

ifinia = idbia0

! On positionne l'indicateur ALE a 1 de maniere a forcer le recalcul
! de la contribution des cellules de bord a chaque appel de GRDCEL
iale = 1

iw1    = idbra0
iw2    = iw1    + ncelet
iw3    = iw2    + ncelet
ifinra = iw3    + ncelet

CALL RASIZE('TESTEL',IFINRA)
!==========

!===============================================================================
! 1. FONCTION ANALYTIQUE SIN(X+2Y+3Z)
!===============================================================================

iphas = 1
iuiph = iu(iphas)
iviph = iv(iphas)
iwiph = iw(iphas)

ivar   = ipr(iphas)
ipclip = iclrtp(ivar,icoef)

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,ivar) = sin(xx+2.d0*yy+3.d0*zz)
enddo

do ifac = 1, nfabor
  xx = cdgfbo(1,ifac)
  yy = cdgfbo(2,ifac)
  zz = cdgfbo(3,ifac)
  coefa(ifac,ipclip) = sin(xx+2.d0*yy+3.d0*zz)
enddo

do ifac = 1, nfabor
  coefb(ifac,ipclip) = 0.d0
enddo

! On active le writer standard

indwri = -1
indact = 1
call pstact(indwri, indact)
!==========

! Options de sorties des variables (gradient non entrelaces)

ipart = -1
idimt = 3
ientla = 0
ivarpr = 1

!===============================================================================
! 2. CALCUL DU GRADIENT DE LA FONCTION ANALYTIQUE

!    NE PAS CHANGER L'ORDRE DE CALCUL DES GRADIENTS:
!      * IMRGRA = 0
!      * IMRGRA = 1 (voisinage standard)
!      * IMRGRA = 2 (voisinage etendu)
!      * IMRGRA = 4 (voisinage etendu)
!      * IMRGRA = 3 (reduction du voisinage etendu)
!===============================================================================

inc = 1
iccocg = 1
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
iphydp = 0

!  2.1 APPEL A GRDCEL AVEC IMRGRA = 0
!  ==================================

imrgra = 0
imligp = -1

call grdcel                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

! On sort le gradient

NAMEVR = 'Grad_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

NAMEVR = 'Err_Grad_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)


!  2.2 APPEL A GRDCEL AVEC IMRGRA = 1
!  ==================================

imrgra = 1
imligp = 1

call grdcel                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )


! On sort le gradient

NAMEVR = 'Grad_LSQ'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

NAMEVR = 'Err_Grad_LSQ'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)


!  2.3 APPEL A GRDCEL AVEC IMRGRA = 2
!  ==================================

imrgra = 2
imligp = 1

call grdcel                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

! On sort le gradient

NAMEVR = 'Grad_LSQ_Ext'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

NAMEVR = 'Err_Grad_LSQ_Ext'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)


!  2.4 APPEL A GRDCEL AVEC IMRGRA = 4
!  ==================================

imrgra = 4
imligp = -1

call grdcel                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

! On sort le gradient

NAMEVR = 'Grad_LSQ_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

NAMEVR = 'Err_Grad_LSQ_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)


!  2.5 APPEL A GRDCEL AVEC IMRGRA = 3
!  ==================================

! Reduction du voisinage etendu

call redvse(anomax)
!==========

imrgra = 3
imligp = 1

call grdcel                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

! On sort le gradient

NAMEVR = 'Grad_LSQ_ExtRed'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

NAMEVR = 'Err_Grad_LSQ_ExtRed'
call psteva(ipart , namevr, idimt, ientla, ivarpr,                &
!==========
            ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)


!----
! FIN
!----

return
end subroutine

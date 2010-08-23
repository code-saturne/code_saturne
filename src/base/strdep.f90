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

subroutine strdep &
!================

 ( idbia0 , idbra0 , itrale , italim , itrfin ,                   &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nvar   ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   flmalf , flmalb , cofale , xprale , depale , rdevel , rtuser , &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! DEPLACEMENT DES STRUCTURES MOBILES EN ALE EN COUPLAGE INTERNE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! itrale           ! e  ! <-- ! numero d'iteration pour l'ale                  !
! italim           ! e  ! <-- ! numero d'iteration couplage implicite          !
! itrfin           ! e  ! <-- ! indicateur de derniere iteration de            !
!                  !    !     !                    couplage implicite          !
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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfabor+1)     !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr)       !    !     !  (optionnel)                                   !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! flmalf(nfac)     ! tr ! --> ! sauvegarde du flux de masse faces int          !
! flmalb(nfabor    ! tr ! --> ! sauvegarde du flux de masse faces brd          !
! cofale           ! tr ! --> ! sauvegarde des cl de p et u                    !
!    (nfabor,8)    !    !     !                                                !
! xprale(ncelet    ! tr ! --> ! sauvegarde de la pression, si nterup           !
!                  !    !     !    est >1                                      !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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
use dimens, only: ndimfb
use ihmpre
use cstphy
use numvar
use optcal
use entsor
use pointe
use albase
use alstru
use alaste
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          itrale , italim , itrfin
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr , nvar
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision depale(nnod,3)
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision flmalf(nfac), flmalb(nfabor), xprale(ncelet)
double precision cofale(nfabor,8)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra, ifinra
integer          istr, ii, iel, ifac, iphas, ntab
integer          iflmas, iflmab, iclp, iclu, iclv, iclw
integer          indast
integer          icvext, icvint, icv

double precision delta

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iflmas = ipprof(ifluma(iu(1)))
iflmab = ipprob(ifluma(iu(1)))
iclp = iclrtp(ipr(1),icoef)
iclu = iclrtp(iu(1),icoef)
iclv = iclrtp(iv(1),icoef)
iclw = iclrtp(iw(1),icoef)

!===============================================================================
! 2.  CALCUL DES EFFORTS SUR LES STRUCTURES
!===============================================================================

do istr = 1, nbstru
  do ii = 1, ndim
    forsta(ii,istr) = forstr(ii,istr)
!-a tester          FORSTA(II,ISTR) = FORSTP(II,ISTR)
    forstr(ii,istr) = 0.d0
  enddo
enddo

iforas = idebra
ifinra = iforas + 3*nbfast
call rasize('strdep',ifinra)
!==========

indast = 0
do ifac = 1, nfabor
  istr = ia(iidfst+ifac-1)
  if (istr.gt.0) then
    do ii = 1, 3
      forstr(ii,istr) = forstr(ii,istr)                           &
           + ra(iforbr+(ifac-1)*ndim+ii-1)
    enddo
  else if (istr.lt.0) then
    indast = indast + 1
    do ii = 1, 3
      ra(iforas+(indast-1)*ndim+ii-1) = asddlf(ii,-istr)          &
           * ra(iforbr+(ifac-1)*ndim+ii-1)
    enddo
  endif
enddo

if (irangp.ge.0) then
  ntab = ndim*nbstru
  call parrsm(ntab,forstr)
endif

!     Calcul de l'effort envoye au structures internes
do istr = 1, nbstru
  do ii = 1, ndim
    forstp(ii,istr) = cfopre*forstr(ii,istr)+                     &
         (1.d0-cfopre)*forsta(ii,istr)
  enddo
enddo

!     Envoi de l'effort applique aux structures externes
if (nbaste.gt.0) then
  call astfor(ntcast, nbfast, ra(iforas))
  !==========
endif

!     Si on est en phase d'initialisation du fluide
if (itrale.le.nalinf) then
  itrfin = -1
  return
endif

!===============================================================================
! 3.  DEFINITION UTILISATEUR DES CARACTERISTIQUES DE LA STRUCTURE
!===============================================================================


if (nbstru.gt.0) then

  ! - Interface Code_Saturne
  !   ======================

  if (iihmpr.eq.1) then

    call uistr2 &
    !==========
 ( xmstru, xcstru, xkstru,     &
   forstr,                     &
   dtref, ttcabs, ntcabs   )

  endif

  call usstr2                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nbstru ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iidfst),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     ,                                                       &
   xmstru , xcstru , xkstru , xstreq , xstr   , xpstr  , forstp , &
   dtstr  ,                                                       &
   rdevel , rtuser ,                                              &
   ra     )

endif

!===============================================================================
! 4.  DEPLACEMENT DES STRUCTURES INTERNES
!===============================================================================


do istr = 1, nbstru

  call newmrk                                                     &
  !==========
 ( istr  , alpnmk  , betnmk          , gamnmk          ,          &
   xmstru(1,1,istr), xcstru(1,1,istr), xkstru(1,1,istr),          &
   xstreq(1,istr)  ,                                              &
   xstr(1,istr)    , xpstr(1,istr)   , xppstr(1,istr)  ,          &
   xsta(1,istr)    , xpsta(1,istr)   , xppsta(1,istr)  ,          &
   forstp(1,istr)  , forsta(1,istr)  , dtstr(istr)     )

enddo

!===============================================================================
! 5.  TEST DE CONVERGENCE
!===============================================================================

icvext = 0
icvint = 0
icv    = 0

iphas = 1

delta = 0.d0
do istr = 1, nbstru
  do ii = 1, 3
    delta = delta + (xstr(ii,istr)-xstp(ii,istr))**2
  enddo
enddo
if (nbstru.gt.0) then
  delta = sqrt(delta)/almax(iphas)/nbstru
  if (delta.lt.epalim) icvint = 1
endif

if (nbaste.gt.0) call astcv1(ntcast, icvext)
                 !==========


if (nbstru.gt.0.and.nbaste.gt.0) then
   icv = icvext*icvint
elseif (nbstru.gt.0.and.nbaste.eq.0) then
   icv = icvint
elseif (nbaste.gt.0.and.nbstru.eq.0) then
   icv = icvext
endif

if (iwarni(iuma).ge.2) write(nfecra,1000) italim, delta

!     si convergence
if (icv.eq.1) then
  if (itrfin.eq.1) then
!       si ITRFIN=1 on sort
    if (iwarni(iuma).ge.1) write(nfecra,1001) italim, delta
    itrfin = -1
  else
!       sinon on refait une derniere iteration pour SYRTHES/T1D/rayonnement
!        et on remet ICV a 0 pour que Code_Aster refasse une iteration aussi
    itrfin = 1
    icv = 0
  endif
elseif (itrfin.eq.0 .and. italim.eq.nalimx-1) then
!       ce sera la derniere iteration
  itrfin = 1
elseif (italim.eq.nalimx) then
!       on a forcement ITRFIN=1 et on sort
  if (nalimx.gt.1) write(nfecra,1100) italim, delta
  itrfin = -1
!       On met ICV a 1 pour que Code_Aster s'arrete lui aussi
  icv = 1
endif

!     On renvoie l'indicateur de convergence final a Code_Aster
call astcv2(ntcast, icv)
!==========

!===============================================================================
! 6.  RETOUR AUX VALEURS ANTERIEURES SI NECESSAIRE
!===============================================================================

!     Si NTERUP    .GT.1, RTPA a ete touche apres NAVSTO, on doit donc
!       revenir a une valeur anterieure
if (itrfin.ne.-1) then
  do ii = 1, nvar
    if (ii.eq.ipr(1) .and. nterup.gt.1) then
      do iel = 1, ncelet
        rtpa(iel,ii) = xprale(iel)
      enddo
    endif
    do iel = 1, ncelet
      rtp(iel,ii) = rtpa(iel,ii)
    enddo
  enddo
  do ifac = 1, nfac
     propfa(ifac,iflmas) = flmalf(ifac)
  enddo
  do ifac = 1, nfabor
     propfb(ifac,iflmab) = flmalb(ifac)
     coefa(ifac,iclp) = cofale(ifac,1)
     coefa(ifac,iclu) = cofale(ifac,2)
     coefa(ifac,iclv) = cofale(ifac,3)
     coefa(ifac,iclw) = cofale(ifac,4)
     coefb(ifac,iclp) = cofale(ifac,5)
     coefb(ifac,iclu) = cofale(ifac,6)
     coefb(ifac,iclv) = cofale(ifac,7)
     coefb(ifac,iclw) = cofale(ifac,8)
  enddo
endif

!----
! FORMATS
!----

#if defined(_CS_LANG_FR)

 1000 format (                                                          &
 '            ALE IMPLICITE : ITER=',I5,' DERIVE=',E12.5     )
 1001 format (                                                          &
 'CONVERGENCE ALE IMPLICITE : ITER=',I5,' DERIVE=',E12.5     )
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ ATTENTION : COUPLAGE IMPLICITE ALE                      ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@  Derive normee :',E12.5                                     ,/,&
'@                                                            '  )

#else

 1000 format (                                                          &
 '            IMPLICIT ALE: ITER=',I5,' DERIVE=',E12.5     )
 1001 format (                                                          &
 'CONVERGENCE IMPLICIT ALE: ITER=',I5,' DERIVE=',E12.5     )
 1100 format (                                                          &
'@                                                            ',/,&
'@ @@ WARNING: IMPLICIT ALE                                   ',/,&
'@    ========                                                ',/,&
'@  Maximum number of iterations ',I10   ,' reached           ',/,&
'@  Normed derive :',E12.5                                     ,/,&
'@                                                            '  )

#endif

!----
! FIN
!----

end subroutine

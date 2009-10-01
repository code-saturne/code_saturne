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

subroutine reseps &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itpsmp ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , grdvit , produc , grarox , graroy , graroz , &
   ckupdc , smcelp , gamma  ,                                     &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR EPSILON DANS LE CAS DU MODELE Rij-epsilon

!   On a ici ISOU   = 7 (7 ieme variable du Rij-epsilon)

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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! ncepdp           ! e  ! <-- ! nombre de cellules avec pdc                    !
! ncesmp           ! e  ! <-- ! nombre de cellules a source de masse           !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! numero de phase                                !
! ivar             ! e  ! <-- ! numero de variable                             !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itpsmp           ! te ! <-- ! type de source de masse pour la                !
! (ncesmp)         !    !     !  variables (cf. ustsma)                        !
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
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! grdvit           ! tr ! --- ! tableau de travail pour terme grad             !
!  (ncelet,3,3)    !    !     !    de vitesse     uniqt pour iturb=31          !
! produc           ! tr ! <-- ! tableau de travail pour production             !
!  (6,ncelet)      !    !     ! (sans rho volume) uniqt pour iturb=30          !
! grarox,y,z       ! tr ! <-- ! tableau de travail pour grad rom               !
!  (ncelet)        !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smcelp(ncesmp    ! tr ! <-- ! valeur de la variable associee a la            !
!                  !    !     !  source de masse                               !
! gamma(ncesmp)    ! tr ! <-- ! valeur du flux de masse                        !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!  (ncelet,*)      !    !     !   lagrangien                                   !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
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

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstnum.h"
include "cstphy.h"
include "pointe.h"
include "period.h"
include "parall.h"
include "lagpar.h"
include "lagran.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , ivar   , isou   , ipp

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itpsmp(ncesmp)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision produc(6,ncelet), grdvit(ncelet,3,3)
double precision grarox(ncelet), graroy(ncelet), graroz(ncelet)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision viscf(nfac), viscb(nfabor)
double precision tslagr(ncelet,*)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra, ifinia
integer          init  , ifac  , iel   , inc   , iccocg
integer          ii    ,jj     , iiun
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          ieiph , iuiph
integer          iclvar, iclvaf
integer          ipcrom, ipcroo, ipcvis, ipcvst, iflmas, iflmab
integer          nswrgp, imligp, iwarnp, iphydp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          idimte, itenso
integer          iptsta
integer          maxelt, ils
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision trprod , trrij ,csteps, rctse
double precision grdpx,grdpy,grdpz,grdsn
double precision surfn2
double precision tseps , kseps , ceps2
double precision tuexpe, thets , thetv , thetap, thetp1

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) nomvar(ipp)
endif

iuiph  = iu(iphas)
ir11ip = ir11(iphas)
ir22ip = ir22(iphas)
ir33ip = ir33(iphas)
ir12ip = ir12(iphas)
ir13ip = ir13(iphas)
ir23ip = ir23(iphas)
ieiph  = iep (iphas)

ipcrom = ipproc(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
iflmas = ipprof(ifluma(iuiph))
iflmab = ipprob(ifluma(iuiph))

iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

!     Constante Ce2, qui vaut CE2 pour ITURB=30 et CSSGE2 pour ITRUB=31
if (iturb(iphas).eq.30) then
  ceps2 = ce2
else
  ceps2 = cssge2
endif

!     S pour Source, V pour Variable
thets  = thetst(iphas)
thetv  = thetav(ivar )

ipcroo = ipcrom
if(isto2t(iphas).gt.0.and.iroext(iphas).gt.0) then
  ipcroo = ipproc(iroma(iphas))
endif
if(isto2t(iphas).gt.0) then
  iptsta = ipproc(itstua(iphas))
else
  iptsta = 0
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 2. TERMES SOURCES  UTILISATEURS
!===============================================================================

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
CALL IASIZE('RESEPS',IFINIA)

call ustsri                                                       &
!==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itpsmp ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smcelp , gamma  , grdvit , produc , &
   smbr   , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , dam    , drtp   ,                   &
   rdevel , rtuser , ra     )

!     Si on extrapole les T.S.
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
!       Sauvegarde pour echange
    tuexpe = propce(iel,iptsta+isou-1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+isou-1) = smbr(iel)
!       Second membre du pas de temps precedent
!       On suppose -ROVSDT > 0 : on implicite
!          le terme source utilisateur (le reste)
    smbr(iel) = rovsdt(iel)*rtpa(iel,ivar) - thets*tuexpe
!       Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel)   = rovsdt(iel)*rtpa(iel,ivar) + smbr(iel)
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

!===============================================================================
! 3.  TERMES SOURCES LAGRANGIEN : COUPLAGE RETOUR
!===============================================================================

!     Ordre 2 non pris en compte
 if (iilagr.eq.2 .and. ltsdyn.eq.1 .and. iphas.eq.ilphas) then

   do iel = 1,ncel
! Ts sur eps
     tseps = -0.5d0 * ( tslagr(iel,itsr11)                        &
                      + tslagr(iel,itsr22)                        &
                      + tslagr(iel,itsr33) )
! rapport k/eps
     kseps = 0.5d0 * ( rtpa(iel,ir11ip)                           &
                     + rtpa(iel,ir22ip)                           &
                     + rtpa(iel,ir33ip) )                         &
                     / rtpa(iel,ieiph)

     smbr(iel)   = smbr(iel) + ce4 *tseps *rtpa(iel,ieiph) /kseps
     rovsdt(iel) = rovsdt(iel) + max( (-ce4*tseps/kseps) , zero)
   enddo

 endif

!===============================================================================
! 4. TERME SOURCE DE MASSE
!===============================================================================


if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t(iphas) , thetv ,    &
   icetsm , itpsmp ,                                              &
   volume , rtpa(1,ivar) , smcelp , gamma  ,                      &
   smbr   ,  rovsdt , w1 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t(iphas).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
      propce(iel,iptsta+isou-1) + w1(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w1(iel)
    enddo
  endif

endif

!===============================================================================
! 5. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
!    ET TERME INSTATIONNAIRE
!===============================================================================

! ---> Calcul de mij

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

! ---> Ajout au second membre

do iel = 1, ncel
  smbr(iel) = smbr(iel)                                           &
              + iconv(ivar)*w1(iel)*rtpa(iel,ivar)
enddo

! ---> Ajout dans la diagonale de la matrice
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
           + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel) &
           - iconv(ivar)*w1(iel)*thetv
enddo


!===============================================================================
! 6. PRODUCTION RHO * Ce1 * epsilon / k * P
!    DISSIPATION RHO*Ce2.epsilon/k*epsilon
!===============================================================================


! ---> Calcul de k pour la suite du sous-programme
!       on utilise un tableau de travail puisqu'il y en a...
do iel = 1, ncel
  w8(iel) = 0.5d0 * (rtpa(iel,ir11ip) + rtpa(iel,ir22ip) +        &
                     rtpa(iel,ir33ip))
enddo
! ---> Calcul de la trace de la production, suivant qu'on est en
!     Rij standard ou en SSG (utilisation de PRODUC ou GRDVIT)
if (iturb(iphas).eq.30) then
  do iel = 1, ncel
    w9(iel) = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
  enddo
else
  do iel = 1, ncel
    w9(iel) = -( rtpa(iel,ir11ip)*grdvit(iel,1,1) +               &
                 rtpa(iel,ir12ip)*grdvit(iel,1,2) +               &
                 rtpa(iel,ir13ip)*grdvit(iel,1,3) +               &
                 rtpa(iel,ir12ip)*grdvit(iel,2,1) +               &
                 rtpa(iel,ir22ip)*grdvit(iel,2,2) +               &
                 rtpa(iel,ir23ip)*grdvit(iel,2,3) +               &
                 rtpa(iel,ir13ip)*grdvit(iel,3,1) +               &
                 rtpa(iel,ir23ip)*grdvit(iel,3,2) +               &
                 rtpa(iel,ir33ip)*grdvit(iel,3,3) )
  enddo
endif


!     Terme explicite (Production)

do iel = 1, ncel
!     Demi-traces
  trprod = w9(iel)
  trrij  = w8(iel)
  w1(iel)   =             propce(iel,ipcroo)*volume(iel)*         &
       ce1*rtpa(iel,ieiph)/trrij*trprod
enddo

!     Si on extrapole les T.S : PROPCE
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
    propce(iel,iptsta+isou-1) =                                   &
    propce(iel,iptsta+isou-1) + w1(iel)
  enddo
!     Sinon : SMBR
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w1(iel)
  enddo
endif

!     Terme implicite (Dissipation)
do iel = 1, ncel
  trrij  = w8(iel)
  smbr(iel) = smbr(iel) - propce(iel,ipcrom)*volume(iel)*         &
                         ceps2*rtpa(iel,ieiph)**2/trrij
enddo

! ---> Matrice

if(isto2t(iphas).gt.0) then
  thetap = thetv
else
  thetap = 1.d0
endif
do iel = 1, ncel
  trrij  = w8(iel)
  rovsdt(iel) = rovsdt(iel) + ceps2*propce(iel,ipcrom)*volume(iel)&
                     *rtpa(iel,ieiph)/trrij*thetap
enddo

!===============================================================================
! 7. TERMES DE GRAVITE
!===============================================================================

if(igrari(iphas).eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijthe                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   , isou   , ipp    ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , grarox , graroy , graroz , w7     ,          &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!     Si on extrapole les T.S. : PROPCE
if(isto2t(iphas).gt.0) then
  do iel = 1, ncel
     propce(iel,iptsta+isou-1) =                                  &
     propce(iel,iptsta+isou-1) + w7(iel)
   enddo
!     Sinon SMBR
 else
   do iel = 1, ncel
     smbr(iel) = smbr(iel) + w7(iel)
   enddo
 endif

endif


!===============================================================================
! 8.a TERMES DE DIFFUSION  A.grad(Eps) : PARTIE EXTRADIAGONALE EXPLICITE
!     RIJ STANDARD
!===============================================================================

if (iturb(iphas).eq.30) then

! ---> Calcul du grad(Eps)


  iccocg = 1
  inc = 1

  nswrgp = nswrgr(ivar )
  imligp = imligr(ivar )
  iwarnp = iwarni(ivar )
  epsrgp = epsrgr(ivar )
  climgp = climgr(ivar )
  extrap = extrag(ivar )
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,ivar )   , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

! ---> Calcul des termes extradiagonaux de A.grad(Eps)

  do iel = 1 , ncel
    trrij  = w8(iel)
    csteps  = propce(iel,ipcroo) * crijep *trrij / rtpa(iel,ieiph)
    w4(iel) = csteps * ( rtpa(iel,ir12ip) * w2(iel) +             &
         rtpa(iel,ir13ip) * w3(iel) )
    w5(iel) = csteps * ( rtpa(iel,ir12ip) * w1(iel) +             &
         rtpa(iel,ir23ip) * w3(iel) )
    w6(iel) = csteps * ( rtpa(iel,ir13ip) * w1(iel) +             &
         rtpa(iel,ir23ip) * w2(iel) )
  enddo

! ---> Assemblage de { A.grad(Eps) } .S aux faces

  call vectds                                                     &
  !==========
 (ndim   , ncelet , ncel   , nfac   , nfabor   ,                  &
  ifacel , ifabor , ia     ,                                      &
  surfac , surfbo , ra(ipond) ,                                   &
  w4     , w5     , w6     ,                                      &
  viscf  , viscb  , ra     )

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                                   ifacel,ifabor,viscf,viscb,w4)

!     Si on extrapole les termes sources
  if(isto2t(iphas).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
           propce(iel,iptsta+isou-1) + w4(iel)
    enddo
!     Sinon
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w4(iel)
    enddo
  endif


!===============================================================================
! 8.b TERMES DE DIFFUSION  A.grad(Eps) : PARTIE DIAGONALE
!     RIJ STANDARD
!===============================================================================
!     Implicitation de (grad(eps).n)n en gradient facette
!     Si IDIFRE=1, terme correctif explicite
!        grad(eps)-(grad(eps).n)n calcule en gradient cellule
!     Les termes de bord sont uniquement pris en compte dans la partie
!        en (grad(eps).n)n
!     (W1,W2,W3) contient toujours le gradient de la variable traitee

!     La parcom/percom-isation du gradient de epsilon a ete faite dans
!       grdcel. Pas utile de recommencer.

  if (idifre(iphas).eq.1) then

    do iel = 1, ncel
      trrij  = w8(iel)
      csteps = propce(iel,ipcroo) * crijep *trrij/rtpa(iel,ieiph)
      w4(iel)=csteps*rtpa(iel,ir11ip)
      w5(iel)=csteps*rtpa(iel,ir22ip)
      w6(iel)=csteps*rtpa(iel,ir33ip)
    enddo

! --->  TRAITEMENT DU PARALLELISME (MEMES DOUTES QUE PERIODICITE ?)

    if(irangp.ge.0) then
      call parcom (w4)
      !==========
      call parcom (w5)
      !==========
      call parcom (w6)
      !==========
    endif

! -->   TRAITEMENT DE LA PERIODICITE
    if(iperio.eq.1) then

      idimte = 21
      itenso = 0

      call percom                                                 &
      !==========
    ( idimte , itenso ,                                           &
      w4     , w4     , w4    ,                                   &
      w5     , w5     , w5    ,                                   &
      w6     , w6     , w6    )

    endif

    do ifac = 1, nfac

      ii=ifacel(1,ifac)
      jj=ifacel(2,ifac)

      surfn2 = ra(isrfan-1+ifac)**2

      grdpx=0.5d0*(w1(ii)+w1(jj))
      grdpy=0.5d0*(w2(ii)+w2(jj))
      grdpz=0.5d0*(w3(ii)+w3(jj))
      grdsn=grdpx*surfac(1,ifac)+grdpy*surfac(2,ifac)             &
           +grdpz*surfac(3,ifac)
      grdpx=grdpx-grdsn*surfac(1,ifac)/surfn2
      grdpy=grdpy-grdsn*surfac(2,ifac)/surfn2
      grdpz=grdpz-grdsn*surfac(3,ifac)/surfn2

      viscf(ifac)= 0.5d0*(                                        &
           (w4(ii)+w4(jj))*grdpx*surfac(1,ifac)                   &
           +(w5(ii)+w5(jj))*grdpy*surfac(2,ifac)                  &
           +(w6(ii)+w6(jj))*grdpz*surfac(3,ifac))

    enddo

    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

    init = 1
    call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,              &
       ifacel,ifabor,viscf,viscb,w1)

!     Si on extrapole les termes sources
    if(isto2t(iphas).gt.0) then
      do iel = 1, ncel
        propce(iel,iptsta+isou-1) =                               &
             propce(iel,iptsta+isou-1) + w1(iel)
      enddo
!     Sinon
    else
      do iel = 1, ncel
        smbr(iel) = smbr(iel) + w1(iel)
      enddo
    endif

  endif

! ---> Viscosite orthotrope pour partie implicite

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      trrij  = w8(iel)
      rctse = propce(iel,ipcrom) * crijep * trrij/rtpa(iel,ieiph)
      w1(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*rctse*rtpa(iel,ir11ip)
      w2(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*rctse*rtpa(iel,ir22ip)
      w3(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*rctse*rtpa(iel,ir33ip)
    enddo

    call visort                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w2     , w3     ,                                     &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

else
!===============================================================================
! 8.c TERMES DE DIFFUSION
!     RIJ SSG
!===============================================================================
! ---> Viscosite

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigmae
    enddo

    call viscfa                                                   &
   !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     ,                                                       &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

endif


!===============================================================================
! 9. RESOLUTION
!===============================================================================

if(isto2t(iphas).gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+isou-1)
  enddo
endif

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ndircp = ndircl(ivar)
ireslp = iresol(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
!MO      IPP    =
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! 10. IMPRESSIONS
!===============================================================================


!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000 format(/,'           SOLVING VARIABLE ',A8           ,/)

#endif

!12345678 : MAX: 12345678901234 MIN: 12345678901234 NORM: 12345678901234
!----
! FIN
!----

return

end subroutine

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

subroutine lwctss &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   izfppp , idevel , ituser , ia     ,                            &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt ,                                              &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME PREMELANGE MODELE LWC
!   ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC.F

! ON RESOUT ROVSDT*D(VAR) = SMBRS

! ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
!  UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
!  ECRASER

! POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
!  QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
!  SMBRS

! DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
! ECRIRE :
!          SMBRS  = SMBRS  + CEXP + CIMP*VAR
!          ROVSDT = ROVSDT + MAX(-CIMP,ZERO)

! ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
!    SMBRS en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!          pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    ROVSDT en kg /s


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
! iscal            ! e  ! <-- ! numero du scalaire                             !
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
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
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
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ivar, iel, iphas, idirac, ivar0
integer          iphydp , itenso , idimte
integer          inc , iccocg
integer          ipcvst
integer          ipcrom, ii

integer          iptscl(ndracm), ipfmal(ndracm)
integer          ipfmel(ndracm), iprhol(ndracm)
double precision sum, epsi
double precision tsgrad, tschim, tsdiss

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0
epsi   = 1.0d-10

! --- Numero du scalaire a traiter : ISCAL

! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! ---
ipcrom = ipproc(irom(iphas))
ipcvst = ipproc(ivisct(iphas))

! --- Numero des grandeurs physiques (voir usclim)
do idirac = 1, ndirac
  iptscl(idirac) = ipproc(itscl(idirac))
  ipfmal(idirac) = ipproc(ifmal(idirac))
  ipfmel(idirac) = ipproc(ifmel(idirac))
  iprhol(idirac) = ipproc(irhol(idirac))
enddo

!===============================================================================
! 2. PRISE EN COMPTE DES TERMES SOURCES
!===============================================================================

if ( ivar.eq.isca(iyfm) ) then

! ---> Terme source pour la fraction massique moyenne de fuel

  do iel = 1, ncel
      sum = zero
      do idirac = 1, ndirac
        sum  = sum + propce(iel,iprhol(idirac))                   &
           *propce(iel,iptscl(idirac))*volume(iel)
      enddo

! terme implicite

      if (rtpa(iel,ivar).gt.epsi) then
        rovsdt(iel) = rovsdt(iel) + max(-sum/rtpa(iel,ivar),zero)
      endif

! terme explicite

       smbrs(iel) =  smbrs(iel) + sum

  enddo

endif

! ---> Terme source pour la variance de la fraction massique moyenne de fuel

if (ivar.eq.isca(iyfp2m)) then

  do iel = 1, ncel
    sum = zero
    do idirac = 1, ndirac
      sum  = sum + (propce(iel,iptscl(idirac))*volume(iel)        &
        *(propce(iel,ipfmal(idirac)) - rtpa(iel,isca(iyfm)))      &
             *propce(iel,iprhol(idirac)))
    enddo
    smbrs(iel) = smbrs(iel) + sum
  enddo

endif

! ---> Terme source pour la covariance

if ( ivar.eq.isca(icoyfp)) then

! --- Calcul du gradient de F
!     =======================

  ii = isca(ifm)
  do iel = 1, ncel
    w10(iel) = rtpa(iel,ii)
  enddo

! En periodique et parallele, echange avant calcul du gradient

!    Parallele
  if(irangp.ge.0) then
    call parcom(w10)
    !==========
  endif

!    Periodique
  if(iperio.eq.1) then
    idimte = 0
    itenso = 0
    call percom                                                   &
    !==========
  ( idimte , itenso ,                                             &
    w10     , w10     , w10     ,                                 &
    w10     , w10     , w10     ,                                 &
    w10     , w10     , w10     )
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0
  iphydp = 0
  inc = 1
  iccocg = 1

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar0  , imrgra , inc    , iccocg , nswrgr(ii) , imligr(ii) ,  &
   iphydp , iwarni(ii) , nfecra ,                                 &
   epsrgr(ii) , climgr(ii) , extrag(ii) ,                         &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w10    , w10    , w10    ,                                     &
   w10    , coefa(1,iclrtp(ii,icoef))  ,                          &
            coefb(1,iclrtp(ii,icoef))  ,                          &
   w1              , w2              , w3     ,                   &
!        d./dx1          , d./dx2          , d./dx3 ,
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

! --- Calcul du gradient de Yfuel
!     ===========================

  ii = isca(iyfm)
  do iel = 1, ncel
    w11(iel) = rtpa(iel,ii)
  enddo

! En periodique et parallele, echange avant calcul du gradient

!    Parallele
  if(irangp.ge.0) then
    call parcom(w11)
    !==========
  endif

!    Periodique
  if(iperio.eq.1) then
    idimte = 0
    itenso = 0
    call percom                                                   &
    !==========
  ( idimte , itenso ,                                             &
    w11     , w11     , w11     ,                                 &
    w11     , w11     , w11     ,                                 &
    w11     , w11     , w11     )
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0
  iphydp = 0
  inc = 1
  iccocg = 1

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar0  , imrgra , inc    , iccocg , nswrgr(ii) , imligr(ii) ,  &
   iphydp , iwarni(ii) , nfecra ,                                 &
   epsrgr(ii) , climgr(ii) , extrag(ii) ,                         &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w11    , w11    , w11    ,                                     &
   w11    , coefa(1,iclrtp(ii,icoef))  ,                          &
            coefb(1,iclrtp(ii,icoef))  ,                          &
   w7              , w8              , w9     ,                   &
!        d./dx1          , d./dx2          , d./dx3 ,
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


! --- Calcul du terme source
!     ======================


! ---> Calcul de K et Epsilon en fonction du modele de turbulence


! ---- TURBULENCE

  if (itytur(iphas).eq.2) then

    do iel = 1, ncel
      w10(iel) = rtpa(iel,ik(iphas))
      w11(iel) = rtpa(iel,iep(iphas))
    enddo

  elseif (itytur(iphas).eq.3) then

    do iel = 1, ncel
      w10(iel) = ( rtpa(iel,ir11(iphas))                          &
                  +rtpa(iel,ir22(iphas))                          &
                  +rtpa(iel,ir33(iphas)) ) / 2.d0
      w11(iel) = rtpa(iel,iep(iphas))
    enddo

  elseif (iturb(iphas).eq.50) then

    do iel = 1, ncel
      w10(iel) = rtpa(iel,ik(iphas))
      w11(iel) = rtpa(iel,iep(iphas))
    enddo

  elseif (iturb(iphas).eq.60) then

    do iel = 1, ncel
      w10(iel) = rtpa(iel,ik(iphas))
      w11(iel) = cmu*rtpa(iel,ik(iphas))*rtpa(iel,iomg(iphas))
    enddo

  endif

  do iel=1,ncel

!  A confirmer :
!   Le terme de dissipation devrait etre implicite
!   Dans le terme de dissipation, il manque une constante Cf
!   Peut-elle etre consideree egale a 1 ?
!   Verifier le signe du terme de production
!-
! terme implicite


    w11(iel) = w11(iel)/(w10(iel)*rvarfl(iscal))                  &
         *volume(iel)*propce(iel,ipcrom)
    rovsdt(iel) = rovsdt(iel) + max(w11(iel),zero)

! terme de gradient

    tsgrad =  (2.0d0                                              &
         * propce(iel,ipcvst)/(sigmas(iscal))                     &
         *(w1(iel)*w7(iel)+w2(iel)*w8(iel)+w3(iel)*w9(iel)))      &
         *volume(iel)


! terme de dissipation

    tsdiss = -w11(iel) * rtpa(iel,ivar)

! terme de chimique

    tschim = zero
    do idirac = 1, ndirac
      tschim =   tschim                                           &
           + (propce(iel,iptscl(idirac))                          &
           *(propce(iel,ipfmel(idirac))-rtpa(iel,isca(ifm)))      &
           *volume(iel))*propce(iel,iprhol(idirac))
    enddo

! --> Somme des termes

    smbrs(iel) = smbrs(iel) + tschim + tsgrad + tsdiss

   enddo

 endif

!----
! FIN
!----

return

end subroutine

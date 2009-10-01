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

subroutine scalai &
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
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   tslagr , coefa  , coefb  ,                                     &
   dtr    , viscf  , viscb  ,                                     &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR LES SCALAIRES SUR UN PAS DE TEMPS

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
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! dtr(ncelet)      ! tr ! --- ! dt*cdtvar                                      !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
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

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "elincl.h"

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
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision tslagr(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision dtr(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbrs(ncelet)
double precision rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra, ifinia
integer          iscal, ivar, iphas, iel
integer          ii, iisc, itspdv, icalc, iappel
integer          ispecf
integer          maxelt, ils

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

ipass = ipass + 1

!===============================================================================
! 2. TRAITEMENT DES SCALAIRES A PHYSIQUE PARTICULIERE
!    Cette section sera dediee aux traitement particuliers.
!    On donne un exemple qui n'est que la recopie (a ISCAPP pres) de
!      la section 3 sur les scalaires utilisateurs.
!===============================================================================

if (ippmod(iphpar).ge.1) then

! ---> Initialisation des RTP  a partir des CL

!   On initialise RTP (et pas RTPA) car RTPA ne sert pas
!      dans codits.

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
CALL IASIZE('SCALAI',IFINIA)

  call ppinv2                                                     &
  !==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   rdevel , rtuser , ra     )

!     On pourra eviter les bugs en initialisant aussi RTPA (sur NCELET)
!     (d'ailleurs, on s'en sert ensuite dans le cpflux etc)
  if(ipass.eq.1.and.isuite.eq.0) then
    if(nscapp.gt.0) then
      do ii = 1, nscapp
        iscal = iscapp(ii)
        ivar  = isca(iscal)
        do iel = 1, ncelet
          rtpa (iel,ivar) = rtp (iel,ivar)
        enddo
      enddo
    endif
  endif

! ---> Calculs TS relatifs a la physique du charbon
!      GMDEV1, GMDEV2, GMHET, GMDCH

  if ( ippmod(icp3pl).ne.-1 ) then

    call cpflux                                                   &
    !==========
   ( idebia , idebra , ncelet , ncel   ,                          &
     rtpa   , propce , volume ,                                   &
     w1     , w2     , w3     ,                                   &
     ra     )

  endif

! ---> Calculs TS relatifs a la physique du fuel
!      GAMEVA, GAMHTF

  if ( ippmod(icfuel).ne.-1 ) then

    call fuflux                                                   &
    !==========
   ( idebia , idebra , ncelet , ncel   ,                          &
     rtpa   , propce , volume ,                                   &
     w1     , w2     , w3     ,                                   &
     ra     )

  endif

!    ATTENTION : POUR LE CLIPPING AVEC ICLP = 1, IL FAUT

!                QUE LES SCALAIRES SOIENT RESOLUS AVANT
!                LEUR VARIANCE ASSOCIEE






! ---> Boucle sur les scalaires physiques particulieres.
!      On peut imaginer a la place des resolutions couplees.
!      Ici, on ne donne qu'un exemple.

  do ii = 1, nscapp

    iscal = iscapp(ii)
    ivar  = isca(iscal)

! ---> Pas de temps (avec facteur multiplicatif eventuel)

    if(cdtvar(ivar).ne.1.d0) then
      do iel = 1, ncel
        dtr(iel) = dt(iel)*cdtvar(ivar)
      enddo
    else
      do iel = 1, ncel
        dtr(iel) = dt(iel)
      enddo
    endif

!     Schema compressible sans choc :
! ---> Traitement special pour la masse volumique,
!                     la temperature et l'energie
!     L'indicateur ISPECF sera non nul si l'on ne doit pas resoudre
!       le scalaire plus bas avec covofi.

    ispecf = 0

    if ( ippmod(icompf).ge.0 ) then

      do iphas = 1, nphas
        if ( iscal.eq.irho(iphas) .or.                            &
             iscal.eq.itempk(iphas) ) then
          ispecf = 1
        elseif ( iscal.eq.ienerg(iphas) ) then
          ispecf = 2
        endif
      enddo

! ---> Masse volumique : deja resolue
! ---> Temperature     : n'est pas une variable transportee
! ---> Enthalpie       : a resoudre

      if(ispecf.eq.2) then

        iphas = iphsca(iscal)

        call cfener                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas)   , ncetsm(iphas)   ,                            &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) , ia(iicesm(iphas)) , ia(iitpsm(iphas)) ,    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ra(ickupd(iphas)) , ra(ismace(iphas)) ,      &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

      endif

    endif

!     Pour le compressible, on ne resout pas celles deja resolues ou
!       qui ne doivent pas l'etre
!     Pour les autres physiques, on resout les scalaires dans l'ordre
!       (pas de tests a faire)

    if(ispecf.eq.0) then

! ---> Variances et scalaires

!   iscavr dira scalaire/variance
!   itspdv dira si calcul des termes de prod et dissip supplementaires
!         ou non (1 : oui, 0 : non)

!         si iscavr = 0
!           scalaire
!           itspdv = 0
!         sinon
!           variance
!           si iscavr > 0 et iscavr < nscal+1
!             itspdv = 1
!           sinon
!             pour le moment, on s'arrete
!             a terme, les combustionnistes pourront donner leur propre
!               grandeur scalaire associee a la variance et
!               eventuellement reconstruite en dehors de covofi.
!               il faudra alors peut etre declarer des tableaux
!               de travail suppl
!           fin si
!         fin si

      iisc = iscal
      if(iscavr(iisc).eq.0) then
        itspdv = 0
      elseif(iscavr(iisc).gt.0.and.iscavr(iisc).le.nscal) then
        itspdv = 1
      else
        write(nfecra,9000)iisc,iisc,nscal,iscavr(iisc)
        call csexit (1)
      endif


! ---> Appel a covofi pour la resolution

      iphas  = iphsca(iisc)

      call covofi                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iisc   , itspdv ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) , ia(iicesm(iphas)) , ia(iitpsm(iphas)) ,    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dtr    , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ra(ickupd(iphas)) , ra(ismace(iphas)) ,      &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )


! ---> Versions Electriques
!             Effet Joule
!             Arc Electrique
!             Conduction ionique

!     On calcule ici j, E, j.E reels et imagimaires

      if ( ippmod(ieljou).ge.1 .or.                               &
           ippmod(ielarc).ge.1 .or.                               &
           ippmod(ielion).ge.1       ) then


!     On utilise le  fait que les scalaires sont dans l'ordre
!       H, PotR, [PotI], [A] pour faire le calcul de j, E et j.E
!       apres  la determination de PotR [et PotI].

        icalc = 0
!            On peut y aller apres PotR si on est en arc
!                                         ou en Joule sans PotI

        if(ippmod(ielarc).ge.1.or.ippmod(ieljou).eq.1             &
             .or.ippmod(ieljou).eq.3) then
          if(iscal.eq.ipotr) then
            icalc = 1
          endif
        endif
!     On y va apres PotI si on est en Joule avec PotI
        if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
          if(iscal.eq.ipoti) then
            icalc = 1
          endif
        endif

        if(icalc.eq.1) then

!     Calcul de j, E et j.E
          iappel = 1

          call elflux                                             &
          !==========
 ( idebia , idebra , iappel ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )


!     Recalage des variables electriques j, j.E (et Pot, E)

          if ( ielcor .eq.1  .and. ntcabs .gt. 1 ) then

            call uselrc                                           &
            !==========
   (idebia , idebra ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr ,                           &
    nvar   , nscal  , nphas  ,                                    &
    nideve , nrdeve , nituse , nrtuse ,                           &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                  &
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    idevel , ituser , ia     ,                                    &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,&
    dt     , rtpa   , rtp    , propce , propfa , propfb ,         &
    coefa  , coefb  , viscf  , viscb  ,                           &
    w1     , w2     , w3     , w4     , w5     ,                  &
    w6     , w7     , w8     , w9     ,                           &
    rdevel , rtuser , ra     )

          endif

        endif

      endif

    endif


! ---> Fin de la Boucle sur les scalaires physiques particulieres.
  enddo
endif


!     On calcule ici A, B, jxB

if ( ippmod(ielarc).ge.1       ) then

!     On utilise le  fait que les scalaires sont dans l'ordre
!       H, PotR, [PotI], [A] pour faire le calcul de A, B, jxB
!       apres la determination et le recalage de j
  iappel = 2

  call elflux                                                     &
  !==========
 ( idebia , idebra , iappel ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

endif

!===============================================================================
! 3. TRAITEMENT DES SCALAIRES UTILISATEURS STANDARD
!     On voit que meme s'ils sont numerotes en premier, on peut les
!       traiter en dernier si on veut.
!     On peut imaginer aussi de by-passer cette phase si le modele a
!       physique particuliere le demande.
!===============================================================================

if(nscaus.gt.0) then

! ---> Boucle sur les scalaires utilisateur.

  do ii = 1, nscaus

    iscal = ii
    ivar  = isca(iscal)

! ---> Pas de temps (avec facteur multiplicatif eventuel)

    if(cdtvar(ivar).ne.1.d0) then
      do iel = 1, ncel
        dtr(iel) = dt(iel)*cdtvar(ivar)
      enddo
    else
      do iel = 1, ncel
        dtr(iel) = dt(iel)
      enddo
    endif


! ---> Variances et scalaires

!   iscavr dira scalaire/variance
!   itspdv dira si calcul des termes de prod et dissip supplementaires
!         ou non (1 : oui, 0 : non)

!         si iscavr = 0
!           scalaire
!           itspdv = 0
!         sinon
!           variance
!           si iscavr > 0 et iscavr < nscal+1
!             itspdv = 1
!           sinon
!             pour le moment, on s'arrete
!             a terme, les combustionnistes pourront donner leur propre
!               grandeur scalaire associee a la variance et
!               eventuellement reconstruite en dehors de covofi.
!               il faudra alors peut etre declarer des tableaux
!               de travail suppl
!           fin si
!         fin si

    iisc = iscal
    if(iscavr(iisc).eq.0) then
      itspdv = 0
    elseif(iscavr(iisc).gt.0.and.iscavr(iisc).le.nscal) then
      itspdv = 1
    else
      write(nfecra,9000)iisc,iisc,nscal,iscavr(iisc)
      call csexit (1)
    endif


! ---> Appel a covofi pour la resolution

    iphas  = iphsca(iisc)

    call covofi                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdc(iphas) , ncetsm(iphas) ,                                &
   nideve , nrdeve , nituse , nrtuse , iisc   , itspdv ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) , ia(iicesm(iphas)) , ia(iitpsm(iphas)) ,    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dtr    , rtp    , rtpa   , propce , propfa , propfb , tslagr , &
   coefa  , coefb  , ra(ickupd(iphas)) , ra(ismace(iphas)) ,      &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   drtp   , smbrs  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )






! ---> Fin de la Boucle sur les scalaires utilisateurs.
  enddo

endif

!===============================================================================
! 4.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA RESOLUTION DES SCALAIRES         ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE NUMERO ',I10                                    ,/,&
'@    ISCAVR(',I10   ,') DOIT ETRE UN ENTIER                  ',/,&
'@      POSITIF OU NUL ET                                     ',/,&
'@      INFERIEUR OU EGAL A NSCAL = ',I10                      ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Si ISCAVR(I) est nul, le scalaire I n est pas une variance',/,&
'@  Si ISCAVR(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J ',/,&
'@    dont le numero est ISCAVR(I)                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE SOLVING SCALARS EQUATIONS          ',/,&
'@    ========                                                ',/,&
'@    SCALAR NUMBER ',I10                                      ,/,&
'@    ISCAVR(',I10   ,') MUST BE A POSITIVE OR NULL INTEGER   ',/,&
'@      AND LOWER OR EQUAL THAN NSCAL = ', I10                 ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculaton will not be run.                           ',/,&
'@                                                            ',/,&
'@  If ISCAVR(I) is null, the scalar I is not a variance      ',/,&
'@  If ISCAVR(I) is positive, the scalar I is a variance:     ',/,&
'@    it is the variance of the fluctuations of the scalar J  ',/,&
'@    whose number is ISCAVR(I)                               ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@  Contact support.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----
return

end subroutine

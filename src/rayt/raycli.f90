!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine raycli &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , isvhb  , isvtb  ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb ,                                     &
   izfrad , isothm ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  , hbord  , tbord  ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rayexp , rayimp ,                                              &
   tparoi , qincid , xlam   , epa    , eps    ,                   &
   flunet , flconv , hfconv , text   , tint   , tempk  ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  1) Calcul des temperatures de paroi
!  2) Mise a jours des conditions aux limites de la variable
!     energetique

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
! isvhb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  coefficients d'echange aux bords              !
! isvtb            ! e  ! <-- ! indicateur de sauvegarde des                   !
!                  !    !     !  temperatures aux bords                        !
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
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! izfrad(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
!   ,nphast)       !    !     !                                                !
! isothm(nfabor    ! te ! <-- ! type de condition de paroi                     !
!   ,nphast)       !    !     !                                                !
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
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! hbord            ! tr ! --> ! coefficients d'echange aux bords               !
! (nfabor)         !    !     !                                                !
! tbord            ! tr ! --> ! temperature aux bords           i              !
! (nfabor)         !    !     !                                                !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! rayexp(ncelet    ! tr ! --> ! terme source radiatif explicite                !
!   ,nphast)       !    !     !                                                !
! rayimp(ncelet    ! tr ! --> ! terme source radiatif implicite                !
!   ,nphast)       !    !     !                                                !
! tempk(ncelet)    ! tr ! --> ! temperature en kelvin                          !
!   ,nphast)       !    !     !                                                !
! tparoi(nfabor    ! tr ! --- ! temperature de paroi en kelvin                 !
!   ,nphast)       !    !     !                                                !
! qincid(nfabor    ! tr ! --> ! densite de flux radiatif aux bords             !
!   ,nphast)       !    !     !                                                !
! xlam (nfabor     ! tr ! --> ! coefficient de conductivite thermique          !
!   ,nphast)       !    !     ! des facettes de paroi (w/m/k)                  !
! epa (nfabor      ! tr ! --> ! epaisseur des facettes de paroi (m)            !
!   ,nphast)       !    !     !                                                !
! eps (nfabor      ! tr ! --> ! emissivite des facettes de bord                !
!   ,nphast)       !    !     !                                                !
! flunet(nfabor    ! tr ! --> ! densite de flux net radiatif aux               !
!   ,nphast)       !    !     ! faces de bord                                  !
! flconv(nfabor    ! tr ! --> ! densite de flux convectif aux faces            !
!   ,nphast)       !    !     ! de bord                                        !
! hfconv(nfabor    ! tr ! --> ! coefficient d'echange fluide aux               !
!   ,nphast)       !    !     ! faces de bord                                  !
! text (nfabor     ! tr ! --> ! temperature de bord externe                    !
!   ,nphast)       !    !     ! en degres celcius                              !
! tint (nfabor     ! tr ! --> ! temperature de bord interne                    !
!   ,nphast)       !    !     ! en degres celcius                              !
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

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "entsor.h"
include "parall.h"
include "radiat.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          isvhb  , isvtb

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          idevel(nideve), ituser(nituse), ia(*)
integer          izfrad(nfabor,nphast),isothm(nfabor,nphast)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision hbord(nfabor),tbord(nfabor)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)

double precision rayexp(ncelet,nphast), rayimp(ncelet,nphast)
double precision tempk(ncelet,nphast)
double precision tparoi(nfabor,nphast), qincid(nfabor,nphast)
double precision xlam(nfabor,nphast), epa(nfabor,nphast)
double precision eps(nfabor,nphast), flunet(nfabor,nphast)
double precision flconv(nfabor,nphast), hfconv(nfabor,nphast)
double precision text(nfabor,nphast), tint(nfabor,nphast)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac, iel, iphas, iph, ideb, ivart, iscat
integer          mode, iok, ifvu, ii, izonem, izone
integer          maxelt, idbia1, ils

double precision tmin , tmax   , tx
double precision cpp, xmtk

integer    ipacli
data       ipacli /0/
save       ipacli

!===============================================================================
!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

maxelt = max(ncelet,nfac,nfabor)

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

!---> NUMERO DE PASSAGE RELATIF

ipacli = ipacli + 1
ideb = 0

!---> VALEURS MIN ET MAX ADMISSIBLES POUR LA TEMPERATURE DE PAROI
!         EN KELVIN

tmin = 0.d0
tmax = grand + tkelvi

!---> COEFF DE RELAX

!      TX est strictement superieur a 0 et inferieur ou egal a 1

!      Pour calculer la temperature de paroi, on calcule un increment
!      de temperature DeltaT entre l'etape courante n et l'etape
!      precedente n-1, puis on calcule :
!           n    n-1                                 n-1
!          T  = T    + DeltaT si le rapport DeltaT/T    =< TX, sinon

!           n    n-1                      n-1             n-1
!          T  = T    * (1 + TX *((DeltaT/T   ) / |DeltaT/T   |))

tx = 0.1d0


!---> INITIALISATIONS PAR DEFAUT BIDON

do iph = 1,nphast
  do ifac = 1,nfabor
    izfrad(ifac,iph) = -1
    isothm(ifac,iph) = -1
    xlam  (ifac,iph) = -grand
    epa   (ifac,iph) = -grand
    eps   (ifac,iph) = -grand
    text  (ifac,iph) = -grand
    tint  (ifac,iph) = -grand
  enddo
enddo


!===============================================================================
! 2. SI PAS DE FICHIER SUITE ALORS INITIALISATION AU PREMIER PASSAGE
!    DE TPAROI ET QINCID :
!      LECTURE DE L'INITIALISATION DE TPAROI A TINT
!      QINCID EST INITIALISE A STEPHN*TINT**4 (SI ON INITIALISE QINCID
!      A ZERO, ON AURA UN DEFICIT SUR LA CONDITION LIMITE DE LUMINANCE
!      AUX PAROIS AU 1er PAS DE TEMPS EN DOM)
!===============================================================================

if (ipacli.eq.1 .and. isuird.eq.0) then

! Indicateur : si non suite et premier pas de temps.
    ideb = 1

    do iph = 1, nphast

      iphas = irapha(iph)

      do iel = 1,ncelet
        rayimp(iel,iph) = zero
        rayexp(iel,iph) = zero
      enddo

      do ifac = 1,nfabor
        hfconv(ifac,iph) = zero
        flconv(ifac,iph) = zero
      enddo

!     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
!       pour être sur que TPAROI ne sera pas modifié
!       (puisqu'on a TBORD libre)
!     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
!       pour être sur que QINCID ne sera pas modifié
!       (puisqu'on a FLUNET libre)

      do ifac = 1,nfabor
        tbord(ifac)      = zero
        flunet(ifac,iph) = zero
      enddo

!         - Interface Code_Saturne
!           ======================

      if (iihmpr.eq.1) then

!---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
        ivart = isca(iscalt(iphas))

        call uiray2                                               &
        !==========
       ( itypfb, iparoi, iparug, ivart, iph, nphast,  izfrad,     &
         isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,          &
         nfabor, nfml, ifmfbr , iprfml , nvar,                    &
         eps, epa, tint, text,                                    &
         xlam, rcodcl)

      endif

      ils    = idebia
      idbia1 = ils + maxelt
      CALL IASIZE('RAYCLI',IDBIA1)

      call usray2                                                 &
      !==========
 ( idbia1 , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   maxelt , ia(ils),                                              &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , izfrad(1,iph) , isothm(1,iph) ,                       &
   tmin   , tmax   , tx     ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &
   tbord  , flunet(1,iph) , hfconv(1,iph) ,  flconv(1,iph) ,      &
   xlam(1,iph)   , epa(1,iph)    , eps(1,iph)    ,                &
   text(1,iph)   , tint(1,iph)   ,                                &
   ra     )

      write(nfecra,1000)

! Tparoi en Kelvin et QINCID en W/m2
      do ifac = 1,nfabor
        if ( itypfb(ifac,iphas).eq.iparoi .or.                    &
             itypfb(ifac,iphas).eq.iparug ) then
          tparoi(ifac,iph) = tint(ifac,iph)
          qincid(ifac,iph) = stephn*tint(ifac,iph)**4
        else
          tparoi(ifac,iph) = 0.d0
          qincid(ifac,iph) = 0.d0
        endif
      enddo

! Fin de la boucle sur NPHAST
    enddo

! Fin détection premier passage
endif

!===============================================================================
! 3. BOUCLE SUR LES PHASES...
!===============================================================================

do iph = 1,nphast

  iphas = irapha(iph)

!---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
  iscat = iscalt(iphas)
  ivart = isca(iscalt(iphas))

!===============================================================================
! 3.1 DONNEES SUR LES FACES FRONTIERES
!===============================================================================

!     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
!       pour être sur que TPAROI ne sera pas modifié
!       (puisqu'on a TBORD libre)
!     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
!       pour être sur que QINCID ne sera pas modifié
!       (puisqu'on a FLUNET libre)

  do ifac = 1,nfabor
    tbord (ifac)     = tparoi(ifac,iph)
    flunet(ifac,iph) = qincid(ifac,iph)
  enddo

!     - Interface Code_Saturne
!       ======================

  if (iihmpr.eq.1) then

    call uiray2                                                   &
    !==========
  ( itypfb, iparoi, iparug, ivart, iph, nphast, izfrad,           &
    isothm, itpimp, ipgrno, iprefl, ifgrno, ifrefl,               &
    nfabor, nfml, ifmfbr , iprfml , nvar,                         &
    eps, epa, tint, text,                                         &
    xlam, rcodcl)

  endif

  ils    = idebia
  idbia1 = ils + maxelt
  CALL IASIZE('RAYCLI',IDBIA1)

  call usray2                                                     &
  !==========
 ( idbia1 , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   maxelt , ia(ils),                                              &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , izfrad(1,iph) , isothm(1,iph) ,                       &
   tmin   , tmax   , tx     ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &
   tbord  , flunet(1,iph) ,  hfconv(1,iph) ,  flconv(1,iph) ,     &
   xlam(1,iph)   , epa(1,iph)    , eps(1,iph)    ,                &
   text(1,iph)   , tint(1,iph)   ,                                &
   ra     )

!===============================================================================
! 3.2 CONTROLE DES DONNEES UTILISATEUR
!===============================================================================

!--> Arret si le numero de zone est non renseigne ou mal renseigne

  iok = 0

  do ifac = 1, nfabor
    if (izfrad(ifac,iph).le.0.or.izfrad(ifac,iph).gt.nozrdm) then
      iok = iok + 1
      write(nfecra,2000)ifac,nozrdm,izfrad(ifac,iph)
    endif
  enddo

  if(iok.ne.0) then
    call csexit (1)
    !==========
  endif

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
!     Stop si depassement.

  nzfrad(iphas) = 0
  do ifac = 1, nfabor
    ifvu = 0
    do ii = 1, nzfrad(iphas)
      if (ilzrad(ii,iphas).eq.izfrad(ifac,iph)) then
        ifvu = 1
      endif
    enddo
    if(ifvu.eq.0) then
      nzfrad(iphas) = nzfrad(iphas) + 1
      if(nzfrad(iphas).le.nbzrdm) then
        ilzrad(nzfrad(iphas),iphas) = izfrad(ifac,iph)
      else
        write(nfecra,2001) nbzrdm
        write(nfecra,2002)(ilzrad(ii,iphas),ii=1,nbzrdm)
        call csexit (1)
        !==========
      endif
    endif
  enddo

! ---> Plus grand numero de zone atteint

  izonem = 0
  do ii = 1, nzfrad(iphas)
    izone = ilzrad(ii,iphas)
    izonem = max(izonem,izone)
  enddo
  if(irangp.ge.0) then
    call parcmx(izonem)
    !==========
  endif
  nozarm(iphas) = izonem




! On verra si ca coute cher ou non.
!   Pour le moment on le fait tout le temps.
!        IF(IWARNI(IVART).GE.-1.OR.IPACLI.LE.3) THEN
  if(1.eq.1) then

    iok = 0

!--> Si en paroi ISOTHM non renseignee : stop
    do ifac = 1, nfabor
      if( (itypfb(ifac,iphas).eq.iparoi  .or.                     &
           itypfb(ifac,iphas).eq.iparug) .and.                    &
           isothm(ifac,iph)  .eq.-1    ) then
        iok = iok + 1
        write(nfecra,2110) iphas,ifac,izfrad(ifac,iph)
      endif
    enddo

!--> Si ISOTHM renseignee en non paroi : stop
    do ifac = 1, nfabor
      if( itypfb(ifac,iphas).ne.iparoi .and.                      &
          itypfb(ifac,iphas).ne.iparug .and.                      &
          isothm(ifac,iph)  .ne.-1         ) then
        iok = iok + 1
        write(nfecra,2111)                                        &
             iphas,ifac,izfrad(ifac,iph),isothm(ifac,iph)
      endif
    enddo

!--> Si valeur physique erronee : stop
    do ifac = 1, nfabor
      if(isothm(ifac,iph).eq.itpimp ) then
        if(eps(ifac,iph) .lt.0.d0.or.eps(ifac,iph).gt.1.d0.or.    &
           tint(ifac,iph).le.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2120) iphas,ifac,izfrad(ifac,iph),         &
               eps(ifac,iph),                                     &
                              tint(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.ipgrno ) then
        if(eps(ifac,iph) .lt.0.d0.or.eps(ifac,iph).gt.1.d0.or.    &
           xlam(ifac,iph).le.0.d0.or.                             &
           epa (ifac,iph).le.0.d0.or.                             &
           text(ifac,iph).le.0.d0.or.                             &
           tint(ifac,iph).le.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2130) iphas,ifac,izfrad(ifac,iph),         &
               eps(ifac,iph),xlam(ifac,iph),epa(ifac,iph),        &
               text(ifac,iph),tint(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.iprefl ) then
        if(xlam(ifac,iph).le.0.d0.or.                             &
           epa (ifac,iph).le.0.d0.or.                             &
           text(ifac,iph).le.0.d0.or.                             &
           tint(ifac,iph).le.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2140) iphas,ifac,izfrad(ifac,iph),         &
                             xlam(ifac,iph),epa(ifac,iph),        &
               text(ifac,iph),tint(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.ifgrno ) then
        if(eps(ifac,iph) .lt.0.d0.or.eps(ifac,iph).gt.1.d0.or.    &
           tint(ifac,iph).le.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2150) iphas,ifac,izfrad(ifac,iph),         &
               eps(ifac,iph),                                     &
                              tint(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.ifrefl ) then
        if(tint(ifac,iph).le.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2160) iphas,ifac,izfrad(ifac,iph),         &
                              tint(ifac,iph)
        endif
      elseif(isothm(ifac,iph).ne.-1) then
          iok = iok + 1
          write(nfecra,2170) iphas,ifac,izfrad(ifac,iph),         &
                             isothm(ifac,iph)
      endif
    enddo

!--> Si valeur renseignee sans raison : stop
    do ifac = 1, nfabor
     if(isothm(ifac,iph).eq.itpimp ) then
        if(xlam(ifac,iph).gt.0.d0.or.epa(ifac,iph).gt.0.d0.or.    &
           text(ifac,iph).gt.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2220) iphas,ifac,izfrad(ifac,iph),         &
               xlam(ifac,iph),epa(ifac,iph),text(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.iprefl ) then
        if(eps(ifac,iph).ge.0.d0                           ) then
          iok = iok + 1
          write(nfecra,2240) iphas,ifac,izfrad(ifac,iph),         &
               eps(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.ifgrno ) then
        if(xlam(ifac,iph).gt.0.d0.or.epa(ifac,iph).gt.0.d0.or.    &
           text(ifac,iph).gt.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2250) iphas,ifac,izfrad(ifac,iph),         &
               xlam(ifac,iph),epa(ifac,iph),text(ifac,iph)
        endif
      elseif(isothm(ifac,iph).eq.ifrefl ) then
        if(eps(ifac,iph).ge.0.d0.or.                              &
           xlam(ifac,iph).gt.0.d0.or.epa(ifac,iph).gt.0.d0.or.    &
           text(ifac,iph).gt.0.d0                          ) then
          iok = iok + 1
          write(nfecra,2260) iphas,ifac,izfrad(ifac,iph),         &
               eps(ifac,iph),                                     &
               xlam(ifac,iph),epa(ifac,iph),text(ifac,iph)
        endif
      endif
    enddo

!--> Stop si erreur
    if(iok.ne.0) then
      call csexit (1)
      !==========
    endif

  endif


!===============================================================================
! 3.2 COMPLETION DES DONNEES UTILISATEUR
!===============================================================================

! ICODCL et EPS (quand il est nul)

  do ifac = 1, nfabor
    if(    isothm(ifac,iph).eq.itpimp ) then
      icodcl(ifac,ivart) = 5
    elseif(isothm(ifac,iph).eq.ipgrno ) then
      icodcl(ifac,ivart) = 5
    elseif(isothm(ifac,iph).eq.iprefl ) then
      icodcl(ifac,ivart) = 5
      eps(ifac,iph) = 0.d0
    elseif(isothm(ifac,iph).eq.ifgrno ) then
      icodcl(ifac,ivart) = 5
    elseif(isothm(ifac,iph).eq.ifrefl ) then
      icodcl(ifac,ivart) = 3
      eps(ifac,iph) = 0.d0
    endif
  enddo


!===============================================================================
! 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL,IPH)
!===============================================================================

  if (abs(iscsth(iscat)).eq.1) then

!---> ON REMPLIT TEMPK

    if (iscsth(iscat).eq.-1) then
      do iel = 1, ncel
        tempk(iel,iph) = rtpa(iel,ivart) + tkelvi
      enddo
    else
      do iel = 1, ncel
        tempk(iel,iph) = rtpa(iel,ivart)
      enddo
    endif

  else if (iscsth(iscat).eq.2) then

!---> LECTURES DES DONNEES UTILISATEURS (TBORD est un auxiliaire)

    mode = 1

    if (ippmod(iphpar).le.1) then

      call usray4                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &

   tparoi(1,iph)   , tbord  , tempk(1,iph)   ,                    &
!                                   Resultat : T en K

   ra     )

    else

      call ppray4                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &

   tparoi(1,iph)   , tbord  , tempk(1,iph)   ,                    &
!                                   Resultat : T en K

   ra     )

    endif

  endif

!===============================================================================
! 5.  CALCUL DES TEMPERATURES DE PAROIS
!===============================================================================


! DANS TOUS LES CAS HFCONV CONTIENT Lambda * Hturb / distance
!   (HFCONV : W/(m2 K) ; Hturb est sans dimension)
!  (au premier passage, il est nul)

!--> CALCUL DU FLUX CONVECTIF
!      Par flux convectif, on entend bien sur
!        flux convectif parallele a la paroi,
!        on suppose que la paroi est etanche...
!      Le flux est calcule dans condli clptur, sauf au premier
!        passage sans suite de calcul, puisque raycli est appele avant.


if (ideb.eq.1) then

  do ifac = 1,nfabor
    if (isothm(ifac,iph).ne.-1) then
      flconv(ifac,iph) =                                          &
      hfconv(ifac,iph)*(tempk(ifabor(ifac),iph)-tparoi(ifac,iph))
    endif
  enddo

endif


!--> Les cas ou il faut calculer TPAROI sont, au premier passage sans suite
!      des cas a temperature imposee TPAROI = TINT

  if (ideb.eq.1) then

    do ifac = 1,nfabor
      if (isothm(ifac,iph).eq.ipgrno .or.                         &
          isothm(ifac,iph).eq.iprefl .or.                         &
          isothm(ifac,iph).eq.ifgrno    ) then
        isothm(ifac,iph) = itpimp
      endif
    enddo

  endif

  if(ideb.eq.0) then

    call raypar                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &

   icodcl , isothm(1,iph)  , izfrad(1,iph) ,                      &

   idevel , ituser , ia     ,                                     &

   tmin   , tmax   , tx     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   rdevel , rtuser ,                                              &

   tparoi(1,iph) , qincid(1,iph) , text(1,iph)   , tint(1,iph)   ,&
   xlam(1,iph)   , epa(1,iph)    , eps(1,iph)    ,                &
   hfconv(1,iph) , flconv(1,iph) , tempk(1,iph)  ,                &

   ra     )

  endif


!===============================================================================
! 6.  CHANGEMENT DES CONDITIONS LIMITES UTILISATEUR
!===============================================================================

!===============================================================================
! 6.1  LA VARIABLE TRANSPORTEE EST LA TEMPERATURE
!===============================================================================

  if (abs(iscsth(iscat)).eq.1) then

    if(iscsth(iscat).eq.-1) then
      xmtk = -tkelvi
    else
      xmtk = 0.d0
    endif

    do ifac = 1,nfabor

      if (isothm(ifac,iph).eq.itpimp .or.                         &
          isothm(ifac,iph).eq.ipgrno .or.                         &
          isothm(ifac,iph).eq.ifgrno    ) then
        rcodcl(ifac,ivart,1) = tparoi(ifac,iph)+xmtk
        rcodcl(ifac,ivart,2) = rinfin
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac,iph).eq.iprefl) then
        rcodcl(ifac,ivart,1) = text(ifac,iph)+xmtk
        rcodcl(ifac,ivart,2) = xlam(ifac,iph)/epa(ifac,iph)
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac,iph).eq.ifrefl) then
        icodcl(ifac,ivart) = 3
        rcodcl(ifac,ivart,1) = 0.d0
        rcodcl(ifac,ivart,2) = rinfin
      endif

    enddo

!===============================================================================
! 6.2  LA VARIABLE TRANSPORTEE EST L'ENTHALPIE
!===============================================================================

  elseif (iscsth(iscat).eq.2) then

!---> LECTURES DES DONNEES UTILISATEURS
!     ON CONVERTIT TPAROI EN ENTHALPIE DE BORD, STOCKEE DANS FLUNET,
!     QUI EST UTILISE COMME AUXILIAIRE

    mode = 0

    do ifac = 1,nfabor
      if (isothm(ifac,iph).eq.itpimp.or.                          &
          isothm(ifac,iph).eq.ipgrno.or.                          &
          isothm(ifac,iph).eq.ifgrno    ) then
        mode = -1
      endif
    enddo

    if (mode.eq.-1) then

      if (ippmod(iphpar).le.1) then

        call usray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &

   tparoi(1,iph)   , flunet(1,iph) , tempk(1,iph)  ,              &
!                          HPAROI
   ra     )

      else

        call ppray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &

   tparoi(1,iph)   , flunet(1,iph) , tempk(1,iph)  ,              &
!                          HPAROI
   ra     )

      endif

    endif

    mode = 0

    do ifac = 1,nfabor
      if (isothm(ifac,iph).eq.iprefl) then
        mode = -1
      endif
    enddo

    if (mode.eq.-1) then

      if (ippmod(iphpar).le.1) then

        call usray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &

   text(1,iph)  , tbord  , tempk(1,iph)  ,                        &
!                       HEXT
   ra     )

      else

        call ppray4                                               &
        !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &

   mode   ,                                                       &

   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb(1,iphas) , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &

   text(1,iph)  , tbord  , tempk(1,iph)  ,                        &
!                       HEXT
   ra     )

      endif

    endif

    do ifac = 1,nfabor

      if (isothm(ifac,iph).eq.itpimp.or.                          &
          isothm(ifac,iph).eq.ipgrno.or.                          &
          isothm(ifac,iph).eq.ifgrno    ) then
        rcodcl(ifac,ivart,1) = flunet(ifac,iph)
        rcodcl(ifac,ivart,2) = rinfin
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac,iph).eq.iprefl) then

        if (icp(iphas).gt.0) then
          iel = ifabor(ifac)
          cpp = propce(iel,ipproc(icp(iphas)))
        else
          cpp = cp0(iphas)
        endif

        rcodcl(ifac,ivart,1) = tbord(ifac)
        rcodcl(ifac,ivart,2) =                                    &
             xlam(ifac,iph)/(epa(ifac,iph)*cpp)
        rcodcl(ifac,ivart,3) = 0.d0

      else if (isothm(ifac,iph).eq.ifrefl) then
        icodcl(ifac,ivart) = 3
        rcodcl(ifac,ivart,1) = 0.d0
        rcodcl(ifac,ivart,2) = rinfin
      endif

    enddo

  endif

!===============================================================================
! 7.  FIN DE LA BOUCLE SUR LES PHASES
!===============================================================================

enddo

!--------
! FORMATS
!--------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT ',/,  &
           3X,'   ------------------------------------------',/,  &
           3X,' Initialisation de la temperature de paroi   ',/,  &
           3X,' (TPAROI) avec le profil utilisateur (TINTP) ',/,  &
           3X,' et du flux incident aux parois (QINCID).    ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,&
'@                                                            ',/,&
'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,&
'@    un entier strictement positif et inferieur ou egal a    ',/,&
'@    NOZRDM = ',I10                                           ,/,&
'@  Ce numero (IZFRDP(IFAC)) vaut ici ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans usray2.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RAYONNEMENT                                 ',/,&
'@    =========                                               ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NBZRDM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites dans usray2.          ',/,&
'@                                                            ',/,&
'@  Les NBZRDM premieres zones frontieres                     ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2002 format(i10)

 2110 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP DOIT ETRE RENSEIGNE SUR TOUTES LES FACES DE PAROI',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Il ne l''a pas ete pour la face ',I10                      ,/,&
'@                    zone         ',I10                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2111 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ISOTHP A ETE RENSEIGNE SUR UNE FACE NON PAROI           ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Sur la face ',I10   ,', zone  ',I10   ,', ISOTHP a ete    ',/,&
'@    renseigne dans usray2 (ISOTHP = ',I10   ,') alors que   ',/,&
'@    la face n''a pas ete declaree de type IPAROI ou IPARUG  ',/,&
'@    dans usclim.                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2 et usclim.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2120 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,&
'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ITPIMP       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2130 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP DOIT ETRE UN REEL INCLUS DANS [0.; 1.]             ',/,&
'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,&
'@                                      STRICTEMENT POSITIFS  ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPGRNO       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5    ,' TINTP = ',E14.5                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2140 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,&
'@                                      STRICTEMENT POSITIFS  ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5    ,' TINTP = ',E14.5                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2150 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,&
'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFGRNO       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2160 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@  TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF               ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,&
'@  TINTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2170 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@   VALEUR NON ADMISSIBLE DE ISOTHP                          ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ',I10         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2220 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP ET TEXTP NE DOIVENT PAS ETRE RENSEIGNES     ',/,&
'@                                     AVEC ISOTHP = ITPIMP   ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ITPIMP       ',/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2240 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    EPSP NE DOIT PAS ETRE RENSEIGNE AVEC ISOTHP = IPREFL    ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,&
'@                                       AVEC ISOTHP = IFGRNO ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFGRNO       ',/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2260 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,&
'@                                       AVEC ISOTHP = IFREFL ',/,&
'@                                                            ',/,&
'@  Phase ',I10                                                ,/,&
'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,&
'@  EPSP  = ',E14.5                                            ,/,&
'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,&
'@  TEXTP = ',E14.5                                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end

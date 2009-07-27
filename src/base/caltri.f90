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

subroutine caltri &
!================

 ( iverif ,                                                       &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,          &
   volume ,                                                       &
   rdevel , rtuser , ra     )

!===============================================================================
!  FONCTION  :
!  ----------

! GESTION DU PROGRAMME (LECTURE, RESOLUTION, ECRITURE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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
! (nfml,nprfml)    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfabor+1)     !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr)       !    !     !  (optionnel)                                   !
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
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! (ndim,nnod)      !    !     !                                                !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________.____._____.________________________________________________.

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
include "dimens.h"
include "pointe.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "albase.h"
include "period.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "coincl.h"
include "cpincl.h"
include "lagpar.h"
include "lagdim.h"
include "lagran.h"
include "vortex.h"
include "ihmpre.h"
include "matiss.h"
include "radiat.h"
include "cplsat.h"

!===============================================================================

! Arguments

integer          nideve , nrdeve , nituse , nrtuse

integer          iverif
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          ipropc , ipropf , ipropb
integer          icoefa , icoefb
integer          irtp   , irtpa
integer          idt
integer          iisstd , ifrcx

integer          idebia , idebra
integer          ifinia , ifinra , idbia1 , idbra1, idbia2
integer          iditia, iditra
integer          ifnia1 , ifnra1 , ifnia2 , ifnia3, ifnra2
integer          jcelbr
integer          iiii

integer          modhis, iappel, modntl, iisuit, iwarn0
integer          ntsdef, nthdef, ntcrel, ntcam1
integer          iphas , ivar

integer          iicoce , iityce
integer          iiitep , iitepa , istatc , istatf
integer          iettp  , iettpa , iauxl  , itslag , istatv
integer          itaup  , iitlag , ipiil  , iindep , iibord
integer          ivagau , itsuf  , itsup  , ibx    , iauxl2
integer          ibrgau , itebru
integer          igradp , igradv , icroul
integer          itepct , itsfex , itsvar
integer          icpgd1 , icpgd2 , icpght
integer          ilagia , ilagra , iiwork
integer          iw1    , iw2    , iw3
integer          inod   , idim
integer          itrale , indact , indwri
integer          maxelt , ils

double precision titer1, titer2
double precision tecrf1, tecrf2


!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

!     Initialisation de la premiere case libre des tableaux IA et RA
idebia = 1
idebra = 1

!     Initialisation du generateur de nombres aleatoires
!     (pas toujours necessaire mais ne coute rien)
call zufalli(0)
!===========

!---> MISE A ZERO DES TABLEAUX UTILISATEUR ET DEVELOPPEUR

if(nituse.gt.0) then
  do iiii = 1, nituse
    ituser(iiii) = 0
  enddo
endif
if(nrtuse.gt.0) then
  do iiii = 1, nrtuse
    rtuser(iiii) = 0.d0
  enddo
endif

if(nideve.gt.0) then
  do iiii = 1, nideve
    idevel(iiii) = 0
  enddo
endif
if(nrdeve.gt.0) then
  do iiii = 1, nrdeve
    rdevel(iiii) = 0.d0
  enddo
endif

!---> Test d'arret mis a 1 si le rayonnement P-1 voit trop de cellules
!       a epaisseur optique superieure a l'unite (voir ppcabs)
istpp1 = 0

!---> Nombre max d'elements pour le selector
maxelt = max(ncelet,nfac,nfabor)

!===============================================================================
! 2. GEOMETRIE
!===============================================================================

!---> CALCULS GEOMETRIQUES

!     (MEMCLG remplit directement pointe.h)
call memclg                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifinia , ifinra )

!---> La memoire sera conservee jusqu'a la fin.
idebia = ifinia
idebra = ifinra


call cregeo                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   rdevel , rtuser , ra     )



!---> CALCUL DE JCELBR (=NCELBR) et REMPLISSAGE DE IA(IICELB)
!         directement en passant par le pointeur et IA.
!         (c'est pas bien, mais c'est mieux qu'avant)

iicelb = idebia
call memcbr                                                       &
!==========
 ( iicelb , ncelet , ncel   , nfabor ,                            &
   jcelbr , ifinia ,                                              &
   ifabor ,                                                       &
   ia     )

!---> La memoire sera conservee jusqu'a la fin.
idebia = ifinia
idebra = ifinra

!===============================================================================
! 3. FIN INITIALISATION DES COMMONS
!===============================================================================

call initi2                                                       &
!==========
 ( idebia , idebra ,                                              &
   jcelbr ,                                                       &
   ia     , ra     )

if (iilagr.gt.0) then

!--> Calcul de LNDNOD (lagran.h)

!   Tableau NCELET de travail entier
  iiwork = idebia
  ifinia = iiwork + ncelet
  call iasize ('lagini',ifinia)

  ifinra = idebra

  call lagini                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ncelet , ncel   , nfac   , nfabor ,                            &
   lndnod ,                                                       &
   ifacel , ifabor ,                                              &
   ia(iiwork) ,                                                   &
   ia     , ra     )

endif


!===============================================================================
! 4. AUTRES TABLEAUX
!===============================================================================

!---> GESTION MEMOIRE

call memtri                                                       &
!==========
 ( idebia , idebra , iverif ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncofab , nproce , nprofa , nprofb ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iisstd , ifrcx  ,                                              &
   idt    , irtp   , irtpa  , ipropc , ipropf , ipropb ,          &
   icoefa , icoefb ,                                              &
   ifinia , ifinra )

!     Reservations complementaires pour Matisse
!       On remplit un pointeur dans matiss.h, pour eviter de charger
!       les arguments par des tableaux specifiques a Matisse.

if (imatis.eq.1) then

  idbia1 = ifinia
  idbra1 = ifinra
  call memmat                                                     &
  !==========
 ( idbia1 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncofab , nproce , nprofa , nprofb ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifinia , ifinra )

endif

!     Reservations memoire pour les tableaux complémentaires
!         nécesaires pour les physiques particulieres

idbia1 = ifinia
idbra1 = ifinra

call memppt                                                       &
!==========
 ( idbia1 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   nvar   , nscal  , nphas  ,                                     &
   ifinia , ifinra )


!===============================================================================
! 4.1 RESERVATION DE LA MEMOIRE POUR LE RAYONNEMENT SEMI-TRANSPARENT
!===============================================================================

if (iirayo.gt.0) then

  idbia1 = ifinia
  idbra1 = ifinra

  call memra1                                                     &
  !==========
 ( idbia1 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   nvar   , nscal  , nphas  ,                                     &
   ifinia , ifinra )

endif

!===============================================================================
! 4.2 RESERVATION DE LA MEMOIRE POUR LE LAGRANGIEN
!===============================================================================

!     Si on ne fait pas de Lagrangien, on initialise
!       quand meme les "pointeurs".

idbia1 = ifinia
idbra1 = ifinra

call memla1                                                       &
!==========
  ( idbia1 , idbra1 ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor ,                  &
    lndnod ,                                                      &
    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
    ntersl , nvlsta , nvisbr ,                                    &
    iiitep , iicoce , iityce ,                                    &
    iettp  , iettpa , iitepa , istatc , istatv, itslag , istatf , &
    ifinia , ifinra )

!===============================================================================
! 4.3 TESTS ELEMENTAIRES : APPEL A TESTEL.F
!===============================================================================

if (iverif.eq.1) then

  write(nfecra, 1000)

  call testel                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  , nvar   ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume , ra(irtp) ,                                   &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )

  goto 200

endif

!===============================================================================
! 5. INITIALISATIONS PAR DEFAUT
!===============================================================================

call iniva0                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncofab ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx ) ,                         &
   rdevel , rtuser , ra     )

!===============================================================================
! 6. CALCUL SUITE EVENTUEL
!===============================================================================

if (isuite.eq.1) then

  call lecamo                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nnod   ,          &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx ) ,                         &
   rdevel , rtuser , ra     )

!     En ALE, il faut recalculer les parametres geometriques
  if (iale.eq.1) then

    do inod = 1, nnod
      do idim = 1, ndim
        xyznod(idim,inod) = ra(ixyzn0+(inod-1)*ndim+idim-1)       &
             + ra(idepal+(idim-1)*nnod+inod-1)
      enddo
    enddo

    call algrma
    !==========
    call calgeo                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   volmin , volmax , voltot ,                                     &
   rdevel , rtuser , ra     )

  endif

endif

!===============================================================================
! 7. INITIALISATIONS (Utilisateur et complementaires)
!    RTP DT ROM ROMB VISCL VISCT VISCLS
!    (TPUCOU en PERIODICITE)
!===============================================================================

call inivar                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncofab ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx ) ,                         &
   rdevel , rtuser , ra     )

!===============================================================================
! 8.1 MODULE DE RAYONNEMENT : CALCUL SUITE EVENTUEL
!===============================================================================

if (iirayo.gt.0 .and. isuird.eq.1) then

  call raylec                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   ra(ipropc) , ra(ipropb) ,                                      &
   rdevel , rtuser , ra     )

endif

!===============================================================================
! 8.2 INITIALISATIONS DES PARTICULES POUR LE LAGRANGIEN
!===============================================================================

if (iilagr.gt.0) then

  call laglec                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   ia(iiitep) , ia ,                                              &
   ra(irtpa)  , ra(ipropc) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatc) , ra(istatv) ,            &
   ra(istatf) , ra(itslag) , ra         )

endif

!===============================================================================
! 8.3 INITIALISATIONS POUR LE MODULE THERMIQUE 1D EN PAROI
!===============================================================================
! On suppose que toutes les phases voient la meme temperature de paroi
! USPT1D a un fonctionnement similaire a USKPDC et USTSMA, mais comme
! on ecrit des infos dans un fichier suite, on a besoin d'une partie de
! la memoire meme apres la boucle en temps -> IFPT1D et TPPT1D
!                                            (IFNIA1 et IFNRA1)

idbia1 = ifinia
idbra1 = ifinra

ils    = idbia1
ifnia2 = ils + maxelt
call iasize('caltri',ifnia2)

iphas = 1

!     Premier appel : definition de NFPT1D et ISUIT1
iappel = 1
call uspt1d                                                       &
!==========
 ( ifnia2 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(idbia1) , ia(idbia1) , ia(idbia1) ,                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idbra1) , ra(idbra1) , ra(idbra1) ,                         &
   ra(idbra1) , ra(idbra1) , ra(idbra1) ,                         &
   ra(idbra1) , ra(idbra1) , ra(idbra1) ,                         &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )

iappel = 1
call vert1d                                                       &
!==========
 (idbia1     , idbra1     ,                                       &
  nfabor     , nfpt1d     , iappel     ,                          &
  ia(idbia1) , ia(idbia1) , ia(idbia1) , ia     ,                 &
  ra(idbra1) , ra(idbra1) ,                                       &
  ra(idbra1) , ra(idbra1) , ra(idbra1) , ra     )

call memt1d                                                       &
!==========
 ( idbia1 , idbra1 , nfabor , ifnia1 , ifnra1 ,ifnia2 , ifnra2 ,  &
   ifinia , ifinra , ia     , ra     )

! On appelle uspt1d lorqu'il y a sur un processeur au moins des faces de
!     bord avec module thermique 1D.

if (nfpt1t.gt.0) then
! Deuxieme appel : remplissage des tableaux de definition de la geometrie
!            et de l'initialisation (IFPT1D,NPPT1D,EPPT1D,RGPT1D,TPPT1D)
  ils    = ifinia
  ifnia3 = ils + maxelt
  call iasize('caltri',ifnia3)

  iappel = 2
  call  uspt1d                                                    &
  !===========
 ( ifnia3 , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iifpt1) , ia(inppt1) , ia(iiclt1) ,                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(itppt1) , ra(irgpt1) , ra(ieppt1) ,                         &
   ra(itept1) , ra(ihept1) , ra(ifept1) ,                         &
   ra(ixlmt1) , ra(ircpt1) , ra(idtpt1) ,                         &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )

  iappel = 2
  call vert1d                                                     &
  !==========
 (ifinia     , ifinra     ,                                       &
  nfabor     , nfpt1d     , iappel     ,                          &
  ia(iifpt1) , ia(inppt1) , ia(iiclt1) , ia     ,                 &
  ra(irgpt1) , ra(ieppt1) ,                                       &
  ra(ixlmt1) , ra(ircpt1) , ra(idtpt1) , ra     )

!     Calcul du max des NPPT1D (pour les fichiers suite)
  nmxt1d = 0
  do iiii = 1, nfpt1d
    nmxt1d = max(ia(inppt1+iiii-1),nmxt1d)
  enddo
  if (irangp.ge.0) call parcmx(nmxt1d)
                   !==========

  if (isuit1.eq.1) then

    call lect1d                                                   &
    !==========
 ( ficmt1     , len(ficmt1), nfpt1d     , nfpt1t    ,             &
   nmxt1d     , nfabor     , ia(inppt1) , ia(iifpt1) , ra(ieppt1),&
   ra(irgpt1) , ra(itppt1))

  else
!     Creation du maillage, initialisation de la temperature.

    call mait1d                                                   &
    !==========
 ( nfpt1d, ia(inppt1), ra(ieppt1), ra(irgpt1),ra(itppt1))

  endif

endif
!     Les infos sur l'epaisseur de la paroi, le nombre de points de
!     discretisation et la raison geometrique ont ete transmises a
!     la structure C. Elles sont maintenant inutiles dans le Fortran.
!     -> on libere la memoire.
ifinia = ifnia2
ifinra = ifnra2

!===============================================================================
! 9. TABLEAUX POUR BLC EN TEMPS MAIS A OUBLIER ENSUITE
!===============================================================================

!  En fin de bloc en temps on doit retrouver IFNIA1 et IFNRA1
iditia = ifnia1
iditra = ifnra1

idbia1 = ifinia
idbra1 = ifinra


do iphas = 1, nphas

  iappel = 1

  if (imatis.eq.1) then

!     Noter que uskpdc n'est pas permis avec Matisse
!       (si necessaire, on pourrait regrouper toutes les pertes de
!        charge de Matisse dans un unique sous-programme et reactiver
!        uskpdc)

    call  mtkpdc                                                  &
    !===========
 ( idbia1 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncepdc(iphas) , iphas  , iappel ,                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(idbia1),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(idbra1) ,                         &
   rdevel , rtuser , ra     )

  else

    ils    = idbia1
    idbia2 = ils + maxelt
    call iasize('caltri',idbia2)

    call  uskpdc                                                  &
    !===========
 ( idbia2 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncepdc(iphas) , iphas  , iappel ,                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(idbia1),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(idbra1) ,                         &
   rdevel , rtuser , ra     )

  endif

enddo

call mempdc                                                       &
!==========
 ( idbia1, idbra1, ncelet, ncel,  nphas, ndim, ifinia, ifinra)


! On appelle uskpdc lorqu'il y a sur un processeur au moins des cellules
!     avec terme source de masse.
!     On ne fait que remplir le tableau d'indirection des cellules
!     On appelle cependant uskpdc avec tous les processeurs, au cas ou
!     l'utilisateur aurait mis en oeuvre des operations globales.

do iphas = 1, nphas

  if(ncpdct(iphas).gt.0) then

    iappel = 2

    if (imatis.eq.1) then

      call  mtkpdc                                                &
      !===========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncepdc(iphas) , iphas  , iappel ,                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)),                                             &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(ickupd(iphas)) ,                  &
   rdevel , rtuser , ra     )

    else

      ils    = ifinia
      ifnia2 = ils + maxelt
      call iasize('caltri',ifnia2)

      call  uskpdc                                                &
      !===========
 ( ifnia2 , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncepdc(iphas) , iphas  , iappel ,                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)),                                             &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(ickupd(iphas)) ,                  &
   rdevel , rtuser , ra     )

    endif

  endif

enddo

idbia1 = ifinia
idbra1 = ifinra

ils    = idbia1
idbia2 = ils + maxelt
call iasize('caltri',idbia2)

do iphas = 1, nphas

  iappel = 1
  call  ustsma                                                    &
  !===========
 ( idbia2 , idbra1 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdc(iphas)   ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncetsm(iphas) ,   iphas  , iappel ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) ,                                            &
   ia(idbia1) , ia(idbia1),                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(ickupd(iphas))       , ra(idbra1),&
   rdevel , rtuser , ra     )

enddo

call memtsm                                                       &
!==========
     ( idbia1 , idbra1 ,                                          &
       ncelet , ncel   , nvar   , nphas  ,                        &
       ifinia , ifinra )

! On appelle ustsma lorqu'il y a sur un processeur au moins des cellules
!     avec terme source de masse.
!     On ne fait que remplir le tableau d'indirection des cellules
!     On appelle cependant ustsma avec tous les processeurs, au cas ou
!     l'utilisateur aurait mis en oeuvre des operations globales.

do iphas = 1, nphas

  if(nctsmt(iphas).gt.0) then

    ils    = ifinia
    ifnia2 = ils + maxelt
    call iasize('caltri',ifnia2)

    iappel = 2
    call  ustsma                                                  &
    !===========
 ( ifnia2 , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdc(iphas)   ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncetsm(iphas) ,   iphas  , iappel ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicepd(iphas)) ,                                            &
   ia(iicesm(iphas)) , ia(iitpsm(iphas)),                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(ickupd(iphas)), ra(ismace(iphas)),&
   rdevel , rtuser , ra     )

  endif

enddo


! -- Methode des vortex pour la L.E.S.
!    (dans verini on s'est deja assure que ITYTUR=4 si IVRTEX=1)

if (ivrtex.eq.1) then

  idbia1 = ifinia
  idbra1 = ifinra

  iphas  = 1
  iappel = 1

!  On met une valeur factice a certains parametres non utilise en IAPPEL=1

  call memvor                                                     &
  !==========
 ( idbia1 , idbra1 , iappel , nfabor , ifinia , ifinra )

  call vorin0( nfabor , ia(iirepv) )
  !==========

  ils    = ifinia
  ifnia2 = ils + maxelt
  call iasize('caltri',ifnia2)

  call usvort                                                     &
  !==========
 ( ifnia2 , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , iappel ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iirepv)      ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  ,                                      &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )

  call vorver ( nfabor , ia(iirepv)  , iappel )
  !==========

  idbia1 = ifinia
  idbra1 = ifinra

! Attention, vorpre reserve de la memoire qu'il faut garder ensuite
!           (-> on utilise IFINIA/IFINRA ensuite)

  call vorpre                                                     &
  !==========
 ( idbia1 , idbra1 , ifinia , ifinra ,                            &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nvar   , nscal  , nphas  , &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iirepv),                                                    &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,          &
   volume , ra(ipropc) , ra(ipropf) , ra(ipropb) ,                &
   rdevel , rtuser , ra     )

endif

! -- Fin de zone Methode des vortex pour la L.E.S.

! -- Structures mobiles en ALE

if (iale.eq.1) then

  idbia1 = ifinia
  idbra1 = ifinra

! Attention, strini reserve de la memoire qu'il faut garder ensuite
!           (-> on utilise IFINIA/IFINRA ensuite)
  call strini                                                     &
  !==========
 ( idbia1 , idbra1 , ifinia , ifinra ,                            &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume , ra(idt),                                     &
   rdevel , rtuser , ra     )

endif

! -- Fin de zone Structures mobiles en ALE

! -- Couplage Code_Saturne/Code_Saturne

call cscini                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,          &
   rdevel , rtuser , ra     )


!===============================================================================
! 10. DEBUT DE LA BOUCLE EN TEMPS
!===============================================================================

write(nfecra,2000)

ntcabs = ntpabs
ttcabs = ttpabs

iwarn0 = 1
do ivar = 1, nvar
  iwarn0 = max(iwarn0,iwarni(ivar))
enddo

if(iwarn0.gt.0) then
  write(nfecra,3000)
endif

if(inpdt0.eq.1) then
  ntmabs = ntcabs
endif

!     Nb d'iter ALE (nb relatif a l'execution en cours)
!     Si ITALIN=1, on fait une iteration d'initialisation
!     (si ITALIN=-999, c'est qu'on a fait une suite de calculs
!      sans relire lecamx -> comme si on faisait une suite
!      d'un calcul sans ALE)
if (italin.eq.-999) italin = 1
itrale = 1
if (italin.eq.1) then
  itrale = 0
  write(nfecra,3002) ttcabs
endif

!     En cas de couplage avec SYRTHES, on lit des maintenant l'entete
!     du premier message, au cas ou il s'agisse d'un message de
!     terminaison.

if (itrale.gt.0) then

  ntcam1 = ntcabs - 1

  call tstsyr (ntmabs, ntcam1)
  !==========

  if (ntmabs .eq. ntcam1) then
    call csexit (0)
    !==========
  endif

endif

 100  continue

if(inpdt0.eq.0 .and. itrale.gt.0) then
  ntcabs = ntcabs + 1
  if(idtvar.eq.0.or.idtvar.eq.1) then
    ttcabs = ttcabs + ra(idt)
  else
    ttcabs = ttcabs + dtref
  endif
  if(iwarn0.gt.0) then
    write(nfecra,3001) ttcabs,ntcabs
  endif
endif


!===============================================================================
! 11. AVANCEE EN TEMPS
!===============================================================================


!     On teste la presence de ficstp pour modifier NTMABS le cas echeant
call modpar(ntcabs,ntmabs)
!==========

call dmtmps(titer1)
!==========

!     Synchronisation Syrthes 3, si ITRALE>0
if (itrale.gt.0) then
  call itdsyr(ntcabs,ntmabs)
  !==========
endif

call tridim                                                       &
!==========
 ( ifinia , ifinra , itrale ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iisstd),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(itslag) , ra(icoefa) , ra(icoefb) ,                         &
   ra(ifrcx)  ,                                                   &
   rdevel , rtuser , ra     )

!===============================================================================
! 12. CALCUL DES MOYENNES TEMPORELLES (ACCUMULATION)
!===============================================================================


if(inpdt0.eq.0 .and. itrale.gt.0) then
call calmom                                                       &
!==========
     ( ifinia , ifinra , ncel   , ncelet ,                        &
       nideve , nrdeve , nituse , nrtuse ,                        &
       idevel , ituser , ia     ,                                 &
       ra(irtp  ) , ra(idt   ) , ra(ipropc) ,                     &
       rdevel , rtuser , ra     )
endif


!===============================================================================
! 13. APPEL DU MODULE LAGRANGIEN
!===============================================================================

if (iilagr.gt.0 .and. inpdt0.eq.0 .and. itrale.gt.0) then

  call memla2                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   nfabor , ncelet , nfac   ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   iindep , iibord , iettpa , iauxl  , iauxl2 ,                   &
   itaup  , iitlag , ipiil  ,                                     &
   ivagau , itsuf  , itsup  , ibx    ,                            &
   igradp , igradv , icroul ,                                     &
   itepct , itsfex , itsvar ,                                     &
   icpgd1 , icpgd2 , icpght ,                                     &
   ibrgau , itebru ,                                              &
   iw1    , iw2    , iw3    ,                                     &
   ilagia , ilagra  )

  call lagune                                                     &
  !==========
 ( ilagia , ilagra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicoce) , ia(iityce) , ia(iifrla) , ia(iiitep) ,            &
   ia(iindep) , ia(iibord) ,                                      &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iettpa) , ra(iitepa) , ra(istatc) , ra(istatv),&
   ra(itslag),                                                    &
   ra(istatf) , ra(itaup)  , ra(iitlag) , ra(ipiil)  , ra(ibx  ) ,&
   ra(ivagau) , ra(itsuf ) , ra(itsup ) , ra(itsvar) ,            &
   ra(itepct) , ra(itsfex) ,                                      &
   ra(icpgd1) , ra(icpgd2) , ra(icpght) ,                         &
   ra(igradp) , ra(igradv) , ra(icroul) ,                         &
   ra(ibrgau) , ra(itebru) ,                                      &
   ra(iw1   ) , ra(iw2   ) , ra(iw3   ) , ra(iauxl)  , ra(iauxl2),&
   rdevel , rtuser , ra     )

!--> Ici on libere la memoire reserve par MEMLA2
!      (i.e. on oublie ILAGIA et ILAGRA)

endif

!===============================================================================
! 14. BRANCHEMENT UTILISATEUR POUR MODIF DES VARIABLES EVENTUELLES
!===============================================================================

!     Appel pour Matisse d'une routine de bilans, d'impression, ...
if (imatis.eq.1 .and. itrale.gt.0) then

  call mtproj                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iiitep),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iettpa) , ra(iitepa) , ra(istatc) , ra(itslag),&
   ra(istatf) ,                                                   &
   rdevel , rtuser , ra     )

endif

if (itrale.gt.0) then

!       Sortie postprocessing de profils 1D

  if (iihmpr.eq.1) then
    call uiprof                                                   &
    !==========
  ( ncelet , ncel,                                                &
    ntmabs, ntcabs, ttcabs,                                       &
    xyzcen, ra(irtp), ra(ipropc) )
  endif

  ils    = ifinia
  ifnia2 = ils + maxelt
  call iasize('caltri',ifnia2)

  call usproj                                                     &
  !==========
 ( ifnia2 , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iiitep),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iettpa) , ra(iitepa) , ra(istatc) , ra(istatv),&
   ra(itslag) , ra(istatf) ,                                      &
   rdevel , rtuser , ra     )

endif

!===============================================================================
! 15. MISE A JOUR DU MAILLAGE (ALE)
!===============================================================================

if (iale.eq.1 .and. inpdt0.eq.0) then

  if (itrale.eq.0 .or. itrale.gt.nalinf) then

    call alemaj                                                   &
    !==========
 ( ifinia , ifinra , itrale ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iimpal)      ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(idepal) , ra(ixyzn0) ,                                      &
   rdevel , rtuser , ra     )

  endif

endif

!===============================================================================
! 16. TEST D'ARRET PAR MANQUE DE TEMPS
!===============================================================================

call armtps(ntcabs,ntmabs)
!==========

!===============================================================================
! 17. TEST D'ARRET ISSU DU MODULE DE RAYONNEMENT P-1
!===============================================================================
if (istpp1.eq.1) then
  ntmabs = ntcabs
endif

!===============================================================================
! 18. TEST D'ARRET PAR DEMANDE DE SYRTHES
!===============================================================================

!     En cas de couplage, on lit des maintenant l'entete du premier
!     message du pas de temps suivant, au cas ou il s'agisse d'un
!     message de terminaison (pas de test sur ITRALE ici, car
!     il serait sur ITRALE + 1, toujours > 0).

call tstsyr (ntmabs, ntcabs)
!==========

!===============================================================================
! 19. SORTIE EVENTUELLE DU FICHIER SUITE
!      (SAUF SI ON EST AU DERNIER PAS DE TEMPS :
!       ON LE FERA APRES)
!===============================================================================

iisuit = 0
if(ntcabs.lt.ntmabs) then
  if(ntsuit.eq.0) then
    ntsdef = max((ntmabs-ntpabs)/4,10)
    if(ntsdef.gt.0) then
      if(mod(ntcabs-ntpabs,ntsdef).eq.0) then
        iisuit = 1
      endif
    endif
  elseif(ntsuit.gt.0) then
    if(mod(ntcabs,ntsuit).eq.0) then
      iisuit = 1
    endif
  endif
endif
if (itrale.eq.0) iisuit = 0

if(iisuit.eq.1) then
  call ecrava                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nnod   ,          &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen     , surfac     , surfbo     , cdgfac     , cdgfbo    ,&
   ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb),&
   ra(icoefa) , ra(icoefb) , ra(ifrcx)  ,                         &
   rdevel , rtuser , ra     )

  if (nfpt1t.gt.0) then
    call ecrt1d                                                   &
    !==========
 ( ficvt1   , len(ficvt1), nfpt1d   ,  nmxt1d  ,                  &
   nfabor   , ra(itppt1) , ia(iifpt1))
  endif

  if (ippmod(iaeros).ge.0) then
     call ecrctw ( ficvct , len(ficvct) )
     !==========
  endif

  if (iilagr.gt.0) then

    call lagout                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicoce) , ia(iityce) , ia(iiitep) ,                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatf) ,                         &
   ra(istatc) , ra(istatv) , ra(itslag) ,                         &
   rdevel , rtuser , ra     )

  endif

  if (iirayo.gt.0) then
    call rayout                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml,  &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )
  endif

  if(iwarn0.gt.0) then
    write(nfecra,3020)ntcabs,ttcabs
  endif

endif

!===============================================================================
! 20. TEST POUR SAVOIR SI ON SORT UN FICHIER POST OU NON
!===============================================================================

call usnpst                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) , ra(istatc) ,                         &
   rdevel , rtuser , ra     )

!===============================================================================
! 21. SORTIE DES FICHIERS POST STANDARDS
!===============================================================================

!     Si ITRALE=0 on desactive tous les writers (car la geometrie n'a pas ete
!       ecrite)
if (itrale.eq.0) then
  indwri = 0
  indact = 0
  call pstact(indwri, indact)
  !==========
endif

call pstvar                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ntcabs , ncelet , ncel   , nfac   , nfabor ,          &
   nfml   , nprfml , nnod   , lndfac , lndfbr , ncelbr ,          &
   nvar   , nscal  , nphas  , nvlsta , nvisbr ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   ttcabs , xyzcen , surfac , surfbo , cdgfac , cdgfbo ,          &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(istatc) , ra(istatv) , ra(istatf) ,                         &
   rdevel , rtuser , ra     )

!===============================================================================
! 22. HISTORIQUES
!===============================================================================

! ON STOCKE SUR TMP SI ON A DECIDE DE LE FAIRE ET QUE C'EST LE MOMENT
! ON SAUVE SI ON A DECIDE DE LE FAIRE, QUE C'EST LE MOMENT ET
!   QUE CE N'EST PAS LE DERNIER PAS DE TEMPS

if(nthist.gt.0 .and. itrale.gt.0) then
  if(mod(ntcabs,nthist).eq.0) then
    modhis = 0
    if(nthsav.gt.0) then
      if(mod(ntcabs,nthsav).lt.nthist ) modhis = 1
    elseif(nthsav.eq.0) then
      ntcrel =  ntcabs-ntpabs
      nthdef = (ntmabs-ntpabs)/4
      if( (ntcrel-10).ge.0 .and. (ntcrel-10).lt.nthist) then
        modhis = 1
      elseif(nthdef.gt.0) then
        if(mod(ntcrel,nthdef).lt.nthist) then
          modhis = 1
        endif
      endif
    endif
    if (ntcabs.eq.ntmabs) modhis=0

    call ecrhis                                                   &
    !==========
   (ifinia , ifinra , ndim   , ncelet , ncel,                     &
    nideve , nrdeve , nituse , nrtuse ,                           &
    modhis ,                                                      &
    idevel , ituser , ia     ,                                    &
    xyzcen ,                                                      &
    rdevel , rtuser , ra     )

    if (iilagr.gt.0) then
      call laghis                                                 &
      !==========
     (ifinia , ifinra , ndim   , ncelet , ncel,                   &
      nideve , nrdeve , nituse , nrtuse ,                         &
      modhis , nvlsta ,                                           &
      idevel , ituser , ia     ,                                  &
      xyzcen , volume ,                                           &
      ra(istatc) , ra(istatv) ,                                   &
      rdevel , rtuser , ra     )
    endif

    if (ihistr.eq.1) then
      call strhis                                                 &
      !==========
     (ifinia , ifinra , ncelet , ncel,                            &
      nideve , nrdeve , nituse , nrtuse ,                         &
      modhis ,                                                    &
      idevel , ituser , ia     ,                                  &
      rdevel , rtuser , ra     )
    endif

  endif
endif

call ushist                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )

!===============================================================================
! 23. ECRITURE LISTING TOUTES LES NTLIST ITERATIONS
!===============================================================================


if(ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif
if(modntl.eq.0) then
  call ecrlis                                                     &
  !==========
  ( ifinia , ifinra ,                                             &
    nvar   , nphas  , ndim   , ncelet , ncel   ,                  &
    nideve , nrdeve , nituse , nrtuse , irtp   ,                  &
    idevel , ituser , ia     ,                                    &
    ra(irtp  ) , ra(irtpa ) , ra(idt ) , volume , xyzcen,         &
    rdevel , rtuser , ra     )

  if (iilagr.gt.0) then

    call laglis                                                   &
    !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iiitep),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatc) ,  ra(istatv) ,           &
   ra(itslag) , ra(istatf),                                       &
   rdevel , rtuser , ra     )

  endif

endif

call dmtmps(titer2)
!==========

if(iwarn0.gt.0) then
  if (itrale.gt.0) then
    write(nfecra,3010)ntcabs,titer2-titer1
  else
    write(nfecra,3012)titer2-titer1
  endif
endif


!===============================================================================
! 24. FIN DE LA BOUCLE EN TEMPS
!===============================================================================

itrale = itrale + 1

if(ntcabs.lt.ntmabs) goto 100


! LIBERATION DES TABLEAUX INTERMEDIAIRES (PDC+TSM)

ifinia = iditia
ifinra = iditra


!===============================================================================
! 25. ECRITURE DES SUITES + FINALISATION HISTORIQUES
!===============================================================================

call dmtmps(tecrf1)
!==========

if(iwarn0.gt.0) then
  write(nfecra,4000)
endif

call ecrava                                                       &
!==========
( ifinia , ifinra ,                                               &
  ndim   , ncelet , ncel   , nfac   , nfabor , nnod   ,           &
  nvar   , nscal  , nphas  ,                                      &
  nideve , nrdeve , nituse , nrtuse ,                             &
  idevel , ituser , ia     ,                                      &
  xyzcen     , surfac     , surfbo     , cdgfac     , cdgfbo    , &
  ra(idt)    , ra(irtp)   , ra(ipropc) , ra(ipropf) , ra(ipropb), &
  ra(icoefa) , ra(icoefb) , ra(ifrcx)  ,                          &
  rdevel , rtuser , ra     )

if (nfpt1t.gt.0) then

  call ecrt1d                                                     &
  !==========
 ( ficvt1   , len(ficvt1), nfpt1d   ,  nmxt1d  ,                  &
   nfabor   , ra(itppt1) , ia(iifpt1))
endif

if (ippmod(iaeros).ge.0) then
  call ecrctw ( ficvct , len(ficvct) )
  !==========
endif

if (iilagr.gt.0) then

  call lagout                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ia(iicoce) , ia(iityce) , ia(iiitep) ,                         &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   ra(iettp)  , ra(iitepa) , ra(istatf) ,                         &
   ra(istatc) , ra(istatv) , ra(itslag) ,                         &
   rdevel , rtuser , ra     )

endif

! Ecriture du fichier suite
if (iirayo.gt.0) then

  call rayout                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml,                    &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   ra(idt)    , ra(irtpa)  , ra(irtp)   ,                         &
   ra(ipropc) , ra(ipropf) , ra(ipropb) ,                         &
   ra(icoefa) , ra(icoefb) ,                                      &
   rdevel , rtuser , ra     )

endif

! ICI ON SAUVE LES HISTORIQUES (SI ON EN A STOCKE)

modhis = 2
call ecrhis                                                       &
!==========
(  ifinia , ifinra , ndim   , ncelet , ncel,                      &
   nideve , nrdeve , nituse , nrtuse ,                            &
   modhis ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen ,                                                       &
   rdevel , rtuser , ra     )

if (iilagr.gt.0) then
    call laghis                                                   &
    !==========
   (ifinia , ifinra , ndim   , ncelet , ncel,                     &
    nideve , nrdeve , nituse , nrtuse ,                           &
    modhis , nvlsta ,                                             &
   idevel , ituser , ia     ,                                     &
    xyzcen , volume ,                                             &
    ra(istatc) , ra(istatv) ,                                     &
   rdevel , rtuser , ra     )
endif
if (ihistr.eq.1) then
  call strhis                                                     &
  !==========
 (ifinia , ifinra , ncelet , ncel,                                &
  nideve , nrdeve , nituse , nrtuse ,                             &
  modhis ,                                                        &
  idevel , ituser , ia     ,                                      &
  rdevel , rtuser , ra     )
endif

!     LE CAS ECHEANT, ON LIBERE LES STRUCTURES C DU MODULE THERMIQUE 1D
!     ET/OU ON FERME LE LISTING LAGRANGIEN

if (nfpt1d.gt.0) then
  call lbrt1d
  !==========
endif

if (iale.gt.0) then
  call lbrale
  !==========
endif

if (iilagr.gt.0) close(implal)

call dmtmps(tecrf2)
!==========

if(iwarn0.gt.0) then
  write(nfecra,4010)tecrf2-tecrf1
endif

 200  continue

!===============================================================================
! 26. MEMOIRE UTILISEE
!===============================================================================

!     Libération des structures liées à la lecture du fichier xml

if (iihmpr.eq.1) then

  call memui2
  !==========

  call memui1(ncharb)
  !==========
endif


write(nfecra,7000)

!----
! FORMATS
!----

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'===============================================================',&
                                                              /,/,&
'              FONCTIONNEMENT EN MODE VERIFICATION            ',/,&
'              ===================================            ',/,&
'                                                             ',/,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)

 2000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                       CORPS DU CALCUL                       ',/,&
'                       ===============                       ',/,&
                                                                /,&
                                                                /,&
'===============================================================',&
                                                              /,/,&
                                                                /)
 3000 format(/,                                                   &
'===============================================================',&
 /)
 3001 format(/,' INSTANT ',E18.9,        '   ITERATION NUMERO ',I15,/,  &
' ============================================================= ',&
 /,/)
 3002 format(/,' INSTANT ',E18.9,        '   INITIALISATION ALE ',/,    &
' ============================================================= ',&
 /,/)
 3010 format(/,' TEMPS CPU POUR L''ITERATION ',I15,' :    ',E14.5,/,/,  &
'===============================================================',&
 /)
 3012 format(/,' TEMPS CPU POUR L''INITIALISATION ALE :    ',E14.5,/,/, &
'===============================================================',&
 /)
 3020 format(/,/,                                                 &
 ' Sortie intermediaire d''un fichier suite ',/,            &
 '   Sauvegarde a l''iteration ',    I10,                         &
                                    ', Temps physique ',E14.5,/,/)

 4000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                   ETAPES FINALES DU CALCUL                  ',/,&
'                   ========================                  ',/,&
                                                                /,&
                                                                /,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)
 4010 format(                                                         /,&
 3X,'** TEMPS CPU POUR LES SORTIES FINALES : ',E14.5           ,/,&
 3X,'   ----------------------------------                    ',/)
 7000 format(/,/,                                                 &
' =========================================================== ',/,&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                 FIN DE L''EXECUTION DU CALCUL               ',/,&
'                 ============================                ',/,&
                                                                /,&
                                                                /,&
'===============================================================')

#else

 1000 format(/,                                                   &
'===============================================================',&
                                                              /,/,&
'              RUNNING IN VERIFICATION MODE                   ',/,&
'              ============================                   ',/,&
'                                                             ',/,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)

 2000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                       MAIN CALCULATION                      ',/,&
'                       ================                      ',/,&
                                                                /,&
                                                                /,&
'===============================================================',&
                                                              /,/,&
                                                                /)
 3000 format(/,                                                   &
'===============================================================',&
 /)
 3001 format(/,' INSTANT ',E18.9,        '   TIME STEP NUMBER ',I15,/,  &
' ============================================================= ',&
 /,/)
 3002 format(/,' INSTANT ',E18.9,        '   ALE INITIALIZATION ',/,    &
' ============================================================= ',&
 /,/)
 3010 format(/,' CPU TIME FOR THE TIME STEP  ',I15,':     ',E14.5,/,/,  &
'===============================================================',&
 /)
 3012 format(/,' CPU TIME FOR ALE INITIALIZATION:          ',E14.5,/,/, &
'===============================================================',&
 /)
 3020 format(/,/,                                                 &
 ' Write an intermediate restart file',/,                   &
 '   Backup at iteration ',    I10,  ', Physical time ',E14.5,/,/)

 4000 format(/,/,                                                 &
'===============================================================',&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                 FINAL STAGE OF THE CALCULATION              ',/,&
'                 ==============================              ',/,&
                                                                /,&
                                                                /,&
' =========================================================== ',/,&
                                                                /,&
                                                                /)
 4010 format(                                                         /,&
 3X,'** CPU TIME FOR FINAL WRITING: ',E14.5                    ,/,&
 3X,'   ---------------------------                           ',/)
 7000 format(/,/,                                                 &
' =========================================================== ',/,&
                                                              /,/,&
                                                                /,&
                                                                /,&
'                      END OF CALCULATION                     ',/,&
'                      ==================                     ',/,&
                                                                /,&
                                                                /,&
'===============================================================')

#endif

!===============================================================================
! 26. FIN
!===============================================================================


return
end

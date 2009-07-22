!-------------------------------------------------------------------------------

!VERS


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

subroutine ushist &
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
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! SORTIE D'HISTORIQUES NON STD LIVREE A L'UTILISATEUR

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
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "parall.h"
include "period.h"

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
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ii, kk, node, ndrang, nvarpp, numcel, lng
double precision xx, yy, zz, xyztmp(3)

!   Numero des noeuds ou on sort des historiques
integer          ncapmx
parameter       (ncapmx=100)
integer          icapt(ncapmx)
save             icapt
integer          ircapt(ncapmx)
save             ircapt

!   Nombre de noeuds ou on sort des historiques
integer          ncapts
save             ncapts

!   Numero du passage actuel dans ce ss pgm
integer          ipass
data             ipass /0/
save             ipass

!   Tableau de valeurs temporaires
double precision vacapt(ncapmx)

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. INITIALISATION
!===============================================================================

! ---> Gestion memoire

idebia = idbia0
idebra = idbra0

! ---> Numero du passage actuel dans ce ss pgm

ipass = ipass + 1

!===============================================================================
! 2. RECHERCHE DES CAPTEURS
!===============================================================================
! Les numeros stockes dans IRCAPT donnent le rang du processeur sur
!   lequel se trouve la sonde. L'utilisateur n'a pas a s'en preoccuper
!   specialement tant qu'il utilise bien la fonction FINDPT pour reperer
!   les sondes.


!  Au premier passage : reperage des numeros de cellule dont le centre est
!    le plus proche des coordonnees XX YY ZZ.
!    En parallelisme, le numero de cellule ICAPT(II) est local au processeur
!    dont le rang est donne par IRCAPT(II) (de 0 a nombre de processeurs-1).
!    NCAPTS donne le nombre de sondes total.

if (ipass.eq.1) then

  ii = 0

  xx = 0.20d0
  yy = 0.15d0
  zz = 0.01d0
  call findpt                                                     &
  !==========
 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  xx = 0.70d0
  yy = 0.15d0
  zz = 0.01d0
  call findpt                                                     &
  !==========
 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  xx = 0.20d0
  yy = 0.75d0
  zz = 0.01d0
  call findpt                                                     &
  !==========
 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  xx = 0.70d0
  yy = 0.75d0
  zz = 0.01d0
  call findpt                                                     &
  !==========
 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  ncapts = ii

  if(ii.gt.ncapmx) then
    WRITE(NFECRA,*) ' USHIST : NCAPMX = ',II,' AU MINIMUM '
    call csexit (1)
  endif

endif


!===============================================================================
! 3. OUVERTURE DES FICHIERS
!     EXEMPLE D'UNE VARIABLE PAR FICHIER
!===============================================================================

! ---> Nombre de variables = nombre de fichiers

nvarpp = nvar


! ---> Au premier passage, on ouvre les fichiers et on ecrit une entete

if(ipass.eq.1) then

!   --> Test du nombre max de fichiers

  if(nvarpp.gt.nushmx) then
    write(nfecra,*)                                               &
  ' USHIST : PAS DROIT A PLUS DE ',NUSHMX,' FICHIERS HISTORIQUES'
    call csexit (1)
  endif

  do ii = 1, nvarpp

!   --> Ouverture des fichiers avec les unites disponibles

    if (irangp.le.0) then
      open(file=ficush(ii),unit=impush(ii))
    endif

!   --> On imprime le numero (global) de la cellule et les coordonnees
!        du centre

    do kk = 1, ncapts
!           Numero de cellule (en parallele : local au processeur courant)
      numcel    = icapt(kk)
      if (irangp.lt.0 .or. irangp.eq.ircapt(kk)) then
!           Coordonnees de la cellule (en parallele, c'est le processeur
!             qui la contient qui travaille)
        xyztmp(1) = xyzcen(1,numcel)
        xyztmp(2) = xyzcen(2,numcel)
        xyztmp(3) = xyzcen(3,numcel)
      else
!           Valeurs bidons sur les autres processeurs
        xyztmp(1) = 0.d0
        xyztmp(2) = 0.d0
        xyztmp(3) = 0.d0
      endif
!           En parallele, le processeur qui a trouve la cellule
!           envoie son numero global et ses coordonnees aux autres.
      if (irangp.ge.0) then
        call parcel(icapt(kk), ircapt(kk), numcel)
        !==========
        lng = 3
        call parbcr(ircapt(kk), lng , xyztmp)
        !==========
      endif
!           On ecrit les informations (seul le processeur 0
!           travaille en parallele : on n'a pasa besoin de
!           plusieurs exemplaires du fichier)
      if (irangp.le.0) then
        WRITE(IMPUSH(II),1000) '#',' Cellule ',NUMCEL,            &
            ' Coord ',XYZTMP(1),XYZTMP(2),XYZTMP(3)
      endif

    enddo

  enddo

endif

 1000 format(a,a9,i10,a7,3e14.5)

!===============================================================================
! 4. ECRITURE
!     EXEMPLE D'UNE VARIABLE PAR FICHIER
!===============================================================================

! Ecriture du numero du pas de temps,
!          de la valeur du temps physique
!          de la variable en tous les points d'historiques
! En sequentiel, la valeur a ecrire est simplement RTP(ICAPT(KK),II)
! en parallele, la valeur a ecrire peut etre sur un autre processeur
!  et il faut la determiner dans  VACAPT(KK) avec PARHIS.

do ii = 1 , nvarpp
  do kk = 1, ncapts
    if (irangp.lt.0) then
      vacapt(kk) = rtp(icapt(kk),ii)
    else
      call parhis(icapt(kk), ircapt(kk), rtp(1,ii), vacapt(kk))
      !==========
    endif
  enddo
  if (irangp.le.0) then
    write (impush(ii),1010) ntcabs,ttcabs,                        &
                            (vacapt(kk),kk=1,ncapts)
  endif
enddo



!   ATTENTION : IL FAUT ADAPTER LE FORMAT POUR PLUS DE  9 CAPTEURS


 1010 format(i10,10e17.9)

!===============================================================================
! 4. FERMETURE
!===============================================================================

if(ntcabs.eq.ntmabs .and. irangp.le.0) then
  do ii = 1, nvarpp
    close(impush(ii))
  enddo
endif


!===============================================================================
! 5. SORTIE
!===============================================================================

return
end

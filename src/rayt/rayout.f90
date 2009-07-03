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

subroutine rayout &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
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
! ----------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  1) ECRITURE FICHIER SUITE,
!  2) Ecriture des fichiers Ensight pour les sorties sur les
!     frontieres du domaine de calcul

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
! lndnod           ! e  ! <-- ! longueur du tableau icocel (optionnel          !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
!                  !    !     ! le module lagrangien                           !
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
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"
include "ppincl.h"
include "radiat.h"
include "parall.h"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndnod , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse
integer          nvar   , nscal  , nphas  , iph

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
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

character        rubriq*64
character        cphase(nphsmx)*2
integer          idebia , idebra
integer          ifinia , ifinra
integer          itrav1 , ip
integer          ierror , nberro , irtyp , itysup , nbval
integer          ivers  , ilecec
integer          impavr

!===============================================================================
!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. ECRITURE DU FICHIER SUITE DU MODULE DE RAYONNEMENT
!===============================================================================


! ---> Ouverture (et on saute si erreur)
!     ILECEC = 2 : ecriture

write(nfecra,6010)

ilecec = 2
call opnsui(ficavr, len(ficavr), ilecec, impavr, ierror)
!==========
if (ierror.ne.0) then
  write(nfecra,9020)
  goto 9998
endif

write(nfecra,6011)

! Entete et Dimensions ou on saute si erreur
!     On inclut une rubrique destinee a distinguer ce fichier
!       d'un autre fichier suite
!     Pour le moment, IVERS n'est pas utilise

nberro = 0

ivers  = 111
itysup = 0
nbval  = 1
irtyp  = 1
RUBRIQ = 'version_fichier_suite_rayonnement'
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,ivers,   &
            ierror)
nberro=nberro+ierror

itysup = 0
nbval  = 1
irtyp  = 1

if(nberro.ne.0) then
  write(nfecra,9120)
  goto 9998
endif

write(nfecra,6012)

! Temps (par securite)

nberro = 0

RUBRIQ = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
irtyp  = 1
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,ntcabs,  &
            ierror)
nberro=nberro+ierror

RUBRIQ = 'instant_precedent'
itysup = 0
nbval  = 1
irtyp  = 2
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,ttcabs,  &
            ierror)
nberro=nberro+ierror

if(nberro.ne.0) then
  write(nfecra,8121)
endif

! Donnees

nberro = 0

!     Aux faces de bord

  itysup = 3
  nbval  = 1
  irtyp  = 2

  RUBRIQ = 'tparoi_fb'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propfb(1,ipprob(itparo)),ierror)
  nberro=nberro+ierror

  RUBRIQ = 'qincid_fb'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propfb(1,ipprob(iqinci)),ierror)
  nberro=nberro+ierror

  RUBRIQ = 'hfconv_fb'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propfb(1,ipprob(ihconv)),ierror)
  nberro=nberro+ierror

  RUBRIQ = 'flconv_fb'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propfb(1,ipprob(ifconv)),ierror)
  nberro=nberro+ierror


!     Aux cellules

  itysup = 1
  nbval  = 1
  irtyp  = 2

  RUBRIQ = 'rayimp_ce'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipproc(itsri(1))),ierror)
  nberro=nberro+ierror

  RUBRIQ = 'rayexp_ce'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipproc(itsre(1))),ierror)
  nberro=nberro+ierror

!  ---> Si pb : on saute

if(nberro.ne.0) then
  write(nfecra,9100)
  goto 9998
endif

write(nfecra,6013)

! ---> Fermeture du fichier suite
call clssui(impavr,ierror)

if (ierror.ne.0) then
  write(nfecra,8011) ficavr
endif

write(nfecra,6014)

! ---> En cas d'erreur, on continue quand meme
 9998 continue


return


!--------
! FORMATS
!--------

 6010 FORMAT (/, 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT ',/,  &
           3X,'   ------------------------------------------',/,  &
           3X,' Ecriture d''un fichier suite                ',/)

 6011 FORMAT (   3X,'   Debut de l''ecriture                      ',/)
 6012 FORMAT (   3X,'   Fin de l''ecriture des dimensions         ',/)
 6013 FORMAT (   3X,'   Fin de l''ecriture des donnees            ',/)
 6014 FORMAT (   3X,' Fin de l''ecriture du fichier suite         ',/)

 9020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@    ERREUR A L''OUVERTURE DU FICHIER SUITE RAYONNEMENT      ',/,&
'@                                                            ',/,&
'@  Le calcul continue mais                                   ',/,&
'@            ne fournira pas de fichier suite rayonnement.   ',/,&
'@                                                            ',/,&
'@  Verifier que le repertoire de travail est accessible en   ',/,&
'@    ecriture et que le fichier suite peut y etre cree.      ',/,&
'@  Voir le sous-programme rayout.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9120 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE L''ECRITURE DES DIMENSIONS             ',/,&
'@                                                            ',/,&
'@  Le calcul continue mais                                   ',/,&
'@            ne fournira pas de fichier suite rayonnement.   ',/,&
'@                                                            ',/,&
'@  Voir le sous-programme rayout.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE L''ECRITURE DU PAS DE TEMPS ET DU TEMPS',/,&
'@                                                            ',/,&
'@    Le calcul continue...                                   ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme rayout.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: A L''ECRITURE DU FICHIER SUITE RAYONNEMENT   ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERREUR LORS DE L''ECRITURE DES DONNEES                ',/,&
'@                                                            ',/,&
'@  Le calcul continue mais                                   ',/,&
'@            ne fournira pas de fichier suite rayonnement.   ',/,&
'@                                                            ',/,&
'@  Voir le sous-programme rayout.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========                              AVAL RAYONNMEMENT',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom (',A13,')                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7000 format(                                                           &
'                                                             ',/,&
' Module de rayonnement :                                     ',/,&
'     Ecriture du fichier Ensight de bord non disponible      ',/,&
'     en parallele.                                           ',/,&
'                                                             ',/)

!----
! FIN
!----

end

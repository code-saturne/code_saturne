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

subroutine lagaff &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , statis , &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!       SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!       -----------------------------------

!  ECRITURE SUR FICHIERS DES INFORMATIONS SUR LE NOMBRE DE PARTICULES
!       - nombre de particules dans le domaine
!       - nombre de particules entrantes
!       - nombre de particules sorties



!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa(nbpmax,     ! tr ! <-- ! caracteristiques des particules                !
!       nvep)      !    !     !  aux particules (poids, ...)                   !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! statis(ncelet    ! tr ! <-- ! cumul des statistiques volumiques              !
!    nvlsta)       !    !     !                                                !
! w1..w3(ncelet    ! tr ! --- ! tableaux de travail                            !
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
use numvar
use cstnum
use optcal
use pointe
use entsor
use parall
use lagpar
use lagran
use cstphy
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision statis(ncelet,nvlsta)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision ra(*)

! Local variables

integer          idebia , idebra

integer          iphas
double precision dnbpr

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================

iphas = 1
ipass = ipass + 1

!===============================================================================
! 2. OUVERTURE DU FICHIER DE STOCKAGE
!===============================================================================

!     Seul le premier processeur ecrit les informations
if (irangp.le.0) then

  if (ipass.eq.1 ) then

    if ( iroule .ge. 1 .and.                                      &
         (iphyla.eq.2 .and. iencra.eq.1) ) then
      write(implal,1000)
    elseif ( iroule .ge. 1 .and.                                  &
           (iphyla.ne.2 .or. iencra.ne.1) ) then
      write(implal,1001)
    elseif ( iroule .ne. 1 .and.                                  &
           (iphyla.eq.2 .and. iencra.eq.1) ) then
      write(implal,1002)
    else
      write(implal,1003)
    endif

  endif

!===============================================================================
! 2 - Ecriture des INFORMATIONS
!===============================================================================

  if (nbptot.gt.0) then
    dnbpr = (nbpert*100.d0)/dble(nbptot)
  else
    dnbpr = 0
  endif

  if ( iroule.ge.1 .and.                                          &
       (iphyla.eq.2 .and. iencra.eq.1) ) then

    write(implal,2000) iplas,(dtp*iplas),                         &
         nbpart        , dnbpar        ,                          &
         nbpnew        , dnbpnw        ,                          &
         nbpout-nbperr , dnbpou-dnbper ,                          &
         nbperr        , dnbper        ,                          &
         dnbpr         ,                                          &
         npcsup        , dnpcsu        ,                          &
         npclon        , dnpclo        ,                          &
         npkill        , dnpkil        ,                          &
         npencr        , dnpenc

  elseif ( iroule.ge.1 .and.                                      &
         (iphyla.ne.2 .or. iencra.ne.1) ) then

    write(implal,2001) iplas,(dtp*iplas),                         &
         nbpart        , dnbpar        ,                          &
         nbpnew        , dnbpnw        ,                          &
         nbpout-nbperr , dnbpou-dnbper ,                          &
         nbperr        , dnbper        ,                          &
         dnbpr         ,                                          &
         npcsup        , dnpcsu        ,                          &
         npclon        , dnpclo        ,                          &
         npkill        , dnpkil

  elseif ( iroule.lt.1 .and.                                      &
         (iphyla.eq.2 .and. iencra.eq.1) ) then

    write(implal,2002) iplas,(dtp*iplas),                         &
         nbpart        , dnbpar        ,                          &
         nbpnew        , dnbpnw        ,                          &
         nbpout-nbperr , dnbpou-dnbper ,                          &
         nbperr        , dnbper        ,                          &
         dnbpr         ,                                          &
         npencr        , dnpenc

  else

    write(implal,2003) iplas,(dtp*iplas),                         &
         nbpart        , dnbpar        ,                          &
         nbpnew        , dnbpnw        ,                          &
         nbpout-nbperr , dnbpou-dnbper ,                          &
         nbperr        , dnbper        ,                          &
         dnbpr

  endif

endif

!===============================================================================


!--------
! FORMATS
!--------

 1000 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     ',/, &
       '#                                         ',/,      &
       '# colonne  1 : numero de pas de temps      ',/,     &
       '# colonne  2 : temps physique              ',/,     &
       '# colonne  3 : nbre inst. de part.    ',/,          &
       '# colonne  4 : nbre inst. de part. (avec poids)',/, &
       '# colonne  5 : nbre inst. de part. injectees',/,    &
       '# colonne  6 : nbre inst. de part. injectees',            &
       ' (avec poids)',/,                                   &
       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,&
       '# colonne  8 : nbre inst. de part. sorties ou deposees',  &
       ' (avec poids)',/,                                   &
       '# colonne  9 : nbre inst. de part. perdues (reperage) ',/,&
       '# colonne 10 : nbre inst. de part. perdues',              &
        ' (reperage, avec poids)',/,                        &
       '# colonne 11 : % de part. perdues',/,               &
       '# colonne 12 : nbre inst. de part. qui ont subi le',      &
       ' clonage',/,                                        &
       '# colonne 13 : nbre inst. de part. qui ont subi le',      &
       ' clonage (avec poids)',/,                           &
       '# colonne 14 : nbre inst. de nouvel. part. par clonage',/,&
       '# colonne 15 : nbre inst. de nouvel. part. par clonage',  &
       ' (avec poids)',/,                                   &
       '# colonne 16 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe ',/,                                &
       '# colonne 17 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe (avec poids)',/,                    &
       '# colonne 18 : nbre inst. de part encrassees',            &
       ' (Charbon)) '/,                                     &
       '# colonne 19 : nbre inst. de part encrassees',            &
       ' (Charbon, avec poids))',/,                         &
       '# ')

 1001 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     ',/, &
       '#                                         ',/,      &
       '# colonne  1 : numero de pas de temps      ',/,     &
       '# colonne  2 : temps physique              ',/,     &
       '# colonne  3 : nbre inst. de part.    ',/,          &
       '# colonne  4 : nbre inst. de part. (avec poids)',/, &
       '# colonne  5 : nbre inst. de part. injectees',/,    &
       '# colonne  6 : nbre inst. de part. injectees',            &
       ' (avec poids)',/,                                   &
       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,&
       '# colonne  8 : nbre inst. de part. sorties ou deposees',  &
       ' (avec poids)',/,                                   &
       '# colonne  9 : nbre inst. de part. perdues (reperage) ',/,&
       '# colonne 10 : nbre inst. de part. perdues',              &
        ' (reperage, avec poids)',/,                        &
       '# colonne 11 : % de part. perdues',/,               &
       '# colonne 12 : nbre inst. de part. qui ont subi le',      &
       ' clonage',/,                                        &
       '# colonne 13 : nbre inst. de part. qui ont subi le',      &
       ' clonage (avec poids)',/,                           &
       '# colonne 14 : nbre inst. de nouvel. part. par clonage',/,&
       '# colonne 15 : nbre inst. de nouvel. part. par clonage',  &
       ' (avec poids)',/,                                   &
       '# colonne 16 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe ',/,                                &
       '# colonne 17 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe (avec poids)',/,                    &
       '# ')

 1002 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     ',/, &
       '#                                         ',/,      &
       '# colonne  1 : numero de pas de temps      ',/,     &
       '# colonne  2 : temps physique              ',/,     &
       '# colonne  3 : nbre inst. de part.    ',/,          &
       '# colonne  4 : nbre inst. de part. (avec poids)',/, &
       '# colonne  5 : nbre inst. de part. injectees',/,    &
       '# colonne  6 : nbre inst. de part. injectees',            &
       ' (avec poids)',/,                                   &
       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,&
       '# colonne  8 : nbre inst. de part. sorties ou deposees',  &
       ' (avec poids)',/,                                   &
       '# colonne  9 : nbre inst. de part. perdues (reperage)',/, &
       '# colonne 10 : nbre inst. de part. perdues',              &
        ' (reperage, avec poids)',/,                        &
       '# colonne 11 : % de part. perdues ',/,              &
       '# colonne 12 : nbre inst. de part. encrassees',           &
       ' (Charbon)) '/,                                     &
       '# colonne 13 : nbre inst. de part. encrassees',           &
       ' (Charbon, avec poids))',/,                         &
       '# ')

 1003 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     ',/, &
       '#                                         ',/,      &
       '# colonne  1 : numero de pas de temps      ',/,     &
       '# colonne  2 : temps physique              ',/,     &
       '# colonne  3 : nbre inst. de part.    ',/,          &
       '# colonne  4 : nbre inst. de part. (avec poids)',/, &
       '# colonne  5 : nbre inst. de part. injectees',/,    &
       '# colonne  6 : nbre inst. de part. injectees',            &
       ' (avec poids)',/,                                   &
       '# colonne  7 : nbre inst. de part. sorties ou deposees',/,&
       '# colonne  8 : nbre inst. de part. sorties ou deposees',  &
       ' (avec poids)',/,                                   &
       '# colonne  9 : nbre inst. de part. perdues (reperage)',/, &
       '# colonne 10 : nbre inst. de part. perdues',              &
        ' (reperage, avec poids)',/,                        &
       '# colonne 11 : % de part. perdues ',/,              &
       '# ')

 2000 format(1x,i8,2x,e10.4,2x,4(i8,2x,e10.4),2x,e10.4,4(i8,2x,e10.4))
 2001 format(1x,i8,2x,e10.4,2x,4(i8,2x,e10.4),2x,e10.4,3(i8,2x,e10.4))
 2002 format(1x,i8,2x,e10.4,2x,4(i8,2x,e10.4),2x,e10.4,1(i8,2x,e10.4))
 2003 format(1x,i8,2x,e10.4,2x,4(i8,2x,e10.4),2x,e10.4)

!====
! FIN
!====

return

end subroutine

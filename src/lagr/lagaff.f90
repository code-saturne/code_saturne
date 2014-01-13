!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine lagaff
!================

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

! Local variables

double precision dnbpr

integer nbpartall, nbpoutall, nbperrall, nbpdepall, npencrall, nbpresall

double precision dnbparall, dnbperall, dnbpouall
double precision dnbdepall, dnpencall, dnbpnwall, dnbresall

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass



!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

ipass = ipass + 1

! Parallelism management

nbpartall = nbpart
nbpoutall = nbpout
nbperrall = nbperr
nbpdepall = nbpdep
npencrall = npencr
nbpresall = nbpres

dnbparall = dnbpar
dnbpouall = dnbpou
dnbperall = dnbper
dnbdepall = dnbdep
dnpencall = dnpenc
dnbpnwall = dnbpnw
dnbresall = dnbres


if (irangp.ge.0) then

   call parcpt(nbpartall)
   call parcpt(nbpoutall)
   call parcpt(nbperrall)
   call parcpt(nbpdepall)
   call parcpt(npencrall)
   call parcpt(nbpresall)

   call parsom(dnbparall)
   call parsom(dnbpouall)
   call parsom(dnbperall)
   call parsom(dnbdepall)
   call parsom(dnpencall)
   call parsom(dnbpnwall)
   call parsom(dnbresall)

endif


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
    elseif (ireent.gt.0) then
      write(implal,1004)
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

    write(implal,2000) iplas,ttcabs,                              &
         nbpartall        , dnbparall        ,                    &
         nbpnew        ,dnbpnwall        ,                        &
         nbpoutall-nbperrall-npencrall                            &
                       , dnbpouall-dnbperall-dnpencall ,          &
         nbpdepall        , dnbdepall        ,                    &
         nbperrall        , dnbperall        ,                    &
         dnbpr         ,                                          &
         npcsup        , dnpcsu        ,                          &
         npclon        , dnpclo        ,                          &
         npkill        , dnpkil        ,                          &
         npencrall        , dnpencall

  elseif ( iroule.ge.1 .and.                                      &
         (iphyla.ne.2 .or. iencra.ne.1) ) then

    write(implal,2001) iplas,ttcabs,                              &
         nbpartall     , dnbparall        ,                       &
         nbpnew        ,dnbpnwall        ,                        &
         nbpoutall-nbperrall , dnbpouall-dnbperall ,              &
         nbpdepall        , dnbdepall        ,                    &
         nbperrall        , dnbperall        ,                    &
         dnbpr         ,                                          &
         npcsup        , dnpcsu        ,                          &
         npclon        , dnpclo        ,                          &
         npkill        , dnpkil

  elseif ( iroule.lt.1 .and.                                      &
         (iphyla.eq.2 .and. iencra.eq.1) ) then

    write(implal,2002) iplas,ttcabs,                              &
         nbpartall     , dnbparall        ,                       &
         nbpnew        ,dnbpnwall        ,                        &
         nbpoutall-nbperrall-npencrall                            &
                       , dnbpouall-dnbperall-dnpencall ,          &
         nbpdepall        , dnbdepall        ,                    &
         nbperrall        , dnbperall        ,                    &
         dnbpr         ,                                          &
         npencrall     , dnpencall

  elseif (ireent.gt.0) then

    write(implal,2004) iplas,ttcabs,                              &
         nbpartall     , dnbparall        ,                       &
         nbpnew        ,dnbpnwall        ,                        &
         nbpoutall-nbperrall , dnbpouall-dnbperall ,              &
         nbpdepall        , dnbdepall        ,                    &
         nbpresall        , dnbresall        ,                    &
         nbperrall        , dnbperall        ,                    &
         dnbpr

  else

    write(implal,2003) iplas  ,ttcabs,                            &
         nbpartall     , dnbparall        ,                       &
         nbpnew        ,dnbpnwall        ,                        &
         nbpoutall-nbperrall , dnbpouall-dnbperall ,              &
         nbpdepall        , dnbdepall        ,                    &
         nbperrall        , dnbperall        ,                    &
         dnbpr

  endif

endif

!===============================================================================


!--------
! FORMATS
!--------

 1000 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     '      ,/, &
       '#                                              '      ,/, &
       '# colonne  1 : numero de pas de temps          '      ,/, &
       '# colonne  2 : temps physique                  '      ,/, &
       '# colonne  3 : nbre inst. de part.             '      ,/, &
       '# colonne  4 : nbre inst. de part. (avec poids)'      ,/, &
       '# colonne  5 : nbre inst. de part. injectees   '      ,/, &
       '# colonne  6 : nbre inst. de part. injectees   '      ,   &
       ' (avec poids)'                                        ,/, &
       '# colonne  7 : nbre inst. de part. sorties, ou deposees et supprimees',/,&
       '# colonne  8 : nbre inst. de part. sorties, ou deposees et supprimees',  &
       ' (avec poids)'                                        ,/, &
       '# colonne  9 : nbre inst. de part. deposees '         ,/, &
       '# colonne 10 : nbre inst. de part. deposees '         ,   &
       ' (avec poids)'                                        ,/, &
       '# colonne 11 : nbre inst. de part. perdues (reperage)',/, &
       '# colonne 12 : nbre inst. de part. perdues'           ,   &
       ' (reperage, avec poids)'                              ,/, &
       '# colonne 13 : % de part. perdues'                    ,/, &
       '# colonne 14 : nbre inst. de part. qui ont subi le'   ,   &
       ' clonage'                                             ,/, &
       '# colonne 15 : nbre inst. de part. qui ont subi le'   ,   &
       ' clonage (avec poids)'                                ,/, &
       '# colonne 16 : nbre inst. de nouvel. part. par clonage',/,&
       '# colonne 17 : nbre inst. de nouvel. part. par clonage',  &
       ' (avec poids)'                                        ,/, &
       '# colonne 18 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe '                                     ,/, &
       '# colonne 19 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe (avec poids)'                         ,/, &
       '# colonne 20 : nbre inst. de part encrassees'         ,   &
       ' (Charbon) '                                          ,/, &
       '# colonne 21 : nbre inst. de part encrassees'         ,   &
       ' (Charbon, avec poids)'                               ,/, &
       '# ')

 1001 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     '      ,/, &
       '#                                              '      ,/, &
       '# colonne  1 : numero de pas de temps          '      ,/, &
       '# colonne  2 : temps physique                  '      ,/, &
       '# colonne  3 : nbre inst. de part.             '      ,/, &
       '# colonne  4 : nbre inst. de part. (avec poids)'      ,/, &
       '# colonne  5 : nbre inst. de part. injectees   '      ,/, &
       '# colonne  6 : nbre inst. de part. injectees   '        , &
       ' (avec poids)'                                        ,/, &
       '# colonne  7 : nbre inst. de part. sorties, ou deposees et supprimees',/,&
       '# colonne  8 : nbre inst. de part. sorties, ou deposees et supprimees',  &
       ' (avec poids)'                                        ,/, &
       '# colonne  9 : nbre inst. de part. deposees    '      ,/, &
       '# colonne 10 : nbre inst. de part. deposees    '        , &
       ' (avec poids)'                                        ,/, &
       '# colonne 11 : nbre inst. de part. perdues (reperage) ',/,&
       '# colonne 12 : nbre inst. de part. perdues'             , &
       ' (reperage, avec poids)'                              ,/, &
       '# colonne 13 : % de part. perdues'                    ,/, &
       '# colonne 14 : nbre inst. de part. qui ont subi le'     , &
       ' clonage'                                             ,/, &
       '# colonne 15 : nbre inst. de part. qui ont subi le'     , &
       ' clonage (avec poids)'                                ,/, &
       '# colonne 16 : nbre inst. de nouvel. part. par clonage',/,&
       '# colonne 17 : nbre inst. de nouvel. part. par clonage'  ,&
       ' (avec poids)'                                        ,/, &
       '# colonne 18 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe '                                     ,/, &
       '# colonne 19 : nbre inst. de nouvel. part. eliminees par',&
       ' roulette russe (avec poids) '                        ,/, &
       '# ')

 1002 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     '      ,/, &
       '#                                              '      ,/, &
       '# colonne  1 : numero de pas de temps          '      ,/, &
       '# colonne  2 : temps physique                  '      ,/, &
       '# colonne  3 : nbre inst. de part.             '      ,/, &
       '# colonne  4 : nbre inst. de part. (avec poids)'      ,/, &
       '# colonne  5 : nbre inst. de part. injectees   '      ,/, &
       '# colonne  6 : nbre inst. de part. injectees   '        , &
       ' (avec poids)'                                        ,/, &
       '# colonne  7 : nbre inst. de part. sorties, ou deposees et supprimees',/,&
       '# colonne  8 : nbre inst. de part. sorties, ou deposees et supprimees',  &
       ' (avec poids)'                                        ,/, &
       '# colonne  9 : nbre inst. de part. deposees'          ,/, &
       '# colonne 10 : nbre inst. de part. deposees'            , &
       ' (avec poids)'                                        ,/, &
       '# colonne 11 : nbre inst. de part. perdues (reperage)',/, &
       '# colonne 12 : nbre inst. de part. perdues'             , &
       ' (reperage, avec poids)'                              ,/, &
       '# colonne 13 : % de part. perdues'                    ,/, &
       '# colonne 14 : nbre inst. de part. encrassees'          , &
       ' (Charbon) '                                          ,/, &
       '# colonne 15 : nbre inst. de part. encrassees'          , &
       ' (Charbon, avec poids)'                               ,/, &
       '# ')

 1003 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     '      ,/, &
       '#                                              '      ,/, &
       '# colonne  1 : numero de pas de temps          '      ,/, &
       '# colonne  2 : temps physique                  '      ,/, &
       '# colonne  3 : nbre inst. de part.             '      ,/, &
       '# colonne  4 : nbre inst. de part. (avec poids)'      ,/, &
       '# colonne  5 : nbre inst. de part. injectees   '      ,/, &
       '# colonne  6 : nbre inst. de part. injectees   '        , &
       ' (avec poids)'                                        ,/, &
       '# colonne  7 : nbre inst. de part. sorties, ou deposees et supprimees',/,&
       '# colonne  8 : nbre inst. de part. sorties, ou deposees et supprimees',  &
       ' (avec poids)'                                        ,/, &
       '# colonne  9 : nbre inst. de part. deposees'          ,/, &
       '# colonne 10 : nbre inst. de part. deposees'            , &
       ' (avec poids)'                                        ,/, &
       '# colonne 11 : nbre inst. de part. perdues (reperage) ',/,&
       '# colonne 12 : nbre inst. de part. perdues'             , &
        ' (reperage, avec poids)'                             ,/, &
       '# colonne 13 : % de part. perdues'                    ,/, &
       '# ')


 1004 format('# ** INFORMATIONS SUR LE CALCUL LAGRANGIEN     ',/, &
       '#    -------------------------------------     '      ,/, &
       '#                                              '      ,/, &
       '# colonne  1 : numero de pas de temps          '      ,/, &
       '# colonne  2 : temps physique                  '      ,/, &
       '# colonne  3 : nbre inst. de part.             '      ,/, &
       '# colonne  4 : nbre inst. de part. (avec poids)'      ,/, &
       '# colonne  5 : nbre inst. de part. injectees   '      ,/, &
       '# colonne  6 : nbre inst. de part. injectees   '        , &
       ' (avec poids)'                                        ,/, &
       '# colonne  7 : nbre inst. de part. sorties, ou deposees et supprimees',/,&
       '# colonne  8 : nbre inst. de part. sorties, ou deposees et supprimees',  &
       ' (avec poids)'                                        ,/, &
       '# colonne  9 : nbre inst. de part. deposees'          ,/, &
       '# colonne 10 : nbre inst. de part. deposees'            , &
       ' (avec poids)'                                        ,/, &
       '# colonne 11 : nbre inst. de part. resuspendues'      ,/, &
       '# colonne 12 : nbre inst. de part. resuspendues'        , &
       ' (avec poids)'                                        ,/, &
       '# colonne 13 : nbre inst. de part. perdues (reperage) ',/,&
       '# colonne 14 : nbre inst. de part. perdues'             , &
        ' (reperage, avec poids)'                             ,/, &
       '# colonne 15 : % de part. perdues'                    ,/, &
       '# ')

 2000 format(1x,i8,2x,e10.4,2x,5(i8,2x,e10.4),2x,e10.4,4(i8,2x,e10.4))
 2001 format(1x,i8,2x,e10.4,2x,5(i8,2x,e10.4),2x,e10.4,3(i8,2x,e10.4))
 2002 format(1x,i8,2x,e10.4,2x,5(i8,2x,e10.4),2x,e10.4,1(i8,2x,e10.4))
 2003 format(1x,i8,2x,e10.4,2x,5(i8,2x,e10.4),2x,e10.4)
 2004 format(1x,i8,2x,e10.4,2x,6(i8,2x,e10.4),2x,e10.4)

!====
! FIN
!====

return

end subroutine

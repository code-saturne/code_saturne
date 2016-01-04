!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine ellecd
!================
!===============================================================================
!  FONCTION  :
!  ---------

! LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
!            RELATIF AU MODULE ELECTRIQUE

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
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use ppincl
use elincl

!===============================================================================

implicit none

! Arguments

! Local variables

integer          it, ios , iesp, ii , i

!===============================================================================
!===============================================================================
! 0. INITIALISATION
!===============================================================================

ngazg = 1

!===============================================================================
! 1. LECTURE DU FICHIER DONNEES SPECIFIQUES
!===============================================================================

! --> Ouverture du fichier

if ( ippmod(ielarc).ge.1 .or. ippmod(ieljou).eq.3                 &
                         .or. ippmod(ieljou).eq.4 ) then
  open ( unit=impfpp, file=ficfpp,                                &
        STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',      &
                                        iostat=ios, err=99 )
  rewind ( unit=impfpp,err=99 )
endif

!==================================================
! 2. LECTURE D'UN FICHIER ARC ELECTRIQUE
!==================================================
!     Il n'y a pas de lecture de fichier prevue pour l'effet Joule



!  on transporte l'enthalpie, on suppose qu'il y a du rayonnement
!  la conductivite, la masse volumique, l'emissivite et
!  la viscosite dependent de la temperature

if ( ippmod(ielarc).ge.1 ) then

  ngazg = 0

! ----- NB de constituants et Nb de points de tabulation
!     on saute les lignes de commentaires
  do ii = 1, 7
    read (impfpp,*)
  enddo
  read ( impfpp,*,err=999,end=999 ) ngazg,npo

  if ( npo.gt.npot ) then
    write(nfecra,8000) npot
    call csexit (1)
    !==========
  endif

  if ( ngazg.gt. ngazgm) then
    write(nfecra,8001) ngazgm, ngazg
    call csexit (1)
    !==========
  endif

  if ( ngazg.lt. 1) then
    write(nfecra,8002)ngazg
    call csexit (1)
    !==========
  endif

! -----  Lecture de l'indicateur pour savoir ce que represente XKABEL
!      on saute les lignes de commentaires
  do ii = 1, 5
    read (impfpp,*)
  enddo
  read ( impfpp,*,err=999,end=999 ) ixkabe
  if ( ixkabe .lt. 0 .or. ixkabe .ge. 3 ) then
    write(nfecra,8003) ixkabe
    call csexit (1)
    !==========
  endif

! ----- En fonction de la temperature pour chaque espece courante
!          Enthalpie massique
!          Masse volumique
!          Chaleur massique
!          Conductivite electrique
!          Viscosite laminaire
!          Conductivite Thermique
!          Coefficent d'absorption (rayonnement)
!     on saute les lignes de commentaires au debut

  if(ngazg.gt.0.and.npo.gt.0) then
    do ii = 1, 7
      read (impfpp,*)
    enddo
    do iesp = 1, ngazg
      do it = 1, npo
        read (impfpp,*,err=999,end=999 )                          &
             th(it)         ,ehgazg(iesp,it),rhoel(iesp,it),      &
             cpel(iesp,it)  ,sigel(iesp,it) ,                     &
             visel(iesp,it) ,xlabel(iesp,it),                     &
             xkabel(iesp,it)
      enddo
    enddo
  endif

endif

!==================================================
! 3. LECTURE D'UN FICHIER EFFET JOULE
!==================================================

if ( ippmod(ieljou).eq. 3 .or. ippmod(ieljou).eq. 4 ) then

! ----- Lecture du transfo de reference

  read (impfpp,*,err=999,end=999 ) ntfref

! ----- Nombre de transfo

!       on saute 2 lignes de commentaires
  read (impfpp,*)
  read (impfpp,*)
  read (impfpp,*,err=999,end=999 ) nbtrf

!     Boucle sur le nombre de transfo

  do i=1,nbtrf

!         on saute la ligne de commentaire
    read (impfpp,*)

!         Tension primaire
    read (impfpp,*,err=999,end=999 ) tenspr(i)

!         Rapport du nombre de spires
    read (impfpp,*,err=999,end=999 ) rnbs(i)

!         Impedances complexes
    read (impfpp,*,err=999,end=999 ) zr(i),zi(i)

!         Type de branchement primaire
    read (impfpp,*,err=999,end=999 ) ibrpr(i)

!         Type de branchement secondaire
    read (impfpp,*,err=999,end=999 ) ibrsec(i)

  enddo

! ----- Nombre d'electrodes

  read (impfpp,*)
  read (impfpp,*)
  read (impfpp,*,err=999,end=999 ) nbelec

!       Boucle sur le nombre d'electrodes

  do i=1,nbelec

!         Tension primaire
    read (impfpp,*,err=999,end=999 ) ielecc(i),ielect(i),         &
                                     ielecb(i)

  enddo

endif

!==================================================
! 4. LECTURE D'UN FICHIER MIGRATION IONIQUE
!==================================================

!    c'est plus complique (chaque espece peut avoir sa propre mobilite)
!    mais il n'y a pas forcement de rayonnement, ni meme de
!    chauffage significatif par effet Joule

!==============================================
! 5. CALCULS DE DONNEES COMPLEMENTAIRES
!==============================================


return


!============================
! 3. SORTIE EN ERREUR
!============================

  99  continue
write ( nfecra,9998 )
call csexit (1)
!==========

  999 continue
write ( nfecra,9999 )
call csexit (1)
!==========


!--------
! FORMATS
!--------


 8000 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,&
'@                                                            ',/,&
'@  Le nombre de points de tabulation lu dans le fichier de   ',/,&
'@    donnees doit etre un entier inferieur ou egal           ',/,&
'@              a NPOT   = ',I10                               ,/,&
'@    Il vaut ici NPO    = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8001 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes courantes lu dans le fichier de      ',/,&
'@    doit etre un entier inferieur ou egal                   ',/,&
'@              a NGAZGM = ',I10                               ,/,&
'@    Il vaut ici NGAZG  = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8002 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes courantes lu dans le fichier de      ',/,&
'@    doit etre un entier superieur ou egal a 1.              ',/,&
'@    Il vaut ici NGAZG = ',I10                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8003 format (                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,&
'@                                                            ',/,&
'@  La valeur de l''indicateur pour le rayonnement            ',/,&
'@    doit etre comprise entre 0 et 2                         ',/,&
'@    elle vaut ici IXKABE = ',I10                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9998 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,&
'@                                                            ',/,&
'@  Erreur a l''ouverture du fichier parametrique.            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (ELLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (VERSIONS ELECTRIQUES)          ',/,&
'@                                                            ',/,&
'@  Erreur a la lecture du fichier parametrique.              ',/,&
'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,&
'@    ou son format inadapte.                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine



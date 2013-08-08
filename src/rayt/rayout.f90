!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine rayout &
!================

 ( rtpa   , rtp    , propce , propfb )

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
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use pointe
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)

! Local variables

character        rubriq*64
character        cphase*2
character        ficsui*32
integer          itrav1 , ip
integer          ierror , nberro , irtyp , itysup , nbval
integer          ivers  , ilecec
integer          impavr

!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

!===============================================================================
! 1. ECRITURE DU FICHIER SUITE DU MODULE DE RAYONNEMENT
!===============================================================================


! ---> Ouverture (et on saute si erreur)
!     ILECEC = 2 : ecriture

write(nfecra,6010)

ilecec = 2
ficsui = 'radiative_transfer'
call opnsui(ficsui, len(ficsui), ilecec, impavr, ierror)
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

  RUBRIQ = 'luminance'
  call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
              propce(1,ipproc(ilumin)),ierror)
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
  write(nfecra,8011) ficsui
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

!----
! FIN
!----

end subroutine

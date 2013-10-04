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

 ( propce )

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
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use field
!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)
double precision, dimension(:), pointer :: f_val
! Local variables

character        rubriq*64
character        ficsui*32
integer          ierror , irtyp , itysup , nbval
integer          ivers  , ilecec
integer          impavr

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

ivers  = 111
itysup = 0
nbval  = 1
irtyp  = 1
rubriq = 'version_fichier_suite_rayonnement'
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,ivers)

itysup = 0
nbval  = 1
irtyp  = 1

write(nfecra,6012)

! Temps (par securite)

rubriq = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
irtyp  = 1
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,ntcabs)

rubriq = 'instant_precedent'
itysup = 0
nbval  = 1
irtyp  = 2
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,ttcabs)

! Donnees

!     Aux faces de bord

itysup = 3
nbval  = 1
irtyp  = 2

rubriq = 'tparoi_fb'
call field_get_val_s(itparo, f_val)
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,f_val)

rubriq = 'qincid_fb'
call field_get_val_s(iqinci, f_val)
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,f_val)

rubriq = 'hfconv_fb'
call field_get_val_s(ihconv, f_val)
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,f_val)

rubriq = 'flconv_fb'
call field_get_val_s(ifconv, f_val)
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,f_val)

!     Aux cellules

itysup = 1
nbval  = 1
irtyp  = 2

rubriq = 'rayimp_ce'
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
            propce(1,ipproc(itsri(1))))

rubriq = 'rayexp_ce'
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
            propce(1,ipproc(itsre(1))))

rubriq = 'luminance'
call ecrsui(impavr,rubriq,len(rubriq),itysup,nbval,irtyp,       &
            propce(1,ipproc(ilumin)))

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
! Formats
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
! End
!----

end subroutine

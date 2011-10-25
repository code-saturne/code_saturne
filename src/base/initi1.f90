!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine initi1 &
!================

 ( iverif )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES COMMONS

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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
use optcal
use entsor
use ihmpre

!===============================================================================

implicit none

! Arguments

integer          iverif

! Local variables

integer          iok

!===============================================================================

!===============================================================================
! 1. INITIALISATION DES VARIABLES EN COMMON
!      AVANT INTERVENTION UTILISATEUR.
!      . mots-cles dynamiques
!      . entsor    .h
!      . dimens    .h
!      . numvar    .h
!      . pointe    .h
!      . optcal    .h
!      . mltgrd    .h
!      . cstphy    .h
!===============================================================================

call iniini
!==========

!===============================================================================
! 2. ENTREE DES DONNEES PAR L'UTILISATEUR
!      ET POSITIONNEMENT DES VARIABLES (VARPOS)
!===============================================================================

call iniusi(iverif)
!==========

call ppini1
!==========

call rayopt
!==========

call lagopt
!==========

! En mode verification, on positionne IMRGRA a 2 de maniere a creer
! le voisinage etendu (complet) necessaire a certains modes de calcul
! de gradient
if (iverif.eq.1) imrgra = 2


!===============================================================================
! 3. DEFINITION DES COUPLAGES AVEC SYRTHES
!===============================================================================

! Le nombre de couplage SYRTHES doit etre connu avant MODINI a des fins
! de verification de coherence avec la definition des scalaires

if (iihmpr.eq.1) then
  call uisyrc
  !==========
endif

call ussyrc
!==========

call ussatc
!==========

!===============================================================================
! 4. MODIFS APRES USINI1
!===============================================================================

call modini
!==========

!===============================================================================
! 5. VERIFS APRES USINI1
!===============================================================================

iok = 0

call verini (iok)
!==========

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
else
  write(nfecra,9998)
endif

!===============================================================================
! 6. Initial definition of fields
!===============================================================================

call fldini
!==========


#if defined(_CS_LANG_FR)

 9998 format(                                                           &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                               (interface, usini1 et autres).',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier les donnees entrees dans l''interface, usini1 ou ',/,&
'@    les autres sous-programmes d''initialisation.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9998 format(                                                           &
''                                                             ,/,&
' No error detected during the data verification'              ,/,&
'                              (interface, usini1 and others).',/)
 9999 format(                                                           &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION'                ,/,&
'@    ========'                                                ,/,&
'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLET'  ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run (',I10,' errors).'         ,/,&
'@'                                                            ,/,&
'@  See previous impressions for more informations.'           ,/,&
'@  Verify the provided data in the interface, usini1 or'      ,/,&
'@    the other initialization subroutines.'                   ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

!===============================================================================
! 6. IMPRESSIONS
!===============================================================================

if (iverif.eq.1) return

call impini
!==========

return
end subroutine

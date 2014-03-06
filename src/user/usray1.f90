!-------------------------------------------------------------------------------

!VERS

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

subroutine usray1
!================



!===============================================================================
!  FONCTION  :
!  ---------


!   ============================================================
!   ROUTINE UTILISATEUR : INITIALISATION DES COMMONS RAYONNEMENT
!   ============================================================



!  CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS

!  METHODE DES VOLUMES FINIS.


! 1/ DONNEES DES LUMINANCES ENTRANTES AUX LIMITES DU DOMAINE
!       (C    .L : REFLEXION ET EMISSION ISOTROPE)

!                              ->  ->           ->
! 2/ CALCUL DE LA LUMINANCE L( X , S ) AU POINT X

!                                   D L
!    PAR RESOLUTION DE L'EQUATION : --- = -TK.L +TS
!                                   D S
!                       ->                o
!    OU ENCORE : DIV (L.S ) = -TK.L + TK.L

!                                 ->   /    ->  ->  ->
! 3/ CALCUL DES DENSITES DE FLUX  Q = /  L( X , S ).S DOMEGA
!                                    /4.PI

!                                      /    ->  ->
!        ET DE L'ABSORPTION       SA= /  L( X , S ).  DOMEGA
!                                    /4.PI

!    PAR INTEGRATION DES LUMINANCES SUR LES ANGLES SOLIDES.

!    N . B : CA SERT A CALCULER LE TAUX D'ECHAUFFEMENT
!    -----
!                                      /    ->  ->  ->  ->
! 4/ CALCUL DU FLUX INCIDENT QINCID = /  L( X , S ).S . N DOMEGA
!                                    /->->
!       ->                          / S.N >0
!       N NORMALE FLUIDE VERS PAROI




!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use ppppar
use radiat

!===============================================================================

implicit none

integer          ipp

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================



! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


ipass = ipass + 1

!===============================================================================
! 1. UTILISATION DU MODULE DE TRANSFERTS RADIATIFS
!===============================================================================

if (ipass.eq.1 .or. ipass.eq.2) then

  !-->  IIRAYO = 0 : PAS DE TRANSFERTS RADIATIFS (PAR DEFAUT)
  !            = 1 : TRANSFERTS RADIATIFS, METHODE DES ORDONNEES DISCRETES
  !            = 2 : TRANSFERTS RADIATIFS, APPROXIMATION P-1

  iirayo = 1

endif

!===============================================================================
! 2. PARAMETRES DU MODULE DE TRANSFERTS RADIATIFS
!===============================================================================

if (ipass.eq.2) then

  !-->  INDICATEUR SUITE DE CALCUL (LECTURE DU FICHIER SUITE DE RAYONNEMENT)
  !     (0      : PAS DE LECTURE D'UN FICHIER SUITE DE RAYONNEMENT
  !      1      : RELECTURE D'UN FICHIER SUITE DE RAYONNEMENT
  !      ISUITE : RELECTURE D'UN FICHIER SUITE DE RAYONNEMENT SI LE CALCUL FLUIDE EST
  !               AUSSI UNE SUITE )

  isuird = isuite

  !-->  FREQUENCE DE PASSAGE DANS LE MODULE DE RAYONNEMENT

  nfreqr = 1

  !-->  Quadrature Sn (n(n+2) directions)
  !
  ! 1 : S4 (24 directions)
  ! 2 : S6 (48 directions)
  ! 3 : S8 (80 directions)
  !
  !-->  Quadrature Tn (8n^2 directions)
  !
  ! 4 : T2 (32 directions)
  ! 5 : T2 (128 directions)
  ! 6 : Tn (8*ndirec^2 directions)

  i_quadrature = 4

  ndirec =  3

  !-->  INITIALISATION DU MODE DE CALCUL DU TERME SOURCE RADIATIF EXPLICITE
  !     IDIVER = 0 => CALCUL SEMI-ANALYTIQUE
  !     IDIVER = 1 => CALCUL CONSERVATIF
  !     IDIVER = 2 => CALCUL SEMI-ANALYTIQUE CORRIGE POUR ETRE CONSERVATIF
  !     (EN RAYONNEMENT TRANSPARENT, LE CHOIX EST SANS INFLUENCE)

  idiver = 2

  !--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS TEMPERATURE DE PAROI

  iimpar = 1

  !--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS SOLVEUR

  iimlum = 0

  !--> SI COMBUSTION GAZ OU CHARBON : CALCUL AUTOMATIQUE
  !    DU COEFFICIENT D'ABSORPTION

  !    IMODAK = 0 : PAS DE CALCUL AUTOMATIQUE (DEFAUT)
  !           = 1 : CALCUL AUTOMATIQUE (MODELE DE MODAK)

  imodak = 0

  !nband = 1

  !ngauss = 1

endif

!===============================================================================
! 3. RENSEIGNEMENTS POUR LE POST-PROCESSING
!===============================================================================

if (ipass.eq.3) then

  ! ihisvr( ) = sortie historique (nombre de sondes et numeros)
  ! si ihisvr(    .,1)  = -1 sortie sur toutes les sondes

  !--> LUMINENCE

  ipp = ipppro(ipproc(ilumin))
  ihisvr(ipp,1) = -1

  !--> VECTEUR DENSITE DE FLUX RADIATIF

  !       composante x
  ipp = ipppro(ipproc(iqx))
  ihisvr(ipp,1) = -1

  !       composante y
  ipp = ipppro(ipproc(iqy))
  ihisvr(ipp,1) = -1

  !       composante z
  ipp = ipppro(ipproc(iqz))
  ihisvr(ipp,1) = -1

endif

!===============================================================================

return

end subroutine usray1

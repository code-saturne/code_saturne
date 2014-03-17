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

subroutine elini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!   INIT DES OPTIONS DES VARIABLES POUR LE MODULE ELECTRIQUE
!      EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1

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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl
use elincl
use ihmpre
use mesh

!===============================================================================

implicit none

! Local variables

integer          idimve
integer          iok
integer          isc , ivar

!===============================================================================

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! --> Conduction ionique (a developper)

! 1.4 Donnees physiques ou numeriques propres aux scalaires ELECTRIQUES
! =====================================================================

! --> Conditions associees aux potentiels
!     (les autres variables ont des comportements par defaut)
ivar = isca(ipotr)
iconv (ivar) = 0
istat (ivar) = 0
idiff (ivar) = 1
idifft(ivar) = 0
idircl(ivar) = 1
imgr  (ivar) = 0

if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ivar = isca(ipoti)
  iconv (ivar) = 0
  istat (ivar) = 0
  idiff (ivar) = 1
  idifft(ivar) = 0
  idircl(ivar) = 1
  imgr  (ivar) = 0
endif

if (ippmod(ielarc).ge.2) then
  do idimve = 1, ndimve
    ivar = isca(ipotva(idimve))
    iconv (ivar) = 0
    istat (ivar) = 0
    idiff (ivar) = 1
    idifft(ivar) = 0
    idircl(ivar) = 1
    imgr  (ivar) = 0
  enddo
endif

! --> "Viscosite" associee au potentiel vecteur
!     (c'est la seule qui est constante)
if (ippmod(ielarc).ge.2) then
  visls0(ipotva(1)) = 1.d0
  visls0(ipotva(2)) = 1.d0
  visls0(ipotva(3)) = 1.d0
endif

! --> Schmidt ou Prandtl turbulent
!     (pour les potentiels, c'est inutile puisque IDIFFT=0)

do isc = 1, nscapp
  sigmas(iscapp(isc)) = 0.7d0
enddo

! ---> Pour tous les scalaires

do isc = 1, nscapp

! ----- Niveau de detail des impressions pour les variables et
!          donc les scalaires (valeurs 0 ou 1)
!          Si = -10000 non modifie par l'utilisateur -> niveau 1

  ivar = isca(iscapp(isc))
  if (iwarni(ivar).eq.-10000) then
    iwarni(ivar) = 1
  endif

! ----- Informations relatives a la resolution des scalaires

!       - Facteur multiplicatif du pas de temps

  cdtvar(ivar) = 1.d0

!         - Schema convectif % schema 2ieme ordre
!           = 0 : upwind
!           = 1 : second ordre
  blencv(ivar) = 1.d0

!         - Type de schema convectif second ordre (utile si BLENCV > 0)
!           = 0 : Second Order Linear Upwind
!           = 1 : Centre
  ischcv(ivar) = 1

!         - Test de pente pour basculer d'un schema centre vers l'upwind
!           = 0 : utlisation automatique du test de pente
!           = 1 : calcul sans test de pente
  isstpc(ivar) = 0

!         - Reconstruction des flux de convection et de diffusion aux faces
!           = 0 : pas de reconstruction
  ircflu(ivar) = 1

enddo

!===============================================================================
! 2. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! --> Coefficient de relaxation de la masse volumique
!      a partir du 2ieme pas de temps, on prend :
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = 0.d0

! --> Recalage des variables electriques
!      IELCOR = 0 : pas de correction
!      IELCOR = 1 : correction
ielcor = 0

!     Intensite de courant imposee (arc electrique) ou
!                Puissance imposee (Joule)
couimp = 0.d0
puisim = 0.d0

!     Differentiel de potentiel Initial en arc (et Joule)
dpot = 0.d0

!     Coefficient pour la correction en Joule
coejou = 1.d0

! ---> Masse volumique variable et viscosite variable (pour les suites)
irovar = 1
ivivar = 1

! ---> Modele pour le recalage de l'intensite (arc electrique)
!       MODREC = 1 : modele standard
!       MODREC = 2 : modele avec un plan de recalage
modrec = 1

!===============================================================================
! 3. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

if (iihmpr.eq.1) then
    call uicpi1 (srrom, diftl0)
    ! gas number and radiatif transfer are read in dp_ELE
    call uieli1 (ippmod(ieljou), ippmod(ielarc), ielcor, couimp, puisim, &
                 modrec, idreca, crit_reca)

    ! Initial value for dpot is set to 1000 V.
    dpot = 1000.d0

endif

call useli1(iihmpr)
!==========

!===============================================================================
! 4. VERIFICATION DES DONNEES ELECTRIQUES
!===============================================================================

iok = 0

call elveri (iok)
!==========

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
  !==========
else
  write(nfecra,9998)
endif

 9998 format(                                                     &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (useli1).',/)
 9999 format(                                                     &
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
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine

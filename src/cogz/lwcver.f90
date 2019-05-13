!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine lwcver &
!================

 ( iok    )

!===============================================================================
!  FONCTION  :
!  ---------

! VERIFICATION DES PARAMETRES DE CALCUL
!   COMBUSTION GAZ : FLAMME DE PREMELANGE - MODELE LWC
!     APRES INTERVENTION UTILISATEUR
!       (COMMONS)
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
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

!===============================================================================
!===============================================================================
! 1. OPTIONS DU CALCUL : TABLEAUX DE ppincl.h : formats 2000
!===============================================================================

! --> Coefficient de relaxation de la masse volumique

if( srrom.lt.0d0 .or. srrom.ge.1d0) then
  WRITE(NFECRA,2000)'SRROM ', SRROM
  iok = iok + 1
endif

!===============================================================================
! 2. TABLEAUX DE cstphy.h et ppthch.F : formats 3000
!===============================================================================

! --> Masse volumique

if( ro0.lt.0d0) then
  WRITE(NFECRA,3000)'RO0   ', RO0
  iok = iok + 1
endif

! --> Diffusivite dynamique en kg/(m s) : DIFTL0

if( diftl0.lt.0d0) then
  WRITE(NFECRA,3010)'DIFTL0', DIFTL0
  iok = iok + 1
else
  visls0(iscalt) = diftl0
endif

! --> Constante du modele LWC

if( vref.lt.0d0) then
  WRITE(NFECRA,3020)'VREF', VREF
  iok = iok + 1
endif
if( lref.lt.0d0) then
  WRITE(NFECRA,3020)'LREF', LREF
  iok = iok + 1
endif
if( ta.lt.0d0) then
  WRITE(NFECRA,3020)'TA', TA
  iok = iok + 1
endif
if( tstar.lt.0d0) then
  WRITE(NFECRA,3020)'TSTAR', TSTAR
  iok = iok + 1
endif

!===============================================================================
! 3. FORMATS VERIFICATION
!===============================================================================
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,                            ' DOIT ETRE UN REEL    ',/,&
'@    SUPERIEUR OU EGAL A ZERO ET INFERIEUR STRICTEMENT A 1   ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uslwc1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uslwc1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uslwc1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A4,' DOIT ETRE UN REEL POSITIF                        ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uslwc1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 6. SORTIE
!===============================================================================

return
end subroutine

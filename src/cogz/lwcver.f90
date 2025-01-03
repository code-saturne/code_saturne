!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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
use ppincl

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

!===============================================================================
! 2. TABLEAUX DE cstphy.h et ppthch.F : formats 3000
!===============================================================================

! --> Constante du modele LWC

if (vref.lt.0d0) then
  write(nfecra,3020)'vref', vref
  iok = iok + 1
endif
if (lref.lt.0d0) then
  write(nfecra,3020)'lref', lref
  iok = iok + 1
endif
if (ta.lt.0d0) then
  write(nfecra,3020)'ta', ta
  iok = iok + 1
endif
if (tstar.lt.0d0) then
  write(nfecra,3020)'tstar', tstar
  iok = iok + 1
endif

!===============================================================================
! 3. FORMATS VERIFICATION
!===============================================================================

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

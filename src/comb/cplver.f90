!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine cplver &
!================

 ( iok    )

!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      VERIFICATION DES PARAMETRES DE CALCUL
!        APRES INTERVENTION UTILISATEUR
!        (COMMONS)
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

if( srrom.lt.0d0 .or. srrom.gt.1d0) then
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

!===============================================================================
! 3. FORMATS VERIFICATION
!===============================================================================

 2000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP WHILE DEFINING INPUT DATA'                 ,/,&
'@    ========'                                                ,/,&
'@    SPECIFIC PHYSICS (PULVERIZED COAL)'                      ,/,&
'@'                                                            ,/,&
'@    ', A6, ' MUST BE A REAL BETWEEN 0 AND 1'                 ,/,&
'@    ITS VALUE HERE IS ',E14.5                                ,/,&
'@'                                                            ,/,&
'@  The calculation can NOT be run.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier user_coal_ini_1.'                                 ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 3000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP WHILE DEFINING INPUT DATA'                 ,/,&
'@    ========'                                                ,/,&
'@    SPECIFIC PHYSICS (PULVERIZED COAL)'                      ,/,&
'@'                                                            ,/,&
'@    ', A6, ' MUST BE A POSITIVE REAL'                        ,/,&
'@    ITS VALUE HERE IS ',E14.5                                ,/,&
'@'                                                            ,/,&
'@  The calculation can NOT be run.'                           ,/,&
'@'                                                            ,/,&
'@  Check user_coal_ini_1.'                                    ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP WHILE DEFINING INPUT DATA'                 ,/,&
'@    ========                                                ',/,&
'@    SPECIFIC PHYSICS (PULVERIZED COAL)'                     ,/,&
'@'                                                            ,/,&
'@    ',A6,' MUST BE A POSITIVE REAL'                          ,/,&
'@    ITS VALUE HERE IS ',E14.5                                ,/,&
'@'                                                            ,/,&
'@  The calculation can run.'                                  ,/,&
'@'                                                            ,/,&
'@  Check user_coal_ini_1.'                                    ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

!===============================================================================
! 4. SORTIE
!===============================================================================

return
end subroutine

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine ppini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES OPTIONS DES VARIABLES SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
!   EN COMPLEMENT DE CE QUI A ETTE DEJA FAIT DANS USINI1

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
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ihmpre

!===============================================================================

implicit none

! Local variables

integer :: ii

!===============================================================================

!===============================================================================
! 0. VERIFICATION ISCACP
!===============================================================================

! L'utilisateur ne doit pas y avoir touche.

do ii = 1, nscapp
  if(iscacp(iscapp(ii)).ne.-10) then
    write(nfecra,1001) ii, iscapp(ii), iscapp(ii), iscacp(iscapp(ii))
    call csexit (1)
    !==========
  endif
enddo

if (itherm .eq. 1) then
  iscacp(iscalt) = 1
endif

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! ---> Combustion gaz : Flamme de diffusion  (Chimie 3 points)
!                       Flamme de premelange (Modele EBU)
!                       Flamme de premelange (Modele LWC)

if ( ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                 &
                         .or. ippmod(icolwc).ge.0 ) then
  call coini1
endif

! ---> Physique particuliere : Combustion Charbon Pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_param
endif

! ---> Physique particuliere : Combustion Eulerienne Charbon Pulverise
!                              Couplee Transport Lagrangien

if ( ippmod(icpl3c).ge.0 ) then
  call cplin1
endif

! ---> Physique particuliere : Combustion fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_param
endif

! ---> Physique particuliere : Compressible

if ( ippmod(icompf).ge.0) then
  call cfini1
endif

! ---> Physique particuliere : Versions electriques

if (ippmod(ieljou).ge.1.or.ippmod(ielarc).ge.1) then
  call elini1 (visls0, diftl0, idircl, isca)
endif

! ---> Physique particuliere : Ecoulements atmospheriques

if ( ippmod(iatmos).ge.0 ) then
  call atini1
endif

! ---> Physique particuliere : Aerorefrigerants

if ( ippmod(iaeros).ge.0) then
  call ctini1
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1001 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@'                                                            ,/,&
'@  Les valeurs de ISCACP pour les scalaires de modele'        ,/,&
'@  (i.e. non_utilisateur) sont renseignees automatiquement.'  ,/,&
'@'                                                            ,/,&
'@  L''utilisateur ne doit pas les renseigner, or'             ,/,&
'@    pour le scalaire ', i10,' correspondant au scalaire'     ,/,&
'@    de modele ', i10,' on a'                                 ,/,&
'@    ISCACP(', i10,') = ', i10                                ,/,&
'@'                                                            ,/,&
'@  Le calcul ne sera pas execute.'                            ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres.'                                  ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 1001 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP WHILE READING INPUT DATA'                  ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The values of ISCACP are set automatically for model'      ,/,&
'@  (i.e. non-user) scalars.'                                  ,/,&
'@'                                                            ,/,&
'@  The user should not set a value for them, however'         ,/,&
'@    for the scalar ', i10,' corresponding to the model'      ,/,&
'@    scalar ', i10,' we have'                                 ,/,&
'@    iscacp(', i10,') = ', i10                                ,/,&
'@'                                                            ,/,&
'@  The calculation could NOT run.'                            ,/,&
'@'                                                            ,/,&
'@  Check parameters.'                                         ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

return
end subroutine

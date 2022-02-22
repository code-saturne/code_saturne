!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!   EN COMPLEMENT DE CE QUI A ETE DEJA FAIT DANS USINI1

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
use field

!===============================================================================

implicit none

! Local variables

integer :: ii, iscacp

!===============================================================================

!===============================================================================
! Check "is_temperature"
!===============================================================================

! L'utilisateur ne doit pas y avoir touche.

do ii = 1, nscapp
  call field_get_key_int(ivarfl(isca(iscapp(ii))), kscacp, iscacp)
  if (iscacp.ne.-1) then
    write(nfecra,1001) iscapp(ii), iscapp(ii), iscacp
    call csexit(1)
  endif
enddo

if (itherm .eq. 1) then
  call field_set_key_int(ivarfl(isca(iscalt)), kscacp, 1)
endif

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! ---> Combustion gaz : Flamme de diffusion  (Chimie 3 points)
!                       Flamme de diffusion  (Steady laminar flamelet)
!                       Flamme de premelange (Modele EBU)
!                       Flamme de premelange (Modele LWC)

if ( ippmod(icod3p).ge.0 .or. ippmod(islfm ).ge.0 .or.            &
     ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0 ) then
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
  call elini1(diftl0)
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

 1001 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: STOP WHILE READING INPUT DATA'                  ,/,&
'@    ======='                                                 ,/,&
'@'                                                            ,/,&
'@  The values of "is_temperature" are set automatically for'  ,/,&
'@  model (i.e. non-user) scalars.'                            ,/,&
'@'                                                            ,/,&
'@  The user should not set a value for them, however'         ,/,&
'@    for the scalar ', i10,' corresponding to the model'      ,/,&
'@    scalar ', i10,' we have'                                 ,/,&
'@    "is_temperature" = ', i11                                ,/,&
'@'                                                            ,/,&
'@  The calculation could NOT run.'                            ,/,&
'@'                                                            ,/,&
'@  Check parameters.'                                         ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

return
end subroutine

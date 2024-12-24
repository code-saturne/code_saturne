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
!> \file atini1.f90
!> \brief Initialisation of variable options for the atmospheric module in
!>        before what is done in usipsu/cs_user_parameters functions
!>

subroutine atini1 () &
 bind(C, name='cs_f_atini1')

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
use ppincl
use atincl
use atsoil
use atchem
use atimbr
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          ii
double precision turb_schmidt

!===============================================================================

!===============================================================================
! 1. VERIFICATIONS
!===============================================================================

if (ippmod(iatmos).le.1) then
  if (iatra1.eq.1.or.iatsoil.ge.1) then
    write(nfecra, 1003)
    call csexit(1)
  endif
endif

!===============================================================================
! 2. Transported variables for IPPMOD(IATMOS) = 0, 1 or 2
!===============================================================================

if (ippmod(iatmos).eq.0) then

  ! constant density
  irovar = 0

else if (ippmod(iatmos).ge.1) then

  ! for the dry or humid atmosphere case, non constant density
  irovar = 1

endif

!===============================================================================
! 5. Turbulent Schmidt and Prandtl number for atmospheric flows
!===============================================================================

if (nscal.gt.0) then
  do ii = 1, nscal
    turb_schmidt = 0.7d0
    call field_set_key_double(ivarfl(isca(ii)), ksigmas, turb_schmidt)
  enddo
endif

!===============================================================================
! 6. Force Rij Matrix stabilisation for all atmospheric models
!===============================================================================

if (itytur.eq.3.and.irijnu.eq.0) irijnu = 1

!===============================================================================
! 7. Some allocation and mapping for meteo...
!===============================================================================

if (ippmod(iatmos).ge.0) then

  call allocate_map_atmo

  if (ichemistry.ge.1) then
    call init_chemistry
  endif

endif

!--------
! Formats
!--------

 1003 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA'               ,/,&
'@    ========='                                               ,/,&
'@                ATMOSPHERIC  MODULE'                         ,/,&
'@'                                                            ,/,&
'@  Ground model (soil_model)'                                 ,/,&
'@   and radiative model (radiative_model_1d)'                 ,/,&
'@   are only available with humid atmosphere model'           ,/,&
'@   or dry atmosphere model.'                                 ,/,&
'@  Computation CAN NOT run.'                                  ,/,&
'@'                                                            ,/,&
'@  Check the input data given through the User Interface'     ,/,&
'@   or in cs_user_model.'                                     ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

!----
! End
!----

return
end subroutine atini1

!> \brief Finalize initialisation of variable options for the atmospheric module
!>       after usipsu/cs_user_parameters functions

subroutine atini2

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
use ppincl
use atincl
use atsoil
use atchem
use atimbr
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

!===============================================================================

!===============================================================================
! 7. Some initialization for meteo...
!===============================================================================

if (ippmod(iatmos).ge.0) then

  call init_meteo

  if (imbrication_flag) then
    call activate_imbrication
  endif

  call cs_at_data_assim_build_ops

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine atini2

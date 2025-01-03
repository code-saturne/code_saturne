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

subroutine impini () &
  bind(C, name='cs_f_impini')

!===============================================================================
! Purpose:
!  ---------

! Print computation parameters after user changes in cs_user_parameters.cpp

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use entsor
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================
! 1. Introduction
!===============================================================================

write(nfecra,1000)

 1000 format(                                                     &
                                                                /,&
' ===========================================================', /,&
                                                                /,&
'               CALCULATION PARAMETERS SUMMARY',                /,&
'               ==============================',                /,&
                                                                /,&
' -----------------------------------------------------------', /)

 9900 format(                                                     &
                                                                /,&
' -----------------------------------------------------------', /)

!===============================================================================
! DEFINITION GENERALE DU CAS
!===============================================================================

! --- Dimensions

write(nfecra,1500)
write(nfecra,1520) nvar,nscal,nscaus,nscapp

write(nfecra,9900)

 1500 format(                                                     &
                                                                /,&
' ** DIMENSIONS',                                               /,&
'    ----------',                                               /)
 1520 format(                                                     &
' --- Physics',                                                 /,&
'       NVAR   = ',4x,i10,    ' (Nb variables                )',/,&
'       NSCAL  = ',4x,i10,    ' (Nb scalars                  )',/,&
'       NSCAUS = ',4x,i10,    ' (Nb user scalars             )',/,&
'       NSCAPP = ',4x,i10,    ' (Nb specific physics scalars )',/)

return
end subroutine

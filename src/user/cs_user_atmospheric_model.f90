!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_atmospheric_model.f90
!>
!> \brief User subroutines dedicated to the atmospheric model.
!>
!> See \ref cs_user_atmospheric_model for examples.
!-------------------------------------------------------------------------------

!===============================================================================

!===============================================================================

!> \brief Fill in vertical profiles of atmospheric properties prior to solve
!> 1D radiative transfers. Altitudes (\ref zvert array) are defined in
!> \ref usatd.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in,out] preray        pressure vertical profile
!> \param[in,out] temray        real temperature vertical profile
!> \param[in,out] romray        density vertical profile
!> \param[in,out] qvray         water vapor content vertical profile
!> \param[in,out] qlray         water liquid content vertical profile
!> \param[in,out] ncray         droplets density vertical profile
!> \param[in,out] aeroso        aerosol concentration vertical profile
!______________________________________________________________________________!

subroutine cs_user_atmo_1d_rad_prf &
     ( preray, temray, romray, qvray, qlray, ncray, aeroso )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use mesh

!===============================================================================

implicit none

! Arguments

double precision preray(kmx), temray(kmx), romray(kmx), qvray(kmx)
double precision qlray(kmx), ncray(kmx), aeroso(kmx)

! Local variables

!===============================================================================

return
end subroutine cs_user_atmo_1d_rad_prf

!===============================================================================

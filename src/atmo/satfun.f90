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
!> \file satfun.f90
!> \brief Computes the saturation mixing ratio (kg/kg) of water
!>      in the atmosphere.

!> \brief Computes the saturation mixing ratio (kg/kg) of water
!>      in the atmosphere.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]   t   thermodynamic temperature of the air parcel in Kelvin
!> \param[in]   p   pressure of the air parcel in Pascal
!-------------------------------------------------------------------------------
double precision function qsatliq (t,p)
!-------------------------------------------------------------------------------

use paramx ! needed by cstphy
use ppppar ! needed by atincl

use cstphy ! defines tkelvi
use atincl ! defines rvsra

!-------------------------------------------------------------------------------

implicit none

double precision t,p
double precision esat

!***********************************************************************

esat    = 610.78d0*exp(17.269d0*(t - tkelvi)/(t - 35.86d0))
qsatliq = esat/(rvsra*p + esat*(1d0 - rvsra))

end function qsatliq

!-------------------------------------------------------------------------------
!> \brief computes the saturation water vapour pressure function of the
!>      temperature (K)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]   t   thermodynamic temperature of the air parcel in Kelvin
!-------------------------------------------------------------------------------
double precision function esatliq (t)

use paramx !needed by cstphy
use ppppar !needed by atincl

use cstphy !defines tkelvi
use atincl !defines rvsra

!-------------------------------------------------------------------------------

implicit none
double precision t

!***********************************************************************

esatliq = 610.78d0*exp(17.269d0*(t - tkelvi)/(t - 35.86d0))

end function esatliq

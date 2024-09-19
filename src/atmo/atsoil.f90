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

!> \file atsoil.f90
!> \brief Module for the atmospheric soil model adapted to the IGN "land use"
!>      file format

module atsoil

use, intrinsic :: iso_c_binding

!> \defgroup at_soil
!=============================================================================
!> \addtogroup at_soil
!> \{
implicit none
!> Number of boundary faces with soil features
integer, save :: nfmodsol
!> Number of soil model's type
integer, save :: nbrsol


! Initialisation values for soil variables
! (filled in usispu/cs_user_parameters)

!> initial soil surface temperature
!> for Sea, it is also the surface temperature
real(c_double), pointer, save :: tsini
!> initial deep soil temperature
real(c_double), pointer, save :: tprini
!> initial soil specific humidity
real(c_double), pointer, save :: qvsini
!> initial water content of the first reservoir
real(c_double), pointer, save :: w1ini
!> initial water content of the second reservoir
real(c_double), pointer, save :: w2ini

! Array for soil categories

!> Thermal inertia of the soil
double precision, dimension(:), pointer :: csol
!> Dynamic roughness length
double precision, dimension(:), pointer :: rugdyn
!> Thermal roughness length
double precision, dimension(:), pointer :: rugthe
!> Albedo
double precision, dimension(:), pointer :: soil_cat_albedo
!> emissivity
double precision, dimension(:), pointer :: soil_cat_emissi
!> Vegetation index
double precision, dimension(:), pointer :: soil_cat_vegeta
!> maximum water capacity of shallow reservoir
double precision, dimension(:), pointer :: soil_cat_c1w
!> ratio of the maximum water capacity of the shallow reservoir to the deep
!> reservoir [0,1]
double precision, dimension(:), pointer :: soil_cat_c2w
!> Rij value for Rij1
double precision, dimension(:), pointer :: soil_cat_r1
!> Rij value for Rij2
double precision, dimension(:), pointer :: soil_cat_r2

!> \}
!=============================================================================

end module atsoil

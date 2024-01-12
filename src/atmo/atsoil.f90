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
!> \defgroup at_soil
!=============================================================================
!> \addtogroup at_soil
!> \{
implicit none
!> Number of boundary faces with soil features
integer, save :: nfmodsol
!> Number of soil model's type
integer, save :: nbrsol

!> kind of soil (water, forest, urban ...) and associated constantes
type categorie_sol
!> Dynamic roughness length
  double precision  :: rugdyn ! => field
!> Thermal roughness length
  double precision  :: rugthe! => field
!> Albedo
  double precision  :: albedo! => field
!> emissivity
  double precision  :: emissi ! => field
!> Vegetation index
  double precision  :: vegeta !> field
!> maximum water capacity of shallow reservoir
  double precision  :: c1w
!> ratio of the maximum water capacity of the shallow reservoir to the deep
!> reservoir [0,1]
  double precision  :: c2w
!> Thermal inertia of the soil
  double precision  :: csol
!> Rij value for Rij1
  double precision  :: r1
!> Rij value for Rij2
  double precision  :: r2
!> deep soil temperture
  double precision  :: tprof
!> Soil category name
  character(len=10) :: nomcat
end type categorie_sol

! Initialisation values for soil variables
! (filled in usispu/cs_user_parameters)

!> initial soil surface temperature
!> for Sea, it is also the surface temperature
double precision :: tsini = 20.d0
!> initial deep soil temperature
double precision :: tprini = 20.d0
!> initial soil specific humidity
double precision :: qvsini = 0.d0
!> initial water content of the first reservoir
double precision :: w1ini = 0.d0
!> initial water content of the second reservoir
double precision :: w2ini = 0.d0

!> array of the different features of each soil category
type(categorie_sol), dimension(:) , allocatable :: tab_sol

!> Defines the soil constants and variables of the vertical arrays
!> used for the 1D radiative model
type soil_tab
!> albedo
  double precision :: albedo
  !> emissivity
  double precision :: emissi
  !> soil thermo temperature
  double precision :: ttsoil
  !> soil potential temperature
  double precision :: tpsoil
  !> total water content
  double precision :: totwat
  !> surface pressure
  double precision :: pressure
  !> density
  double precision :: density
end type soil_tab

!> Class soilvert dimension
type(soil_tab), dimension(:), allocatable, save :: soilvert
type(soil_tab), save :: soil_mean

!> \}
!=============================================================================

end module atsoil

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
  double precision  :: rugdyn
!> Thermal roughness length
  double precision  :: rugthe
!> Albedo
  double precision  :: albedo
!> emissivity
  double precision  :: emissi
!> Vegetation index
  double precision  :: vegeta
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

!> Class definition of soil_variables
type variables_sol
  type(categorie_sol) :: constantes
!> soil temperature
  double precision :: temp_sol
!> deep soil temperature
  double precision :: tempp
!> total water content
  double precision :: total_water
!> ratio of the shallow reservoir water content to its maximum capacity
  double precision :: w1
!> ratio of the deep reservoir water content to its maximum capacity
  double precision :: w2
end type variables_sol

! Initialisation values for soil variables (filled in usati1)

!> initial soil surface temperature
double precision :: tsini
!> initial deep soil temperature
double precision :: tprini
!> initial soil specific humidity
double precision :: qvsini
!> Sea surface temperature
double precision :: tmer

!> array of the different features of each soil category
type(categorie_sol) , dimension(:) , allocatable :: tab_sol
!> index of boundary faces with soil features
integer , dimension(:) , allocatable, save       :: indsol
!> percentage of soil's category in each boundary face
integer , dimension(:,:) , allocatable           :: pourcent_sol
!> Class soil variable dimension
type(variables_sol) , dimension(:) , allocatable, save :: solution_sol

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
  !> ir downcoming flux
  double precision :: foir
  !> solar radation absorbed by the soil
  double precision :: fos
  end type soil_tab

!> Class soilvert dimension
type(soil_tab), dimension(:), allocatable, save :: soilvert
!> \}
!=============================================================================

end module atsoil

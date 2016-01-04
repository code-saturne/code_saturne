!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!> Module for the atmospheric soil model structures.

module atsoil

!=============================================================================

implicit none

integer, save :: nfmodsol , nbrsol

type categorie_sol
  double precision  :: rugdyn
  double precision  :: rugthe
  double precision  :: albedo
  double precision  :: emissi
  double precision  :: vegeta
  double precision  :: c1w
  double precision  :: c2w
  double precision  :: csol
  double precision  :: r1
  double precision  :: r2
  double precision  :: tprof
  character(len=10) :: nomcat
end type categorie_sol

type variables_sol
  type(categorie_sol) :: constantes
  double precision :: temp_sol      ! itempl
  double precision :: tempp
  double precision :: total_water   ! itotwt
  double precision :: w1            ! sol water content
  double precision :: w2
end type variables_sol

! Initialisation values for soil variables (filled in usati1)
double precision :: tsini   ! TSINI  : initial soil surface temperature
double precision :: tprini  ! TPRINI : initial deep soil temperature
double precision :: qvsini  ! QVSINI : initial soil specific humidity
double precision :: tmer    ! sea temperature

type(categorie_sol) , dimension(:) , allocatable :: tab_sol
integer , dimension(:) , allocatable, save       :: indsol
integer , dimension(:,:) , allocatable           :: pourcent_sol
type(variables_sol) , dimension(:) , allocatable, save :: solution_sol

! Defines the soil constants and variables of the vertical arrays
! used for the 1D radiative model

type soil_tab
  double precision :: albedo   ! albedo
  double precision :: emissi   ! emissivity
  double precision :: ttsoil   ! soil thermo temperature
  double precision :: tpsoil   ! soil potential temperature
  double precision :: totwat   ! total water content
  double precision :: pressure ! surface pressure
  double precision :: density  ! density
  double precision :: foir     ! ir downcoming flux
  double precision :: fos      ! solar radation absorbed by the soil
end type soil_tab

type(soil_tab), dimension(:), allocatable, save :: soilvert

!=============================================================================

end module atsoil

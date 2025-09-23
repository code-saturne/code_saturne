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

!==============================================================================
! Function:
! ---------

!> \file atimbr.f90
!>
!> \brief Atmospheric Imbrication module.
!>  This module contains the data structure and subroutines to
!>  perform atmospheric imbrication or nesting of a CFD domain within
!>  a large scale meteorological field.
!>  Starting from a set of large  scale meteological profiles (in the
!>  format of meteo files) an interpolation is performed for each
!>  boundary face both spatially and temporally (using Cressman method)
!>
!>  This imbrication reads additional meteo files.
!______________________________________________________________________________

module atimbr
!> \defgroup at_imbrication

use entsor
use atincl
implicit none

!> \addtogroup at_imbrication
!> \{
! --------------------------------------------------------------
!> activation flag
! --------------------------------------------------------------
logical(c_bool), pointer :: imbrication_flag
logical(c_bool), pointer, save :: imbrication_verbose

! --------------------------------------------------------------
!> Flags for activating the cressman interpolation for the boundary
!> conditions
! --------------------------------------------------------------
logical(c_bool), pointer :: cressman_u
logical(c_bool), pointer :: cressman_v
logical(c_bool), pointer :: cressman_tke
logical(c_bool), pointer :: cressman_eps
logical(c_bool), pointer :: cressman_theta
logical(c_bool), pointer :: cressman_qw
logical(c_bool), pointer :: cressman_nc

! --------------------------------------------------------------
!> numerical parameters for the cressman interpolation formulas
! --------------------------------------------------------------
real(c_double), pointer :: vertical_influence_radius
real(c_double), pointer :: horizontal_influence_radius

! --------------------------------------------------------------
!> Parameter for "meteo" files
! --------------------------------------------------------------
integer line_len
parameter(line_len = 132)
character(line_len) :: imbrication_files_list
character(line_len), dimension(:), allocatable :: imbrication_files
integer number_of_files
character*(3) skip_chars
parameter (skip_chars = "/#!")
!> Profile dimension variable
integer thermal_profile_dim
integer dynamical_profile_dim
!> Time sections per files
integer sections_per_file
data dynamical_profile_dim /-1/
data thermal_profile_dim /-1/
data sections_per_file /-1/

! --------------------------------------------------------------
!> read data from "meteo" files
! --------------------------------------------------------------
!> Time variables
integer, dimension(:,:), allocatable :: years
integer, dimension(:,:), allocatable :: ordinals
integer, dimension(:,:), allocatable :: hours
integer, dimension(:,:), allocatable :: minutes
double precision, dimension(:,:),allocatable :: seconds
!> Positions
double precision, dimension(:,:), allocatable :: xpos
double precision, dimension(:,:), allocatable :: ypos
double precision, dimension(:,:), allocatable :: ground_pressure
!> Vertical grid for temperature and humidity variables
double precision, dimension(:,:,:), allocatable :: zt
double precision, dimension(:,:,:), allocatable :: tempC
double precision, dimension(:,:,:), allocatable :: qw
double precision, dimension(:,:,:), allocatable :: Nc
!> Vertical grid for wind variables
double precision, dimension(:,:,:), allocatable :: zd
double precision, dimension(:,:,:), allocatable :: u
double precision, dimension(:,:,:), allocatable :: v
double precision, dimension(:,:,:), allocatable :: tke
double precision, dimension(:,:,:), allocatable :: eps

! --------------------------------------------------------------
!> derived data
! --------------------------------------------------------------
double precision, dimension(:,:), allocatable,target :: times
double precision, dimension(:,:,:), allocatable :: pressure
double precision, dimension(:,:,:), allocatable :: theta
double precision, dimension(:,:,:), allocatable :: density

! --------------------------------------------------------------
!> time interpolated profiles
! --------------------------------------------------------------
double precision, dimension(:,:), allocatable :: ti_zt
double precision, dimension(:,:), allocatable :: ti_tempC
double precision, dimension(:,:), allocatable :: ti_qw
double precision, dimension(:,:), allocatable :: ti_Nc
double precision, dimension(:,:), allocatable :: ti_zd
double precision, dimension(:,:), allocatable :: ti_u
double precision, dimension(:,:), allocatable :: ti_v
double precision, dimension(:,:), allocatable :: ti_tke
double precision, dimension(:,:), allocatable :: ti_eps
double precision, dimension(:,:), allocatable :: ti_pressure
double precision, dimension(:,:), allocatable :: ti_theta
double precision, dimension(:,:), allocatable :: ti_density

! --------------------------------------------------------------
!> additional variables
! --------------------------------------------------------------
double precision, dimension(:,:,:), allocatable :: coordinates_th
double precision, dimension(:,:,:), allocatable :: influence_param_th
double precision, dimension(:,:,:), allocatable :: coordinates_dyn
double precision, dimension(:,:,:), allocatable :: influence_param_dyn
integer(c_int), pointer ::  id_u
integer(c_int), pointer ::  id_v
integer(c_int), pointer ::  id_qw
integer(c_int), pointer ::  id_nc
integer(c_int), pointer ::  id_tke
integer(c_int), pointer ::  id_eps
integer(c_int), pointer ::  id_theta

interface

  subroutine cs_f_atmo_get_pointers_imbrication(p_imbrication_flag,       &
    p_imbrication_verbose, p_cressman_u, p_cressman_v, p_cressman_qw,     &
    p_cressman_nc, p_cressman_tke, p_cressman_eps, p_cressman_theta,      &
    p_vertical_influence_radius, p_horizontal_influence_radius,           &
    p_id_u, p_id_v, p_id_nc, p_id_qw, p_id_tke, p_id_eps, p_id_theta)     &
    bind(C, name='cs_f_atmo_get_pointers_imbrication')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_imbrication_flag, p_cressman_qw
      type(c_ptr), intent(out) :: p_imbrication_verbose, p_cressman_u
      type(c_ptr), intent(out) :: p_horizontal_influence_radius, p_cressman_v
      type(c_ptr), intent(out) :: p_cressman_nc, p_cressman_tke, p_cressman_eps
      type(c_ptr), intent(out) :: p_vertical_influence_radius, p_cressman_theta
      type(c_ptr), intent(out) :: p_id_u, p_id_v, p_id_nc, p_id_qw, p_id_tke
      type(c_ptr), intent(out) :: p_id_eps, p_id_theta
    end subroutine cs_f_atmo_get_pointers_imbrication

end interface
!> \}
contains

! ----------------------------------------------------------------
!> \brief  Map Fortran to C variables
! ----------------------------------------------------------------

subroutine atmo_init_imbrication()
  implicit none
  type(c_ptr) :: c_imbrication_flag, c_cressman_qw
  type(c_ptr) :: c_imbrication_verbose, c_cressman_u
  type(c_ptr) :: c_horizontal_influence_radius, c_cressman_v
  type(c_ptr) :: c_cressman_nc, c_cressman_tke, c_cressman_eps
  type(c_ptr) :: c_vertical_influence_radius, c_cressman_theta
  type(c_ptr) :: c_id_u, c_id_v, c_id_qw, c_id_nc, c_id_tke, c_id_eps, c_id_theta

  call cs_f_atmo_get_pointers_imbrication(c_imbrication_flag,                &
       c_imbrication_verbose, c_cressman_u, c_cressman_v, c_cressman_qw,     &
       c_cressman_nc, c_cressman_tke, c_cressman_eps, c_cressman_theta,      &
       c_vertical_influence_radius, c_horizontal_influence_radius,           &
       c_id_u, c_id_v, c_id_qw,c_id_nc, c_id_tke, c_id_eps, c_id_theta)

  call c_f_pointer(c_imbrication_flag, imbrication_flag)
  call c_f_pointer(c_imbrication_verbose, imbrication_verbose)

  call c_f_pointer(c_cressman_u, cressman_u)
  call c_f_pointer(c_cressman_v, cressman_v)
  call c_f_pointer(c_cressman_qw, cressman_qw)
  call c_f_pointer(c_cressman_nc, cressman_nc)
  call c_f_pointer(c_cressman_tke, cressman_tke)
  call c_f_pointer(c_cressman_eps, cressman_eps)
  call c_f_pointer(c_cressman_theta, cressman_theta)

  call c_f_pointer(c_vertical_influence_radius, vertical_influence_radius)
  call c_f_pointer(c_horizontal_influence_radius, horizontal_influence_radius)

  call c_f_pointer(c_id_u, id_u)
  call c_f_pointer(c_id_v, id_v)
  call c_f_pointer(c_id_qw, id_qw)
  call c_f_pointer(c_id_nc, id_nc)
  call c_f_pointer(c_id_tke, id_tke)
  call c_f_pointer(c_id_eps, id_eps)
  call c_f_pointer(c_id_theta, id_theta)

end subroutine atmo_init_imbrication

end module atimbr

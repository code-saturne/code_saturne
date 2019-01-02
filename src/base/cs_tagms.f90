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

!> \file cs_tagms.f90
!> Module for parameters options and physical properties of the 0-D thermal
!> model used by the metal mass structures modelling coupled with specific
!> condensation correlations.

module cs_tagms

  !=============================================================================

  implicit none

  !=============================================================================
  !> \defgroup cs_tagms Module for calculation options

  !> \addtogroup cs_tagms
  !> \{

  !----------------------------------------------------------------------------
  ! Numerical parameters to solve the 0-D thermal model of the metal mass
  !----------------------------------------------------------------------------

  !> \defgroup cs_tagms_numerical_parameters Numerical parameters of 0-D thermal model

  !> \addtogroup cs_tagms_numerical_parameters
  !> \{

  !> \anchor  xem
  !> the wall thickness of the metal mass structures used by the 0-D thermal model
  double precision xem

  !> \}

  !----------------------------------------------------------------------------
  ! Physical parameters to solve the 0-D thermal model of the metal mass
  !----------------------------------------------------------------------------

  !> \defgroup cs_tagms_physical_properties Physical properties of 0-D thermal model

  !> \addtogroup cs_tagms_physcial_properties
  !> \{

  !> \anchor xro
  !> the density (kg.m-3) of the metal mass structures
  double precision xro_m

  !> \anchor xcond
  !> the conductivity coefficient (W.m-1.C-1) of the metal mass structures
  double precision xcond_m

  !> \anchor xcp
  !> the specific heat coefficient (J.kg-1.C-1) of the metal mass structures
  double precision xcp_m

  !> \anchor m_metal
  !> the metal mass (kg) of the metal mass structures
  double precision m_metal

  !> \anchor s_metal
  !> the exchange surface (m2) of the metal mass structures
  double precision s_metal

  !> \anchor tmet0
  !> the initial temperature of the metal mass structures
  double precision tmet0

  !> \anchor t_metal
  !> the wall temperature computed with the 0-D thermal model
  !> associated to the metal mass structure material
  double precision, dimension(:,:),  allocatable :: t_metal

  !> \}

  !> \}

contains

  !=============================================================================

  subroutine init_tagms

  use mesh, only: ncelet

  allocate(t_metal(ncelet,2))

  !---> Array initialization
  t_metal(:,:) = 0.d0

  end subroutine init_tagms

  !=============================================================================

  subroutine finalize_tagms

    deallocate(t_metal)

  end subroutine finalize_tagms

  !=============================================================================


end module cs_tagms

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_tagmr.f90
!> Module for parameters options and physical properties of the 1-D thermal
!> model coupled with condensation.

module cs_tagmr

  !=============================================================================

  implicit none

  !=============================================================================
  !> \defgroup cs_tagmr Module for calculation options

  !> \addtogroup cs_tagmr
  !> \{

  !----------------------------------------------------------------------------
  ! Numerical parameters to solve the 1-D thermal conduction Equation
  !----------------------------------------------------------------------------

  !> \defgroup numerical_parameters of 1-D thermal conduction Equation

  !> \addtogroup numerical_parameters
  !> \{

  !> \anchor nmur
  !> number of discretized points associated to the (ii)th face
  !> with the 1-D thermal model coupled with condensation
  integer, save :: nmur

  !> \anchor dxp
  !> space step associated to the spatial discretization of the 1-D thermal model
  !> coupled with condensation model
  double precision, dimension(:),  allocatable ::  dxp

  !> \anchor  theta-scheme of the 1-D thermal model
  !>    - 0 : explicit scheme
  !>    - 1 : implicit scheme
  double precision theta

  !> \anchor  dxmin
  !> the minimal space step of 1-D thermal model
  !> by default equal to 0 with a homogeneus space step.
  !> this numerical parameter is used to impose a geometric progression
  !> ratio of the mesh refinement associated to (ii)th face with
  !> the 1-D thermal model.
  double precision dxmin

  !> \anchor  epais
  !> the wall thickness associated to the (ii)th face with 1-D thermal module
  double precision epais

  !> \}

  !----------------------------------------------------------------------------
  ! Physical parameters to solve the 1-D thermal conduction equation
  !----------------------------------------------------------------------------

  !> \defgroup physical_properties of 1-D thermal conduction Equation

  !> \addtogroup physcial_properties
  !> \{

  !> \anchor xrob
  !> the concrete density associated to solid material
  double precision rob

  !> \anchor condb
  !> the concrete conductivity coefficient associated to solid material
  double precision condb

  !> \anchor cpb
  !> the concrete specific heat coefficient associated to solid material
  double precision cpb

  !> \anchor hext
  !> the exterior exchange coefficient associated to solid material
  double precision hext

  !> \anchor text
  !> the exterior temperature associated to solid material
  double precision text

  !> \anchor tpar0
  !> the initial temperature associated to solid material
  double precision tpar0

  !> \anchor tmur
  !> the wall temperature computed with the 1-D thermal model
  !> associated to concrete solid material
  double precision, dimension(:,:),  allocatable :: tmur

  !> \}

  !> \}

contains

  !=============================================================================

  subroutine init_tagmr

  use pointe, only:nfbpcd

  allocate(dxp(nmur))
  allocate(tmur(nfbpcd,nmur))

  !---> Array initialization
  dxp(:)    = 0.d0
  tmur(:,:) = 0.d0

  end subroutine init_tagmr

  !=============================================================================

  subroutine finalize_tagmr

    deallocate(dxp)
    deallocate(tmur)

  end subroutine finalize_tagmr

  !=============================================================================

end module cs_tagmr

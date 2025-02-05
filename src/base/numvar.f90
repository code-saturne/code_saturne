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

!> \file numvar.f90
!> \brief Module for variable numbering

module numvar

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup numvar Module for variable numbering

  !> \addtogroup numvar
  !> \{

  !----------------------------------------------------------------------------
  ! Main variables
  !----------------------------------------------------------------------------

  !> \defgroup main_variables Main variables
  !> \brief Main variable field indices (previously stored in rtp, rtpa).

  !> \addtogroup main_variables
  !> \{

  !> \anchor ipr
  !> pressure
  integer, save :: ipr = 0

  !> \anchor iu
  !> velocity component \f$ u_x \f$
  integer, save :: iu = 0

  !> \anchor iv
  !> velocity component \f$ u_y \f$
  integer, save :: iv = 0

  !> \anchor iw
  !> velocity component \f$ u_z \f$
  integer, save :: iw = 0

  !> \anchor ik
  !> turbulent kinetic energy \f$ k \f$
  integer, save :: ik = 0

  !> \anchor iep
  !> turbulent dissipation \f$ \varepsilon \f$
  integer, save :: iep = 0

  !> \anchor irij
  !> Reynolds stress tenso \f$ R_{ij} \f$
  integer, save :: irij = 0

  !> \anchor iomg
  !> variable \f$ \omega \f$ of the \f$ k-\omega \f$ SST
  integer, save :: iomg = 0

  !> \anchor isca
  !> isca(i) is the index of the scalar i
  integer, save :: isca(nscamx)

  !> \anchor iscapp
  !> iscapp(i) is the index of the specific physics scalar i
  integer, save :: iscapp(nscamx)

  !> \anchor nscaus
  !> number of user scalars solutions of an advection equation
  integer, save :: nscaus = 0

  !> \anchor nscapp
  !> number of specific physics scalars
  integer, save :: nscapp = 0

  !> \}

  !----------------------------------------------------------------------------
  ! Physical properties
  !----------------------------------------------------------------------------

  !> \defgroup physical_prop Physical properties
  !> \brief Physical properties are using the field API.
  !> See \ref cs_user_boundary_conditions for some examples.

  !> \addtogroup physical_prop
  !> \{

  !> dynamic molecular viscosity (in kg/(m.s))
  integer, save :: iviscl = -1

  !> dynamic turbulent viscosity
  integer, save :: ivisct = -1

  !> interior and boundary convective mass flux key ids of the variables
  integer, save :: kimasf = -1, kbmasf = -1

  !> constant diffusivity field id key for scalars
  integer, save :: kvisl0 = -1

  !> variable diffusivity field id key for scalars
  integer, save :: kivisl = -1

  !> turbulent schmidt key for scalars
  integer, save :: ksigmas = -1

  !> cell density field ids of the variables
  integer, save :: icrom = -1

  !> boundary density field ids of the variables
  integer, save :: ibrom = -1

  !> \}

  !----------------------------------------------------------------------------
  ! Mapping to field structures
  !----------------------------------------------------------------------------

  !> \defgroup field_map Mapping to field structures

  !> \addtogroup field_map
  !> \{

  !> Field id for variable i
  integer, save :: ivarfl(nvarmx)

  !> \}

  !=============================================================================

  !> \}

end module numvar

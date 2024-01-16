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

!> \file cs_nz_condensation.f90
!> Module for parameters options and physical properties of the condensation
!> model using specific zones with different wall temperatures and material
!> properties.

module cs_nz_condensation

  !=============================================================================

  use, intrinsic :: iso_c_binding
  implicit none

  !=============================================================================
  !> \defgroup cs_nz_condensation Module for calculation options

  !> \addtogroup cs_nz_condensation
  !> \{

  !----------------------------------------------------------------------------
  ! Parameters to compute the condensation source terms
  !----------------------------------------------------------------------------

  !> \defgroup Parameters to choice the exchange coefficients type and the
  !> wall temperature used by the condensation model
  !> (with or without 1D thermal model to the wall).

  !> \addtogroup condensation and thermal exchange coefficients
  !> \{

  !> \anchor nfbpcd
  !> number of faces in which a condensation source terms is imposed.
  !> See \c ifbpcd and the user subroutine \ref cs_user_wall_condensation
  integer(c_int), pointer, save :: nfbpcd

  !> \anchor itypcd
  !> type of condensation source terms for each variable
  !> - 0 for an variable at ambient value,
  !> - 1 for an variable at imposed value.
  !> See the user subroutine \ref cs_user_wall_condensation
  integer, dimension(:,:), pointer, save :: itypcd

  !> \anchor thermal_condensation_flux
  !> value of the thermal flux for the condensation model.
  !> See the user subroutine \ref cs_user_wall_condensation
  double precision, dimension(:), pointer, save :: thermal_condensation_flux

  !> \anchor flthr
  !> external heat flux used as flux conditions
  !> of the 1d thermal model (in unit \f$W.m^{-2}\f$).
  double precision, dimension(:), pointer, save :: flthr

  !> \anchor dflthr
  !> external heat flux derivative used as flux conditions
  !> of the 1d thermal model (in unit \f$W.m^{-3}\f$).
  double precision, dimension(:), pointer, save :: dflthr

  !> \anchor spcond
  !> value of the condensation source terms for pressure.
  !> For the other variables, eventual imposed specific value.
  !> See the user subroutine \ref cs_user_wall_condensation
  double precision, dimension(:,:), pointer, save :: spcond

  !> \anchor twall_cond
  !> Temperature at condensing wall faces (for post-processing purposes)
  double precision, dimension(:), pointer, save :: twall_cond

  !> value of the thermal exchange coefficient associated to
  !> the condensation model used.
  !> See the user subroutine \ref cs_user_wall_condensation
  double precision, dimension(:), pointer, save:: hpcond

  !> list on the nfbpcd faces in which a condensation source terms is imposed.
  !> See \c ifbpcd and the user subroutine \ref cs_user_wall_condensation
  integer, dimension(:), pointer, save :: ifbpcd

  !> \anchor nzones
  !> number of the zones with a specific condensation source terms
  !> depending on the wall temperature and material properties.
  !> by default (nzones = 1) if the user does not specified different zones
  !> in the user subroutine \ref cs_user_wall_condensation.f90 .
  integer(c_int), pointer, save :: nzones

  !> \anchor izzftcd
  !> list of the zones associated to the faces where a condensation source term
  !> is imposed. On each zone a specific wall temperature and material properties
  !> can be specified by the user.
  integer, dimension(:), pointer, save :: izzftcd

  !> \anchor izcophc
  !> choice the way to compute the exchange coefficient of the
  !> condensation source term used by the copain model.
  !>    - 1: the turbulent exchange coefficient of the flow
  !>    - 2: the exchange coefficient of the copain correlation
  !>    - 3: the maximal value between the two previous exchange coefficients
  integer, dimension(:), pointer, save :: izcophc

  !> \anchor izcophg
  !> choice the way to compute the thermal exchange coefficient associated
  !> to the heat transfer to the wall due to the condensation phenomenon.
  !>    - 2: the thermal exchange coefficient of the copain correlation
  !>    - 3: the maximal value between the current and previous thermal
  !>         exchange coefficient evaluated by the copain correlation
  integer, dimension(:), pointer, save :: izcophg

  !> \}

  !> \addtogroup Wall temperature model
  !> \{

  !> \anchor iztag1d
  !> choice the way to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the wall
  !>    - 1: the wall temperature is computed with a 1-D thermal model
  !>         with implicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         exchange coefficient evaluated by the copain correlation
  integer, dimension(:), pointer, save :: iztag1d

  !> \anchor nztag1d
  !> Indicate if the thermal 1D model of severe accident is usde
  !> to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the wall
  !>    - 1: the wall temperature is computed with a 1-D thermal model
  !>         with implicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         exchange coefficient evaluated by the copain correlation
  integer(c_int), pointer, save :: nztag1d

  !> Constant value of the wall temperature given by the user when
  !> the thermal 1D model is not activated for the condensation model with
  !> different zones specified in the user subroutine
  !> \ref cs_user_wall_condensation.
  double precision, dimension(:), pointer, save :: ztpar

  !> \anchor zxref
  !> Coordinates of the reference point for forced and mixed convection regimes
  !> index 1 : coordinate, index 2: zone_id
  double precision, dimension(:, :), pointer, save :: zxrefcond
  double precision, dimension(:, :), pointer, save :: zprojcond


  !> \anchor ncmast
  !> number of the cells in which a condensation source terms is imposed.
  !> See \c lstmast list and the subroutine \ref cs_user_wall_condensation
  integer(c_int), pointer, save :: ncmast

  !> \anchor nvolumes
  !> number of the volume strutures with a specific condensation source terms
  !> depending on the wall temperature and material properties.
  !> by default (nvolumes = 1) if the user does not specified different volumes
  !> in the user subroutine \ref cs_user_wall_condensation.f90 .
  integer(c_int), pointer, save :: nvolumes

  !> \anchor ltmast
  !> list on the ncmast cells in which a condensation source terms is imposed.
  !> See  the user subroutine \ref cs_user_wall_condensation.
  integer, dimension(:), pointer, save :: ltmast

  !> zone type where a condensation source terms is imposed to model
  !> the metal structures condensation on a volumic zone.
  integer, dimension(:), pointer, save :: izmast

  !> type of condensation source terms for each variable
  !> - 0 for a variable at ambient value,
  !> - 1 for a variable at imposed value.
  !> See the user function \ref  cs_user_wall_condensation.
  integer, dimension(:,:), pointer, save :: itypst

  !> value of the condensation source terms for pressure
  !> associated to the metal structures modelling.
  !> For the other variables, eventual imposed specific value.
  !> See the user function \ref cs_user_wall_condensation.
  double precision, dimension(:, :), pointer, save :: svcond

  !> value of the thermal flux for the condensation model
  !> associated to the metal structures modelling.
  !> See the user function \ref cs_user_wall_condensation.
  double precision, dimension(:), pointer, save :: flxmst

  !> choice the way to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the metal mass structures wall
  !>    - 1: the wall temperature is computed with a 0-D thermal model
  !>         with explicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         and past to the copain correlation to evaluate the exchange
  !>         coefficient
  integer, dimension(:), pointer, save :: itagms

  !> \}
  !> \}

interface

  !---------------------------------------------------------------------------
  !> \brief Create wall condensation structure
  !
  !> \param[in]   nfbpcd    Number of faces with wall condensation activated
  !> \param[in]   nvar      Number of variables (?)
  !---------------------------------------------------------------------------

  subroutine cs_f_wall_condensation_create(nfbpcd, nzones, ncmast,  &
                                           nvolumes, nvar) &
    bind(C, name='cs_wall_condensation_create')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: nfbpcd, nzones, ncmast, nvolumes,  nvar
  end subroutine cs_f_wall_condensation_create

  !---------------------------------------------------------------------------
  !> \brief Return pointers to spcond
  !
  !> \param[out]   spcond   Pointer to spcond
  !---------------------------------------------------------------------------

  subroutine cs_f_wall_condensation_get_size_pointers(nfbpcd, nzones, &
                                                      ncmast, nvolumes, &
                                                      nztag1d) &
    bind(C, name='cs_f_wall_condensation_get_size_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: nfbpcd, nzones, ncmast, nvolumes, nztag1d
  end subroutine cs_f_wall_condensation_get_size_pointers

  !---------------------------------------------------------------------------
  !> \brief Return pointers to spcond
  !
  !> \param[out]   spcond   Pointer to spcond
  !---------------------------------------------------------------------------

  subroutine cs_f_wall_condensation_get_pointers(ifbpcd, itypcd, izzftcd, &
                                                 spcond, hpcond, twall_cond, &
                                                 thermflux, flthr, dflthr, &
                                                 izcophc, izcophg, iztag1d, &
                                                 ztpar, zxrefcond, zprojcond, &
                                                 ltmast, itypst, izmast, &
                                                 svcond, flxmst, itagms) &
    bind(C, name='cs_f_wall_condensation_get_pointers')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: ifbpcd, itypcd, izzftcd
    type(c_ptr), intent(out) :: spcond, hpcond, twall_cond
    type(c_ptr), intent(out) :: thermflux, flthr, dflthr, izcophc
    type(c_ptr), intent(out) :: izcophg, iztag1d, ztpar, zxrefcond, zprojcond
    type(c_ptr), intent(out) :: ltmast, itypst, izmast, svcond, flxmst, itagms
  end subroutine cs_f_wall_condensation_get_pointers

end interface

contains

  subroutine init_sizes_pcond()
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr) :: c_nfbpcd, c_nzones, c_ncmast, c_nvolumes, c_nztag1d

    call cs_f_wall_condensation_get_size_pointers(c_nfbpcd, c_nzones, &
                                                  c_ncmast, c_nvolumes, &
                                                  c_nztag1d)
    call c_f_pointer(c_nfbpcd, nfbpcd)
    call c_f_pointer(c_nzones, nzones)
    call c_f_pointer(c_ncmast, ncmast)
    call c_f_pointer(c_nvolumes, nvolumes)
    call c_f_pointer(c_nztag1d, nztag1d)

  end subroutine init_sizes_pcond

  !=============================================================================

  subroutine init_nz_pcond(nvar)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in) :: nvar

    ! Local variables

    type(c_ptr) :: c_spcond, c_hpcond, c_ifbpcd, c_walltemp
    type(c_ptr) :: c_izzftcd, c_itypcd, c_thermflux, c_flthr, c_dflthr
    type(c_ptr) :: c_izcophc, c_izcophg, c_iztag1d, c_ztpar
    type(c_ptr) :: c_zxrefcond, c_zprojcond
    type(c_ptr) :: c_ltmast, c_itypst, c_izmast, c_svcond, c_flxmst, c_itagms

    if (nzones < 1) nzones = 1
    if (nvolumes < 1) nvolumes = 1

    call cs_f_wall_condensation_create(nfbpcd, nzones, ncmast, nvolumes, nvar)
    call cs_f_wall_condensation_get_pointers(c_ifbpcd, c_itypcd, c_izzftcd,  &
                                             c_spcond, c_hpcond, c_walltemp, &
                                             c_thermflux, c_flthr, c_dflthr, &
                                             c_izcophc, c_izcophg, c_iztag1d, &
                                             c_ztpar, c_zxrefcond, c_zprojcond, &
                                             c_ltmast, c_itypst, c_izmast, &
                                             c_svcond, c_flxmst, c_itagms)
    call c_f_pointer(c_ifbpcd, ifbpcd, [nfbpcd])
    call c_f_pointer(c_itypcd, itypcd, [nfbpcd, nvar])
    call c_f_pointer(c_izzftcd, izzftcd, [nfbpcd])
    call c_f_pointer(c_spcond, spcond, [nfbpcd, nvar])
    call c_f_pointer(c_hpcond, hpcond, [nfbpcd])
    call c_f_pointer(c_walltemp, twall_cond, [nfbpcd])
    call c_f_pointer(c_thermflux, thermal_condensation_flux, [nfbpcd])
    call c_f_pointer(c_flthr, flthr, [nfbpcd])
    call c_f_pointer(c_dflthr, dflthr, [nfbpcd])
    call c_f_pointer(c_izcophc, izcophc, [nzones])
    call c_f_pointer(c_izcophg, izcophg, [nzones])
    call c_f_pointer(c_iztag1d, iztag1d, [nzones])
    call c_f_pointer(c_ztpar, ztpar, [nzones])
    call c_f_pointer(c_zxrefcond, zxrefcond, [3, nzones])
    call c_f_pointer(c_zprojcond, zprojcond, [3, nzones])
    call c_f_pointer(c_ltmast, ltmast, [ncmast])
    call c_f_pointer(c_itypst, itypst, [ncmast, nvar])
    call c_f_pointer(c_izmast, izmast, [ncmast])
    call c_f_pointer(c_svcond, svcond, [ncmast, nvar])
    call c_f_pointer(c_flxmst, flxmst, [ncmast])
    call c_f_pointer(c_itagms, itagms, [nvolumes])

  end subroutine init_nz_pcond

  !=============================================================================

end module cs_nz_condensation

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file lagran.f90
!> \brief Module for Lagrangian model.

module lagran

  !===========================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !> \defgroup lagran Module for Lagrangian model

  !> \addtogroup lagran
  !> \{

  !=============================================================================

  !> \defgroup base Base

  !> \addtogroup base
  !> \{

  !> \anchor iilagr
  !> activates (>0) or deactivates (=0) the Lagrangian module
  !> the different values correspond to the following modellings:
  !> - = 1 Lagrangian two-phase flow in one-way coupling (no influence of
  !> the particles on the continuous phase)
  !> - = 2 Lagrangian two-phase flow with two-way coupling (influence of
  !> the particles on the dynamics of the continuous phase).
  !> Dynamics, temperature and mass may be coupled independently.
  !> - = 3 Lagrangian two-phase flow on frozen continuous phase. This option can
  !> only be used in case of a calculation restart. All the
  !> Eulerian fields are frozen (including the scalar fields). This option
  !> automatically implies \ref iccvfg = 1
  integer(c_int), pointer, save :: iilagr

  !> \}

  !=============================================================================

  !> \defgroup lag_st_pointers Lagrangian source term pointers

  !> \addtogroup lag_st_pointers
  !> \{

  !> \anchor ntersl
  integer(c_int), pointer, save :: ntersl

  !> \anchor ptsvar
  double precision, dimension(:,:), pointer, save :: ptsvar

  !> \}

  !> \}

  !=============================================================================

  !> \defgroup specific_physic Specific physics

  !> \addtogroup specific_physic
  !> \{

  !> - 0: no deposition submodel activated,
  !> - 1: deposition submodel used
  integer(c_int), pointer, save :: idepst

  !> - 0: no head losses calculation for influence of the deposit on the flow
  !> - 1: head losses calculation for influence of the deposit on the flow
  integer(c_int), pointer, save :: iflow

  !> - 0: no precipitation/dissolution model
  !> - 1: precipitation/dissolution model
  integer(c_int), pointer, save :: ipreci

  !> \}

  !=============================================================================

  !> \defgroup source_terms Source terms

  !> \addtogroup source_terms
  !> \{

  !> activation (=1) or not (=0) of the two-way coupling on the dynamics
  !> of the continuous phase.
  !> Useful if \ref iilagr = 2 and \ref iccvfg = 0
  integer(c_int), pointer, save ::  ltsdyn

  !> activation (=1) or not (=0) of the two-way coupling on the mass.
  !> Useful if \ref iilagr = 2, \ref physical_model = 1 and \ref impvar = 1
  integer(c_int), pointer, save ::  ltsmas

  !> if \ref physical_model = 1 and \ref itpvar = 1, \ref ltsthe
  !> activates (=1) or not (=0) the two-way coupling on temperature.
  !> if \ref physical_model = 2, \ref ltsthe activates (=1) or not (=0) the
  !> two-way coupling on the eulerian variables related to pulverised
  !> coal combustion.
  !> Useful if \ref iilagr = 2
  integer(c_int), pointer, save ::  ltsthe

  !> implicit source term for the continuous phase velocity and
  !> for the turbulent energy if the \f$k-\varepsilon\f$ model is used
  integer(c_int), pointer, save ::  itsli

  !> explicit source term for the turbulent dissipation and the
  !> turbulent energy if the \f$k-\varepsilon\f$ turbulence model is used
  !> for the continuous phase
  integer(c_int), pointer, save ::  itske

  !> explicit thermal source term for the thermal scalar of the continuous phase
  integer(c_int), pointer, save ::  itste

  !> implicit thermal source term for the thermal scalar of the continuous phase
  integer(c_int), pointer, save ::  itsti

  !> mass source term
  integer(c_int), pointer, save ::  itsmas

  !> source term for the light volatile matters
  integer(c_int), dimension(:), pointer, save ::  itsmv1

  !> source term for the heavy volatile matters
  integer(c_int),  dimension(:), pointer, save ::  itsmv2

  !> source term for the carbon released during heterogeneous combustion
  integer(c_int), pointer, save ::  itsco

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function returning particle attribute pointers
    subroutine cs_f_lagr_dim_pointers(p_ntersl) &
      bind(C, name='cs_f_lagr_dim_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out)     :: p_ntersl
    end subroutine cs_f_lagr_dim_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function returning particle attribute pointers

    subroutine cs_f_lagr_params_pointers(p_iilagr, p_idepst, p_iflow, p_ipreci) &
      bind(C, name='cs_f_lagr_params_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_iilagr, p_idepst, p_iflow, p_ipreci
    end subroutine cs_f_lagr_params_pointers

    subroutine cs_f_lagr_source_terms_pointers(p_ltsdyn, p_ltsmas,             &
                                               p_ltsthe, p_itsli,              &
                                               p_itske, p_itste, p_itsti,      &
                                               p_itsmas, p_itsco,              &
                                               p_itsmv1, p_itsmv2, dim_itsmv1, &
                                               dim_itsmv2) &
      bind(C, name='cs_f_lagr_source_terms_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: dim_itsmv1, dim_itsmv2
      type(c_ptr), intent(out) :: p_ltsdyn, p_ltsmas, p_ltsthe, p_itsli,       &
                                  p_itske, p_itste,                            &
                                  p_itsti, p_itsmas, p_itsmv1, p_itsmv2,       &
                                  p_itsco
    end subroutine cs_f_lagr_source_terms_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function passing specific physics options

    subroutine  cs_f_lagr_specific_physics(iirayo, ncharb, ncharm, diftl0)     &
      bind(C, name='cs_f_lagr_specific_physics')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: iirayo, ncharb, ncharm
      real(c_double) :: diftl0
    end subroutine cs_f_lagr_specific_physics

    !---------------------------------------------------------------------------

    ! Interface to C function passing coal combustion parameters

    subroutine  cs_f_lagr_coal_comb(ih2o, io2, ico, iatc, prefth, trefth,      &
                                    natom, wmolat, ngazem, wmole, iym1,        &
                                    ncharm, a1ch, h02ch, e1ch, a2ch, e2ch,     &
                                    y1ch, y2ch, cp2ch, ahetch, ehetch,         &
                                    rho0ch, xwatch, xashch, thcdch)            &
      bind(C, name='cs_f_lagr_coal_comb')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: ih2o, io2, ico, iatc, natom, ngazem, ncharm
      real(c_double) :: prefth, trefth
      integer(c_int), dimension(ngazem) :: iym1
      real(c_double), dimension(natom) :: wmole, wmolat
      real(c_double), dimension(ncharm) :: a1ch, h02ch, e1ch, a2ch, e2ch,      &
                                           y1ch, y2ch, cp2ch, ahetch, ehetch,  &
                                           rho0ch, xwatch, xashch, thcdch
    end subroutine cs_f_lagr_coal_comb

    !---------------------------------------------------------------------------

    !> Prepare for execution of the Lagrangian model.

    subroutine cs_lagr_solve_initialize(dt)  &
      bind(C, name='cs_lagr_solve_initialize')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: dt
    end subroutine cs_lagr_solve_initialize

    !---------------------------------------------------------------------------

    !> Execute one time step of the Lagrangian model.

    subroutine cs_lagr_solve_time_step(itypfb, dt)  &
      bind(C, name='cs_lagr_solve_time_step')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(*), intent(in) :: itypfb
      real(kind=c_double), dimension(*), intent(in) :: dt
    end subroutine cs_lagr_solve_time_step

    !---------------------------------------------------------------------------

    !> \brief Allocate bound_stat and return fortran compatible pointer

    subroutine cs_lagr_init_c_arrays(dim_tslagr, p_tslagr)         &
      bind(C, name='cs_lagr_init_c_arrays')
      use, intrinsic ::  iso_c_binding

      implicit none
      integer(c_int), dimension(2) :: dim_tslagr
      type(c_ptr), intent(out)     :: p_tslagr
    end subroutine cs_lagr_init_c_arrays

    !---------------------------------------------------------------------------

    subroutine cs_lagr_init_par ()&
      bind(C, name='cs_lagr_init_par')

      use, intrinsic :: iso_c_binding

    end subroutine cs_lagr_init_par

    !---------------------------------------------------------------------------

    !> \brief Write particle data to checkpoint.

    !> \param[in, out]  r  pointer to restart structure

    function lagr_restart_write_particle_data(r) result(n) &
      bind(C, name='cs_lagr_restart_write_particle_data')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      integer(c_int) :: n
    end function lagr_restart_write_particle_data

  !=============================================================================

  end interface

contains

  !=============================================================================

  subroutine lagran_pointers()

    implicit none

    ! lagr_option_base
    type(c_ptr) :: p_iilagr

    ! lagr_params
    type(c_ptr) :: p_idepst, p_iflow, p_ipreci

    ! lagr_option_source_terms
    type(c_ptr) :: p_ltsdyn, p_ltsmas, p_ltsthe, p_itsli, p_itske,     &
                   p_itste, p_itsti,                                   &
                   p_itsmas, p_itsmv1, p_itsmv2, p_itsco
    integer(c_int) :: dim_itsmv1, dim_itsmv2

    call cs_f_lagr_params_pointers(p_iilagr, p_idepst, p_iflow, p_ipreci)
    call c_f_pointer(p_iilagr, iilagr)
    call c_f_pointer(p_idepst, idepst)
    call c_f_pointer(p_iflow , iflow)
    call c_f_pointer(p_ipreci, ipreci)

    call cs_f_lagr_source_terms_pointers(p_ltsdyn, p_ltsmas,           &
                                         p_ltsthe, p_itsli,            &
                                         p_itske,  p_itste, p_itsti,   &
                                         p_itsmas, p_itsco,            &
                                         p_itsmv1, p_itsmv2,           &
                                         dim_itsmv1, dim_itsmv2)
    call c_f_pointer(p_ltsdyn, ltsdyn)
    call c_f_pointer(p_ltsmas, ltsmas)
    call c_f_pointer(p_ltsthe, ltsthe)
    call c_f_pointer(p_itsli , itsli )
    call c_f_pointer(p_itske , itske )
    call c_f_pointer(p_itste , itste )
    call c_f_pointer(p_itsti , itsti )
    call c_f_pointer(p_itsmas, itsmas)
    call c_f_pointer(p_itsco , itsco )
    call c_f_pointer(p_itsmv1, itsmv1, [dim_itsmv1])
    call c_f_pointer(p_itsmv2, itsmv2, [dim_itsmv2])

    return

  end subroutine lagran_pointers

  !=============================================================================

  subroutine init_lagr_dim_pointers()

    implicit none

    type(c_ptr) :: p_ntersl

    call cs_f_lagr_dim_pointers(p_ntersl)

    call c_f_pointer(p_ntersl, ntersl)

    return

  end subroutine init_lagr_dim_pointers

  !=============================================================================

  ! Initialize auxiliary arrays

  subroutine init_lagr_arrays(tslagr)

    implicit none

    double precision, dimension(:,:), pointer  :: tslagr
    integer(c_int),   dimension(2)             :: dim_tslagr
    type(c_ptr)                                :: p_tslagr

    call cs_lagr_init_c_arrays(dim_tslagr, p_tslagr)

    call c_f_pointer(p_tslagr, tslagr, [dim_tslagr])

    return

  end subroutine init_lagr_arrays

  !=============================================================================

  subroutine lagran_init_map

    use ppincl, only: iccoal, icfuel, ieljou, ielarc, icoebu, icod3p,          &
                      icpl3c, iym1, icompf
    use cpincl, only: ncharb, xashch, cp2ch, xwatch, rho0ch, a1ch,             &
                      a2ch, e1ch, e2ch, io2, ih2o, ico,                        &
                      ahetch, ehetch, thcdch, y1ch, y2ch, h02ch
    use ppppar, only: ncharm
    use ppthch, only: diftl0, ngazem, wmole, wmolat, trefth, prefth, iatc,     &
                      natom, wmolat
    use radiat, only: iirayo

    call init_lagr_dim_pointers

    call lagran_pointers

    call cs_f_lagr_specific_physics(iirayo,         &
                                    ncharb,         &
                                    ncharm,         &
                                    diftl0)

    call cs_f_lagr_coal_comb(ih2o,   &
                             io2,    &
                             ico,    &
                             iatc,   &
                             prefth, &
                             trefth, &
                             natom,  &
                             wmolat, &
                             ngazem, &
                             wmole,  &
                             iym1,   &
                             ncharm, &
                             a1ch,   &
                             h02ch,  &
                             e1ch,   &
                             a2ch,   &
                             e2ch,   &
                             y1ch,   &
                             y2ch,   &
                             cp2ch,  &
                             ahetch, &
                             ehetch, &
                             rho0ch, &
                             xwatch, &
                             xashch, &
                             thcdch)

  end subroutine lagran_init_map

  !=============================================================================

end module lagran

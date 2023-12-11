!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file ppthch.f90
!> Module for specific physics thermophysical data

module ppthch

  !===========================================================================

  use, intrinsic :: iso_c_binding

  use cstphy

  implicit none

  !===========================================================================

  !> \defgroup  thermophysical Module for specific physics thermophysical data

  !> \addtogroup thermophysical
  !> \{

  !> reference temperature for the specific physics, in K
  double precision trefth

  !> reference pressure for the specific physics, in Pa
  double precision prefth

  !> molar volume under normal pressure and temperature conditions
  !>  (1 atmosphere, 0 \f$\text{\degresC}\f$) in \f$m^{-3}\f$
  double precision volmol

  parameter ( trefth = 25.d0 + tkelvi ,                             &
              prefth = 1.01325d5      ,                             &
              volmol = 22.41d-3       )

  !--> DONNEES

  !> maximal number of global species
  integer    ngazgm

  !> maximal number of elementary gas components
  integer    ngazem

  !> maximal number of tabulation points
  integer    npot

  !> maximal number of atomic species
  integer    natom

  !> maximal number of global reactions in gas phase
  integer    nrgazm

  parameter( ngazgm = 25 , ngazem = 20 ,                                     &
             npot  = 500 , natom  = 5   , nrgazm = 1 )
  integer    iatc, iath, iato, iatn , iats
  parameter( iatc = 1, iath = 2, iato = 3, iatn = 4 , iats = 5 )

  !> name of global species
  character(len=150) :: nomcog(ngazgm)
  !> name of elementary species
  character(len=12) :: nomcoe(ngazem)

  !> number of tabulation points
  integer, pointer, save :: npo

  !> number of elementary gas components
  integer, pointer, save ::  ngaze
  !> number of global species
  integer, pointer, save ::  ngazg
  !> number of atomic species
  integer, pointer, save ::  nato

  !> number of global reactions in gas phase
  integer, pointer, save ::  nrgaz

  !> rank of O2 in gas composition
  integer, save ::           iio2
  !> rank of H2O in gas composition
  integer, save ::           iih2o
  !> rank of CO2 in gas composition
  integer, save ::           iico2
  !> rank of CO in gas composition
  integer, save ::           iico
  !> rank of C in gas composition
  integer, pointer, save ::  iic

  !> rank of fuel in the r-th reaction
  integer, save ::           igfuel(nrgazm)
  !> rank of oxydiser in the r-th reaction
  integer, save ::           igoxy(nrgazm)
  !> rank of products in the r-th reaction
  integer, save ::           igprod(nrgazm)

  !> stoechiometric coefficient of global species
  double precision, save ::  nreact(ngazgm)

  !> temperature (in K)
  real(c_double), pointer, save :: th(:)

  !> engaze(ij) is the massic enthalpy (J/kg) of the i-th elementary gas component
  !> at temperature  th(j)
  double precision, save ::  ehgaze(ngazem,npot)

  !> engazg(ij) is the massic enthalpy (J/kg) of the i-th global secies
  !> at temperature  th(j)
  double precision, save ::  ehgazg(ngazgm,npot)

  !> cpgazg(ij) is the massic calorific capacity (J/kg/K) of the i-th global secies
  !> at temperature  th(j)
  real(c_double), pointer, save ::  cpgazg(:,:)

  !> molar mass of an elementary gas component
  real(c_double), pointer, save ::  wmole(:)

  !> molar mass of a global species
  real(c_double), pointer, save ::  wmolg(:)

  !> molar mass of atoms
  double precision, save ::  wmolat(natom)

  !> Stoichiometry in reaction global species.  Negative for the reactants,
  !> and positive for the products
  double precision, save ::  stoeg(ngazgm,nrgazm)

  !> Mixing rate at the stoichiometry
  real(c_double), pointer, save ::  fs(:)

  !> Absorption coefficient of global species
  double precision, save ::  ckabsg(ngazgm)

  !> Absorption coefficient of gas mixture
  real(c_double), pointer, save ::  ckabs1

  !> \anchor diftl0
  !> molecular diffusivity for the enthalpy (\f$kg.m^{-1}.s^{-1}\f$)
  !> for gas or coal combustion (the code then automatically sets
  !> \ref diffusivity_ref to \ref diftl0 for the scalar
  !> representing the enthalpy).
  !>
  !> Always useful for gas or coal combustion.
  real(c_double), pointer, save ::  diftl0

  !> Molar coefficient of CO2
  real(c_double), pointer, save ::  xco2
  !> Molar coefficient of H2O
  real(c_double), pointer, save ::  xh2o

  !=============================================================================

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_ppthch_get_pointers(p_ngaze, p_ngazg, p_nato, p_nrgaz,         &
                                        p_iic, p_npo, p_wmole, p_wmolg, p_diftl0,  &
                                        p_xco2, p_xh2o, p_ckabs1, p_fs, p_th,      &
                                        p_cpgazg)      &
      bind(C, name='cs_f_ppthch_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ngaze, p_ngazg, p_nato, p_nrgaz, p_iic,    &
                                  p_wmolg, p_wmole, p_diftl0, p_xco2, p_xh2o,  &
                                  p_ckabs1,  p_fs, p_th, p_npo, p_cpgazg
    end subroutine cs_f_ppthch_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine thch_models_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ngaze, p_ngazg, p_nato, p_nrgaz, p_iic
    type(c_ptr) :: p_wmole, p_wmolg, p_xco2, p_xh2o, p_ckabs1
    type(c_ptr) :: p_fs, p_diftl0, p_th, p_npo, p_cpgazg

    call cs_f_ppthch_get_pointers(p_ngaze, p_ngazg, p_nato, p_nrgaz,     &
                                  p_iic,  p_npo,                         &
                                  p_wmole, p_wmolg, p_diftl0,            &
                                  p_xco2, p_xh2o, p_ckabs1, p_fs,        &
                                  p_th, p_cpgazg)

    call c_f_pointer(p_ngaze, ngaze)
    call c_f_pointer(p_ngazg, ngazg)
    call c_f_pointer(p_nato, nato)
    call c_f_pointer(p_nrgaz, nrgaz)
    call c_f_pointer(p_iic, iic)
    call c_f_pointer(p_wmole, wmole, [ngazem])
    call c_f_pointer(p_wmolg, wmolg, [ngazgm])
    call c_f_pointer(p_xco2, xco2)
    call c_f_pointer(p_diftl0, diftl0)
    call c_f_pointer(p_xh2o, xh2o)
    call c_f_pointer(p_ckabs1, ckabs1)
    call c_f_pointer(p_fs, fs, [nrgazm])
    call c_f_pointer(p_th, th, [npot])
    call c_f_pointer(p_npo, npo)
    call c_f_pointer(p_cpgazg, cpgazg, [ngazgm, npot])

  end subroutine thch_models_init

  !=============================================================================

end module ppthch

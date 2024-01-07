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

!> \file cpincl.f90
!> Module for pulverized coal combustion

module cpincl

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use ppppar
  use ppthch

  implicit none

  !=============================================================================

  ! Combustion variable pointers for pulverized coal

  ! see also ppppar: ncharm, ncpcmx, nclcpm

  ! epsicp: precision for tests

  double precision epsicp
  parameter ( epsicp = 1.d-8 )

  ! Data relative to coal

  !> Number of coals
  integer(c_int), pointer, save :: ncharb

  !> coal with drift (0: without drift (default), 1: with)
  integer(c_int), pointer, save  ::  i_comb_drift

  ! By coal (given quantities)

  ! Granulometric distribution

  !> Number of classes per coal
  integer(c_int), pointer, save :: nclpch(:)

  ! - Proprietes sur charbon sec

  real(c_double), pointer, save :: cch(:), hch (:), och(:), sch(:), nch(:),  &
                                  alpha(:), beta(:), teta(:), omega(:),      &
                                  pcich(:), rho0ch(:), thcdch(:),            &
                                  cck(:), hck(:), ock(:), sck(:), nck(:),    &
                                  gamma(:), delta(:), kappa(:), zeta(:),     &
                                  rhock(:), pcick(:),                        &
                                  cpashc(:), h0ashc(:), h02ch(:), cp2wat(:), &
                                  crepn1(:,:), crepn2(:,:),                  &
                                  cp2ch(:), xashsec(:), xashch(:), xwatch(:)

    integer(c_int), pointer, save :: iy1ch (:), iy2ch (:)
    integer(c_int), pointer, save :: iochet(:) , ioetc2(:), ioetwt(:)
    real(c_double), pointer, save :: y1ch(:), a1ch(:), e1ch(:),         &
                                     y2ch(:), a2ch(:), e2ch(:),         &
                                     ahetch(:), ehetch(:),              &
                                     ahetc2(:), ehetc2(:),              &
                                     ahetwt(:), ehetwt(:)

  integer    nsolim
  parameter( nsolim = 4*ncharm )

  integer, save ::          nsolid

  integer(c_int), pointer, save :: ich(:), ick(:), iash(:), iwat(:)
  real(c_double), pointer, save :: ehsoli(:,:), wmols(:), eh0sol(:)

  ! By class (deduced quantities)

  ! Number of classes
  integer(c_int), pointer, save :: nclacp

  integer(c_int), pointer, save :: ichcor(:)
  real(c_double), pointer, save :: diam20(:), dia2mn(:),                  &
                                   rho20 (:), rho2mn(:),                  &
                                   xmp0(:),   xmash (:)

  !--> Donnees relatives a la combustion des especes gazeuses

  integer(c_int), pointer, save :: ico, ico2, ih2o, io2, in2

  integer(c_int), pointer, save :: ichx1c(:), ichx2c(:), ichx1, ichx2

  real(c_double), pointer, save :: chx1(:), chx2(:),                         &
                                   a1(:), b1(:),c1(:),d1(:), e1(:), f1(:),   &
                                   a2(:), b2(:),c2(:),d2(:), e2(:), f2(:)

  ! Complement Table

  real(c_double), pointer, save :: thc(:)
  integer(c_int), pointer, save :: npoc

  !--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE

  !> \defgroup coal_combustion  Pulverized coal combustion variables

  !> \addtogroup coal_combustion
  !> \{

  ! ---- Variables transportees
  !        Phase continue (melange gazeux)

  !> mean value of the tracer 1 representing the light
  !> volatiles released by the coal \c icha
  integer, save :: if1m(ncharm)

  !> mean value of the tracer 2 representing the heavy
  !> volatiles released by the coal \c icha
  integer, save :: if2m(ncharm)

  !> tracer 4: mass of the oxydant 2 divided by the mass of bulk
  integer, save :: if4m
  !> tracer 5: mass of the oxydant 3 divided by the mass of bulk
  integer, save :: if5m
  !> tracer 6: water coming from drying
  integer, save :: if6m
  !> tracer 7: mass of the carbon from coal oxydized by O2
  !> divided by the mass of bulk
  integer, save :: if7m
  !> tracer 8: mass of the carbon from coal gasified by CO2
  !> divided by the mass of bulk
  integer, save :: if8m
  !> tracer 9: mass of the Carbon from coal gasified by H2O
  !> divided by the mass of bulk
  integer, save :: if9m

  !> f1f2 variance
  integer, save :: ifvp2m

  !        Phase dispersee (classe de particules)
  !> coke mass fraction related to the class icla
  integer, save :: ixck(nclcpm)

  !> reactive coal mass fraction related to the class \c icla
  integer, save :: ixch(nclcpm)

  !> number of particles of the class \c icla per kg of air-coal mixture
  integer, save :: inp(nclcpm)

  !>  mass enthalpy of the coal of class \c icla, if we are in permeatic conditions
  integer, save :: ih2(nclcpm)

  ! TODO absent de la doc utilisateur
  !> transported variable of dispersed phase (particle class)
  integer, save :: ixwt(nclcpm)

  ! ---- Variables d'etat
  !        Phase continue (melange gazeux)

  !> mass fractions:
  !>  - iym1(1): mass fraction of \f$CH_{X1m}\f$ (light volatiles) in the gas mixture
  !>  - iym1(2): mass fraction of \f$CH_{X2m}\f$ (heavy volatiles) in the gas mixture
  !>  - iym1(3): mass fraction of CO in the gas mixture
  !>  - iym1(4): mass fraction of \f$O_2\f$ in the gas mixture
  !>  - iym1(5): mass fraction of \f$CO_2\f$ in the gas mixture
  !>  - iym1(6): mass fraction of \f$H_2O\f$ in the gas mixture
  !>  - iym1(7): mass fraction of \f$N_2\f$ in the gas mixture
  integer, save :: iym1(ngazem)

  ! TODO absent de la doc utilisateur
  !> State variables of continuous phase (gas mixture)
  integer, save :: irom1

  !>  molar mass of the gas mixture
  integer, save :: immel

  !        Phase dispersee (classes de particules)

  !> temperature of the particles of the class \c icla
  integer, save :: itemp2(nclcpm)

  !> density of the particles of the class \c icla
  integer, save :: irom2(nclcpm)

  !> diameter of the particles of the class \c icla
  integer, save :: idiam2(nclcpm)

  !>  solid mass fraction of the class \c icla
  integer, save :: ix2(nclcpm)

  !> disappearance rate of the reactive coal of the class \c icla
  integer, save :: igmdch(nclcpm)

  !> coke disappearance rate of the coke burnout of the class \c icla
  integer, save :: igmhet(nclcpm)

  !> Implicite part of the exchanges to the gas by molecular distribution
  integer, save :: igmtr(nclcpm)

  ! TODO absent de la doc utilisateur
  !> State variables of dispersed phase (particles class)
  integer, save :: ighco2(nclcpm)

  !>  mass transfer caused by the release of light volatiles  of the class \c icla
  integer, save :: igmdv1(nclcpm)

  !>  mass transfer caused by the release of heavy volatiles  of the class \c icla
  integer, save :: igmdv2(nclcpm)

  ! TODO absent de la doc utilisateur
  !> State variables of dispersed phase (particles class)
  integer, save :: igmsec(nclcpm)

  !> Used for bulk balance of Carbon
  integer, save :: ibcarbone
  !> Used for bulk balance of Oxygen
  integer, save :: iboxygen
  !> Used for bulk balance of Hydrogen
  integer, save :: ibhydrogen

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_coal_get_pointers(p_ncharb, p_nclacp,               &
                                             p_nclpch, p_idrift,               &
                                             p_ich, p_ick, p_iash, p_iwat,     &
                                             p_ehsoli, p_wmols, p_eh0sol,      &
                                             p_ichcor, p_diam20, p_dia2mn,     &
                                             p_rho20, p_rho2mn,                &
                                             p_xmp0, p_xmash)                  &
      bind(C, name='cs_f_cpincl_coal_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ncharb, p_nclacp,                          &
                                  p_nclpch, p_idrift,                          &
                                  p_ich, p_ick, p_iash, p_iwat,                &
                                  p_ehsoli, p_wmols, p_eh0sol,                 &
                                  p_ichcor, p_diam20, p_dia2mn,                &
                                  p_rho20, p_rho2mn,                           &
                                  p_xmp0, p_xmash
    end subroutine cs_f_cpincl_coal_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_1(                                     &
         p_cch, p_hch, p_och, p_sch, p_nch, p_alpha, p_beta, p_teta, p_omega,  &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_gamma, p_delta, p_kappa, p_zeta, p_rhock, p_pcick, p_cpashc,        &
         p_h0ashc, p_h02ch, p_cp2wat, p_crepn1, p_crepn2, p_cp2ch,             &
         p_xashsec, p_xashch, p_xwatch)                                        &
      bind(C, name='cs_f_cpincl_get_pointers_1')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_cch, p_hch, p_och, p_sch, p_nch, p_alpha, p_beta, p_teta, p_omega,  &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_gamma, p_delta, p_kappa, p_zeta, p_rhock, p_pcick, p_cpashc,        &
         p_h0ashc, p_h02ch, p_cp2wat, p_crepn1, p_crepn2, p_cp2ch,             &
         p_xashsec, p_xashch, p_xwatch
    end subroutine cs_f_cpincl_get_pointers_1

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_2(                                     &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt)                               &
      bind(C, name='cs_f_cpincl_get_pointers_2')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt
    end subroutine cs_f_cpincl_get_pointers_2

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_3(                                     &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2, p_thc, p_npoc)                    &
      bind(C, name='cs_f_cpincl_get_pointers_3')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2, p_thc, p_npoc
    end subroutine cs_f_cpincl_get_pointers_3

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine cp_model_map_coal() &
    bind(C, name='cs_f_cp_model_map_coal')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ncharb, p_nclacp, p_nclpch, p_idrift,          &
                   p_ich, p_ick, p_iash, p_iwat,                    &
                   p_ehsoli, p_wmols, p_eh0sol, p_ichcor,           &
                   p_diam20, p_dia2mn, p_rho20, p_rho2mn,           &
                   p_xmp0, p_xmash

    type(c_ptr) ::                                                             &
         p_cch, p_hch, p_och, p_sch, p_nch, p_alpha, p_beta, p_teta, p_omega,  &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_gamma, p_delta, p_kappa, p_zeta, p_rhock, p_pcick, p_cpashc,        &
         p_h0ashc, p_h02ch, p_cp2wat, p_crepn1, p_crepn2, p_cp2ch,             &
         p_xashsec, p_xashch, p_xwatch

    type(c_ptr)::                                                              &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt

    type(c_ptr) ::                                                             &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2, p_thc, p_npoc

    call cs_f_cpincl_coal_get_pointers(p_ncharb, p_nclacp,               &
                                       p_nclpch, p_idrift,               &
                                       p_ich, p_ick, p_iash, p_iwat,     &
                                       p_ehsoli, p_wmols, p_eh0sol,      &
                                       p_ichcor, p_diam20, p_dia2mn,     &
                                       p_rho20, p_rho2mn,                &
                                       p_xmp0, p_xmash)

    call c_f_pointer(p_ncharb, ncharb)
    call c_f_pointer(p_nclacp, nclacp)

    call c_f_pointer(p_nclpch, nclpch, [ncharm])

    call c_f_pointer(p_idrift, i_comb_drift)

    call c_f_pointer(p_ich, ich, [ncharm])
    call c_f_pointer(p_ick, ick, [ncharm])
    call c_f_pointer(p_iash, iash, [ncharm])
    call c_f_pointer(p_iwat, iwat, [ncharm])
    call c_f_pointer(p_ehsoli, ehsoli, [nsolim, npot])
    call c_f_pointer(p_wmols, wmols, [nsolim])
    call c_f_pointer(p_eh0sol, eh0sol, [nsolim])
    call c_f_pointer(p_eh0sol, eh0sol, [nsolim])

    call c_f_pointer(p_ichcor, ichcor, [nclcpm])

    call c_f_pointer(p_diam20, diam20, [nclcpm])
    call c_f_pointer(p_dia2mn, dia2mn, [nclcpm])
    call c_f_pointer(p_rho20,  rho20,  [nclcpm])
    call c_f_pointer(p_rho2mn, rho2mn, [nclcpm])
    call c_f_pointer(p_xmp0,   xmp0,   [nclcpm])
    call c_f_pointer(p_xmash,  xmash,  [nclcpm])

    call cs_f_cpincl_get_pointers_1(                                           &
         p_cch, p_hch, p_och, p_sch, p_nch, p_alpha, p_beta, p_teta, p_omega,  &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_gamma, p_delta, p_kappa, p_zeta, p_rhock, p_pcick, p_cpashc,        &
         p_h0ashc, p_h02ch, p_cp2wat, p_crepn1, p_crepn2, p_cp2ch,             &
         p_xashsec, p_xashch, p_xwatch)

    call c_f_pointer(p_cch, cch, [ncharm])
    call c_f_pointer(p_hch, hch, [ncharm])
    call c_f_pointer(p_och, och, [ncharm])
    call c_f_pointer(p_sch, sch, [ncharm])
    call c_f_pointer(p_nch, nch, [ncharm])
    call c_f_pointer(p_alpha, alpha, [ncharm])
    call c_f_pointer(p_beta, beta, [ncharm])
    call c_f_pointer(p_teta, teta, [ncharm])
    call c_f_pointer(p_omega, omega, [ncharm])
    call c_f_pointer(p_pcich, pcich, [ncharm])
    call c_f_pointer(p_rho0ch, rho0ch, [ncharm])
    call c_f_pointer(p_thcdch, thcdch, [ncharm])
    call c_f_pointer(p_cck, cck, [ncharm])
    call c_f_pointer(p_hck, hck, [ncharm])
    call c_f_pointer(p_ock, ock, [ncharm])
    call c_f_pointer(p_sck, sck, [ncharm])
    call c_f_pointer(p_nck, nck, [ncharm])
    call c_f_pointer(p_gamma, gamma, [ncharm])
    call c_f_pointer(p_delta, delta, [ncharm])
    call c_f_pointer(p_kappa, kappa, [ncharm])
    call c_f_pointer(p_zeta, zeta, [ncharm])
    call c_f_pointer(p_pcick, pcick, [ncharm])
    call c_f_pointer(p_rhock, rhock, [ncharm])
    call c_f_pointer(p_cpashc, cpashc, [ncharm])
    call c_f_pointer(p_h0ashc, h0ashc, [ncharm])
    call c_f_pointer(p_h02ch, h02ch, [ncharm])
    call c_f_pointer(p_cp2wat, cp2wat, [ncharm])
    call c_f_pointer(p_crepn1, crepn1, [2, ncharm])
    call c_f_pointer(p_crepn2, crepn2, [2, ncharm])
    call c_f_pointer(p_cp2ch, cp2ch, [ncharm])
    call c_f_pointer(p_xashsec, xashsec, [ncharm])
    call c_f_pointer(p_xashch, xashch, [ncharm])
    call c_f_pointer(p_xwatch, xwatch, [ncharm])

    call cs_f_cpincl_get_pointers_2(                                           &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt)

    call c_f_pointer(p_iy1ch, iy1ch, [ncharm])
    call c_f_pointer(p_iy2ch, iy2ch, [ncharm])
    call c_f_pointer(p_iochet, iochet, [ncharm])
    call c_f_pointer(p_ioetc2, ioetc2, [ncharm])
    call c_f_pointer(p_ioetwt, ioetwt, [ncharm])

    call c_f_pointer(p_y1ch, y1ch, [ncharm])
    call c_f_pointer(p_a1ch, a1ch, [ncharm])
    call c_f_pointer(p_e1ch, e1ch, [ncharm])
    call c_f_pointer(p_y2ch, y2ch, [ncharm])
    call c_f_pointer(p_a2ch, a2ch, [ncharm])
    call c_f_pointer(p_e2ch, e2ch, [ncharm])
    call c_f_pointer(p_ahetch, ahetch, [ncharm])
    call c_f_pointer(p_ehetch, ehetch, [ncharm])
    call c_f_pointer(p_ahetc2, ahetc2, [ncharm])
    call c_f_pointer(p_ehetc2, ehetc2, [ncharm])
    call c_f_pointer(p_ahetwt, ahetwt, [ncharm])
    call c_f_pointer(p_ehetwt, ehetwt, [ncharm])

    call cs_f_cpincl_get_pointers_3(                                           &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2, p_thc, p_npoc)

    call c_f_pointer(p_ico, ico)
    call c_f_pointer(p_ico2, ico2)
    call c_f_pointer(p_ih2o, ih2o)
    call c_f_pointer(p_io2, io2)
    call c_f_pointer(p_in2, in2)
    call c_f_pointer(p_ichx1c, ichx1c, [ncharm])
    call c_f_pointer(p_ichx2c, ichx2c, [ncharm])
    call c_f_pointer(p_ichx1, ichx1)
    call c_f_pointer(p_ichx2, ichx2)

    call c_f_pointer(p_ichx2c, ichx2c, [ncharm])

    call c_f_pointer(p_chx1, chx1, [ncharm])
    call c_f_pointer(p_chx2, chx2, [ncharm])
    call c_f_pointer(p_a1, a1, [ncharm])
    call c_f_pointer(p_b1, b1, [ncharm])
    call c_f_pointer(p_c1, c1, [ncharm])
    call c_f_pointer(p_d1, d1, [ncharm])
    call c_f_pointer(p_e1, e1, [ncharm])
    call c_f_pointer(p_f1, f1, [ncharm])
    call c_f_pointer(p_a2, a2, [ncharm])
    call c_f_pointer(p_b2, b2, [ncharm])
    call c_f_pointer(p_c2, c2, [ncharm])
    call c_f_pointer(p_d2, d2, [ncharm])
    call c_f_pointer(p_e2, e2, [ncharm])
    call c_f_pointer(p_f2, f2, [ncharm])

    call c_f_pointer(p_thc, thc, [npot])
    call c_f_pointer(p_npoc, npoc)

  end subroutine cp_model_map_coal

  !=============================================================================

end module cpincl

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

  !> maximal number of tabulation points
  integer    npotcp
  parameter( npotcp = 8 )

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
                                  pcich(:), rho0ch(:), thcdch(:),            &
                                  cck(:), hck(:), ock(:), sck(:), nck(:),    &
                                  rhock(:), pcick(:),                        &
                                  cpashc(:), h0ashc(:), h02ch(:),            &
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

  integer(c_int), pointer, save :: nsolid

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

  ! POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE

  ! ---- Variables transportees
  !        Phase continue (melange gazeux)

  integer(c_int), pointer, save :: ihgas
  integer(c_int), pointer, save :: if1m(:)
  integer(c_int), pointer, save :: if2m(:)
  integer(c_int), pointer, save :: if4m
  integer(c_int), pointer, save :: if5m
  integer(c_int), pointer, save :: if6m
  integer(c_int), pointer, save :: if7m
  integer(c_int), pointer, save :: if8m
  integer(c_int), pointer, save :: if9m
  integer(c_int), pointer, save :: ifvp2m
  integer(c_int), pointer, save :: ixck(:)
  integer(c_int), pointer, save :: ixch(:)
  integer(c_int), pointer, save :: inp(:)
  integer(c_int), pointer, save :: ih2(:)
  integer(c_int), pointer, save :: ixwt(:)

  ! ---- Variables d'etat
  !        Phase continue (melange gazeux)

  integer(c_int), pointer, save :: iym1(:)
  integer(c_int), pointer, save :: irom1
  integer(c_int), pointer, save :: immel

  !        Phase dispersee (classes de particules)

  integer(c_int), pointer, save :: itemp2(:)
  integer(c_int), pointer, save :: irom2(:)
  integer(c_int), pointer, save :: idiam2(:)
  integer(c_int), pointer, save :: ix2(:)
  integer(c_int), pointer, save :: igmdch(:)
  integer(c_int), pointer, save :: igmhet(:)
  integer(c_int), pointer, save :: igmtr(:)
  integer(c_int), pointer, save :: ighco2(:)
  integer(c_int), pointer, save :: igmdv1(:)
  integer(c_int), pointer, save :: igmdv2(:)
  integer(c_int), pointer, save :: igmsec(:)
  integer(c_int), pointer, save :: ibcarbone
  integer(c_int), pointer, save :: iboxygen
  integer(c_int), pointer, save :: ibhydrogen

  ! Moved from ppthch

  !> engaze(ij) is the massic enthalpy (J/kg) of the i-th elementary gas component
  !> at temperature  th[j]
  real(c_double), pointer, save ::  ehgaze(:,:)

  ! Moved from ppcfu

  ! prise en compte H2  , H2S , SO2 , HCN , NH3
  integer(c_int), pointer, save :: ihy, ih2s, iso2, ihcn, inh3

  integer(c_int), pointer, save :: noxyd

  ! nb de moles de I dans J

  real(c_double), pointer, save :: af3(:),af4(:),af5(:),af6(:)
  real(c_double), pointer, save :: af7(:),af8(:),af9(:)

  ! Equation sur NOX :
  !
  integer(c_int), pointer, save :: ieqnox, imdnox, irb

  ! Combustion heterogene avec le  CO2

  integer(c_int), pointer, save :: ihtco2, ieqco2, iyco2

  !   Scalaires supplementaires : fraction massique de H2, HCN et NO

  integer(c_int), pointer, save ::  iyhcn, iyno, iynh3, ihox

  !   Propriétés supplementaires :
  integer(c_int), pointer, save ::  ighcn1, ighcn2, ignoth, ignh31, ignh32
  !
  !   Affichage des termes source:
  integer(c_int), pointer, save ::  ifhcnd, ifhcnc, ifnh3d, ifnh3c
  integer(c_int), pointer, save ::  ifnohc, ifnonh, ifnoch, ifnoth, ifhcnr
  integer(c_int), pointer, save ::  icnohc, icnonh, icnorb
  !
  ! Constante cinetique (Model de Chen)
  integer(c_int), pointer, save :: igrb

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_coal_get_pointers(p_ncharb, p_nclacp,               &
                                             p_nclpch, p_idrift, p_nsolid,     &
                                             p_ich, p_ick, p_iash, p_iwat,     &
                                             p_ehsoli, p_wmols, p_eh0sol,      &
                                             p_ichcor, p_diam20, p_dia2mn,     &
                                             p_rho20, p_rho2mn,                &
                                             p_xmp0, p_xmash)                  &
      bind(C, name='cs_f_cpincl_coal_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ncharb, p_nclacp,                          &
                                  p_nclpch, p_idrift, p_nsolid,                &
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
         p_cch, p_hch, p_och, p_sch, p_nch,                                    &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_rhock, p_pcick, p_cpashc,                                           &
         p_h0ashc, p_h02ch, p_crepn1, p_crepn2, p_cp2ch,                       &
         p_xashsec, p_xashch, p_xwatch)                                        &
      bind(C, name='cs_f_cpincl_get_pointers_1')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_cch, p_hch, p_och, p_sch, p_nch,                                    &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_rhock, p_pcick, p_cpashc,                                           &
         p_h0ashc, p_h02ch, p_crepn1, p_crepn2, p_cp2ch,                       &
         p_xashsec, p_xashch, p_xwatch
    end subroutine cs_f_cpincl_get_pointers_1

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_2(                                     &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt, p_ehgaze)                     &
      bind(C, name='cs_f_cpincl_get_pointers_2')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt, p_ehgaze
    end subroutine cs_f_cpincl_get_pointers_2

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_3(                                     &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2)                                   &
      bind(C, name='cs_f_cpincl_get_pointers_3')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2
    end subroutine cs_f_cpincl_get_pointers_3

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_4(                                     &
         p_ihgas,                                                              &
         p_if1m, p_if2m, p_if4m, p_if5m, p_if6m, p_if7m, p_if8m, p_if9m,       &
         p_ifvp2m, p_ixck, p_ixch, p_inp, p_ih2, p_ixwt, p_iym1, p_irom1,      &
         p_immel, p_itemp2, p_irom2, p_idiam2, p_ix2, p_igmdch, p_igmhet,      &
         p_igmtr, p_ighco2, p_igmdv1, p_igmdv2, p_igmsec,                      &
         p_ibcarbone, p_iboxygen, p_ibhydrogen)                                &
      bind(C, name='cs_f_cpincl_get_pointers_4')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_ihgas,                                                              &
         p_if1m, p_if2m, p_if4m, p_if5m, p_if6m, p_if7m, p_if8m, p_if9m,       &
         p_ifvp2m, p_ixck, p_ixch, p_inp, p_ih2, p_ixwt, p_iym1, p_irom1,      &
         p_immel, p_itemp2, p_irom2, p_idiam2, p_ix2, p_igmdch, p_igmhet,      &
         p_igmtr, p_ighco2, p_igmdv1, p_igmdv2, p_igmsec,                      &
         p_ibcarbone, p_iboxygen, p_ibhydrogen
    end subroutine cs_f_cpincl_get_pointers_4

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers_5(                                     &
         p_af3, p_af4, p_af5, p_af6, p_af7, p_af8, p_af9,                      &
         p_ihy, p_ih2s, p_iso2, p_ihcn, p_inh3,                                &
         p_ihtco2, p_ieqco2, p_iyco2, p_ieqnox, p_imdnox, p_irb,               &
         p_iyhcn, p_iyno, p_iynh3, p_ihox, p_igrb, p_noxyd,                    &
         p_ighcn1, p_ighcn2, p_ignoth, p_ignh31, p_ignh32,                     &
         p_ifhcnd, p_ifhcnc, p_ifnh3d, p_ifnh3c, p_ifnohc, p_ifnonh, p_ifnoch, &
         p_ifnoth, p_ifhcnr, p_icnohc, p_icnonh, p_icnorb)                                     &
      bind(C, name='cs_f_cpincl_get_pointers_5')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) ::                                              &
         p_af3, p_af4, p_af5, p_af6, p_af7, p_af8, p_af9,                      &
         p_ihy, p_ih2s, p_iso2, p_ihcn, p_inh3,                                &
         p_ihtco2,  p_ieqco2, p_iyco2, p_ieqnox, p_imdnox, p_irb,              &
         p_iyhcn, p_iyno, p_iynh3, p_ihox, p_igrb, p_noxyd,                    &
         p_ighcn1, p_ighcn2, p_ignoth, p_ignh31, p_ignh32,                     &
         p_ifhcnd, p_ifhcnc, p_ifnh3d, p_ifnh3c, p_ifnohc, p_ifnonh, p_ifnoch, &
         p_ifnoth, p_ifhcnr, p_icnohc, p_icnonh, p_icnorb

    end subroutine cs_f_cpincl_get_pointers_5

    ! Interface to C function
    ! Defines the source terms for scalars which are part of
    ! specific physics models. Source terms are defined over one time step.

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
                   p_nsolid, p_ich, p_ick, p_iash, p_iwat,          &
                   p_ehsoli, p_wmols, p_eh0sol, p_ichcor,           &
                   p_diam20, p_dia2mn, p_rho20, p_rho2mn,           &
                   p_xmp0, p_xmash

    type(c_ptr) ::                                                             &
         p_cch, p_hch, p_och, p_sch, p_nch,                                    &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_rhock, p_pcick, p_cpashc,                                           &
         p_h0ashc, p_h02ch, p_crepn1, p_crepn2, p_cp2ch,                       &
         p_xashsec, p_xashch, p_xwatch

    type(c_ptr)::                                                              &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt, p_ehgaze

    type(c_ptr) ::                                                             &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2

    type(c_ptr) ::                                                             &
         p_ihgas,                                                              &
         p_if1m, p_if2m, p_if4m, p_if5m, p_if6m, p_if7m, p_if8m, p_if9m,       &
         p_ifvp2m, p_ixck, p_ixch, p_inp, p_ih2, p_ixwt, p_iym1, p_irom1,      &
         p_immel, p_itemp2, p_irom2, p_idiam2, p_ix2, p_igmdch, p_igmhet,      &
         p_igmtr, p_ighco2, p_igmdv1, p_igmdv2, p_igmsec,                      &
         p_ibcarbone, p_iboxygen, p_ibhydrogen

    type(c_ptr) ::                                                             &
         p_af3, p_af4, p_af5, p_af6, p_af7, p_af8, p_af9,                      &
         p_ihy, p_ih2s, p_iso2, p_ihcn, p_inh3,                                &
         p_ihtco2,  p_ieqco2, p_iyco2, p_ieqnox, p_imdnox, p_irb,              &
         p_iyhcn, p_iyno, p_iynh3, p_ihox, p_igrb, p_noxyd,                    &
         p_ighcn1, p_ighcn2, p_ignoth, p_ignh31, p_ignh32,                     &
         p_ifhcnd, p_ifhcnc, p_ifnh3d, p_ifnh3c, p_ifnohc, p_ifnonh, p_ifnoch, &
         p_ifnoth, p_ifhcnr, p_icnohc, p_icnonh, p_icnorb

    call cs_f_cpincl_coal_get_pointers(p_ncharb, p_nclacp,               &
                                       p_nclpch, p_idrift, p_nsolid,     &
                                       p_ich, p_ick, p_iash, p_iwat,     &
                                       p_ehsoli, p_wmols, p_eh0sol,      &
                                       p_ichcor, p_diam20, p_dia2mn,     &
                                       p_rho20, p_rho2mn,                &
                                       p_xmp0, p_xmash)

    call c_f_pointer(p_ncharb, ncharb)
    call c_f_pointer(p_nclacp, nclacp)
    call c_f_pointer(p_nsolid, nsolid)

    call c_f_pointer(p_nclpch, nclpch, [ncharm])

    call c_f_pointer(p_idrift, i_comb_drift)

    call c_f_pointer(p_ich, ich, [ncharm])
    call c_f_pointer(p_ick, ick, [ncharm])
    call c_f_pointer(p_iash, iash, [ncharm])
    call c_f_pointer(p_iwat, iwat, [ncharm])
    call c_f_pointer(p_ehsoli, ehsoli, [nsolim, npotcp])
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
         p_cch, p_hch, p_och, p_sch, p_nch,                                    &
         p_pcich, p_rho0ch, p_thcdch, p_cck, p_hck, p_ock, p_sck, p_nck,       &
         p_rhock, p_pcick, p_cpashc,                                           &
         p_h0ashc, p_h02ch, p_crepn1, p_crepn2, p_cp2ch,                       &
         p_xashsec, p_xashch, p_xwatch)

    call c_f_pointer(p_cch, cch, [ncharm])
    call c_f_pointer(p_hch, hch, [ncharm])
    call c_f_pointer(p_och, och, [ncharm])
    call c_f_pointer(p_sch, sch, [ncharm])
    call c_f_pointer(p_nch, nch, [ncharm])
    call c_f_pointer(p_pcich, pcich, [ncharm])
    call c_f_pointer(p_rho0ch, rho0ch, [ncharm])
    call c_f_pointer(p_thcdch, thcdch, [ncharm])
    call c_f_pointer(p_cck, cck, [ncharm])
    call c_f_pointer(p_hck, hck, [ncharm])
    call c_f_pointer(p_ock, ock, [ncharm])
    call c_f_pointer(p_sck, sck, [ncharm])
    call c_f_pointer(p_nck, nck, [ncharm])
    call c_f_pointer(p_pcick, pcick, [ncharm])
    call c_f_pointer(p_rhock, rhock, [ncharm])
    call c_f_pointer(p_cpashc, cpashc, [ncharm])
    call c_f_pointer(p_h0ashc, h0ashc, [ncharm])
    call c_f_pointer(p_h02ch, h02ch, [ncharm])
    call c_f_pointer(p_crepn1, crepn1, [2, ncharm])
    call c_f_pointer(p_crepn2, crepn2, [2, ncharm])
    call c_f_pointer(p_cp2ch, cp2ch, [ncharm])
    call c_f_pointer(p_xashsec, xashsec, [ncharm])
    call c_f_pointer(p_xashch, xashch, [ncharm])
    call c_f_pointer(p_xwatch, xwatch, [ncharm])

    call cs_f_cpincl_get_pointers_2(                                           &
         p_iy1ch, p_iy2ch, p_iochet, p_ioetc2, p_ioetwt,                       &
         p_y1ch, p_a1ch, p_e1ch, p_y2ch, p_a2ch, p_e2ch, p_ahetch, p_ehetch,   &
         p_ahetc2, p_ehetc2, p_ahetwt, p_ehetwt, p_ehgaze)

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

    call c_f_pointer(p_ehgaze, ehgaze, [ngazem, npotcp])

    call cs_f_cpincl_get_pointers_3(                                           &
         p_ico, p_ico2, p_ih2o, p_io2, p_in2, p_ichx1c, p_ichx2c,              &
         p_ichx1, p_ichx2, p_chx1, p_chx2, p_a1, p_b1, p_c1, p_d1, p_e1, p_f1, &
         p_a2, p_b2, p_c2, p_d2, p_e2, p_f2)

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

    call cs_f_cpincl_get_pointers_4(                                           &
         p_ihgas,                                                              &
         p_if1m, p_if2m, p_if4m, p_if5m, p_if6m, p_if7m, p_if8m, p_if9m,       &
         p_ifvp2m, p_ixck, p_ixch, p_inp, p_ih2, p_ixwt, p_iym1, p_irom1,      &
         p_immel, p_itemp2, p_irom2, p_idiam2, p_ix2, p_igmdch, p_igmhet,      &
         p_igmtr, p_ighco2, p_igmdv1, p_igmdv2, p_igmsec,                      &
         p_ibcarbone, p_iboxygen, p_ibhydrogen)

    call c_f_pointer(p_ihgas, ihgas)
    call c_f_pointer(p_if1m, if1m, [ncharm])
    call c_f_pointer(p_if2m, if2m, [ncharm])
    call c_f_pointer(p_if4m, if4m)
    call c_f_pointer(p_if5m, if5m)
    call c_f_pointer(p_if6m, if6m)
    call c_f_pointer(p_if7m, if7m)
    call c_f_pointer(p_if8m, if8m)
    call c_f_pointer(p_if9m, if9m)
    call c_f_pointer(p_ifvp2m, ifvp2m)
    call c_f_pointer(p_ixck, ixck, [nclcpm])
    call c_f_pointer(p_ixch, ixch, [nclcpm])
    call c_f_pointer(p_inp, inp, [nclcpm])
    call c_f_pointer(p_ih2, ih2, [nclcpm])
    call c_f_pointer(p_ixwt, ixwt, [nclcpm])
    call c_f_pointer(p_iym1, iym1, [ngazem])
    call c_f_pointer(p_irom1, irom1)
    call c_f_pointer(p_immel, immel)
    call c_f_pointer(p_itemp2, itemp2, [nclcpm])
    call c_f_pointer(p_irom2, irom2, [nclcpm])
    call c_f_pointer(p_idiam2, idiam2, [nclcpm])
    call c_f_pointer(p_ix2, ix2, [nclcpm])
    call c_f_pointer(p_igmdch, igmdch, [nclcpm])
    call c_f_pointer(p_igmhet, igmhet, [nclcpm])
    call c_f_pointer(p_igmtr, igmtr, [nclcpm])
    call c_f_pointer(p_ighco2, ighco2, [nclcpm])
    call c_f_pointer(p_igmdv1, igmdv1, [nclcpm])
    call c_f_pointer(p_igmdv2, igmdv2, [nclcpm])
    call c_f_pointer(p_igmsec, igmsec, [nclcpm])
    call c_f_pointer(p_ibcarbone, ibcarbone)
    call c_f_pointer(p_iboxygen, iboxygen)
    call c_f_pointer(p_ibhydrogen, ibhydrogen)

    call cs_f_cpincl_get_pointers_5(                                          &
         p_af3, p_af4, p_af5, p_af6, p_af7, p_af8, p_af9,                     &
         p_ihy, p_ih2s, p_iso2, p_ihcn, p_inh3,                               &
         p_ihtco2,  p_ieqco2, p_iyco2, p_ieqnox, p_imdnox, p_irb,             &
         p_iyhcn, p_iyno, p_iynh3, p_ihox, p_igrb, p_noxyd,                   &
         p_ighcn1, p_ighcn2, p_ignoth, p_ignh31, p_ignh32,                    &
         p_ifhcnd, p_ifhcnc, p_ifnh3d, p_ifnh3c, p_ifnohc, p_ifnonh, p_ifnoch, &
         p_ifnoth, p_ifhcnr, p_icnohc, p_icnonh, p_icnorb)

    call c_f_pointer(p_af3, af3, [ngazgm])
    call c_f_pointer(p_af4, af4, [ngazgm])
    call c_f_pointer(p_af5, af5, [ngazgm])
    call c_f_pointer(p_af6, af6, [ngazgm])
    call c_f_pointer(p_af7, af7, [ngazgm])
    call c_f_pointer(p_af8, af8, [ngazgm])
    call c_f_pointer(p_af9, af9, [ngazgm])

    call c_f_pointer(p_ihy, ihy)
    call c_f_pointer(p_ih2s, ih2s)
    call c_f_pointer(p_iso2, iso2)
    call c_f_pointer(p_ihcn, ihcn)
    call c_f_pointer(p_inh3, inh3)

    call c_f_pointer(p_ihtco2, ihtco2)
    call c_f_pointer(p_ieqco2, ieqco2)
    call c_f_pointer(p_iyco2, iyco2)

    call c_f_pointer(p_ieqnox, ieqnox)
    call c_f_pointer(p_imdnox, imdnox)
    call c_f_pointer(p_irb, irb)

    call c_f_pointer(p_iyhcn, iyhcn)
    call c_f_pointer(p_iyno, iyno)
    call c_f_pointer(p_iynh3, iynh3)
    call c_f_pointer(p_ihox, ihox)
    call c_f_pointer(p_igrb, igrb)
    call c_f_pointer(p_noxyd, noxyd)

    call c_f_pointer(p_ighcn1, ighcn1)
    call c_f_pointer(p_ighcn2, ighcn2)
    call c_f_pointer(p_ignoth, ignoth)
    call c_f_pointer(p_ignh31, ignh31)
    call c_f_pointer(p_ignh32, ignh32)
    call c_f_pointer(p_ifhcnd, ifhcnd)
    call c_f_pointer(p_ifhcnc, ifhcnc)
    call c_f_pointer(p_ifnh3d, ifnh3d)
    call c_f_pointer(p_ifnh3c, ifnh3c)
    call c_f_pointer(p_ifnohc, ifnohc)
    call c_f_pointer(p_ifnonh, ifnonh)
    call c_f_pointer(p_ifnoch, ifnoch)
    call c_f_pointer(p_ifnoth, ifnoth)
    call c_f_pointer(p_ifhcnr, ifhcnr)
    call c_f_pointer(p_icnohc, icnohc)
    call c_f_pointer(p_icnonh, icnonh)
    call c_f_pointer(p_icnorb, icnorb)

  end subroutine cp_model_map_coal

  !=============================================================================

end module cpincl

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

  !      - Proprietes sur charbon sec
  !        cch(ch)      --> Composition elementaire en C, H, O , S , N sur sec (%)
  !        hch(ch)          du charbon
  !        och(ch)
  !        sch(ch)
  !        nch(ch)
  !        alpha(ch)    --> Composition du charbon reactif
  !        beta(ch)         sous la forme ch(alpha)o(beta)s(gamma)
  !        gamma(ch)        alpha(ch) = hch(ch)/cch(ch)
  !        omega(ch)        beta (ch)  = och(ch)/cch(ch)
  !                         teta (ch)  = sch(ch)/cchl(ch)
  !                         omega(ch)  = nch(ch)/cch(ch)
  !        pcich(ch)    --> pci (J/kg) charbon
  !        rho0ch(ch)   --> Masse volumique initiale (kg/m3)
  !        thcdch(ch)   --> Conductivite thermique du charbon (W/m/K)
  !      - Proprietes sur charbon sec du coke
  !        cck(ch)      --> Composition elementaire en C, H, O , S , N sur sec (%)
  !        hck(ch)          du coke
  !        ock(ch)
  !        sck(ch)
  !        nck(ch)
  !        gamma(ch)    --> Composition du coke
  !        delta(ch)        sous la forme ch(gamma)o(delta)s(kappa)n(zeta)
  !        kappa(ch)        gamma(ch) = hck(ch)/cck(ch)
  !        zeta (ch)        delta(ch) = ock(ch)/cck(ch)
  !                         kappa(ch) = sck(ch)/cck(ch)
  !                         zeta(ch)  = nck(ch)/cck(ch)
  !        pcick(ch)    --> PCI (J/kg) coke
  !        rhock(ch)    --> Masse volumique coke
  !      - Proprietes sur charbon sec des cendres (ou humide)
  !        cpashc(ch)   --> Cp des cendres (J/kg/K)
  !        h0ashc(ch)   --> Enthalpie de formation des cendres (J/kg)
  !        h02ch        --> H0 du Charbon
  !        crepn1(2,ch) --> repartition de l'azote en HCN etNo reaction 1
  !        crepn2(2,ch) --> repartition de l'azote en HCN etNo reaction 2

  double precision, save :: cch   (ncharm), hch   (ncharm), och   (ncharm),  &
                            sch(ncharm)  , nch(ncharm),                      &
                            alpha (ncharm), beta  (ncharm), teta(ncharm)  ,  &
                            omega(ncharm),                                   &
                            pcich (ncharm), rho0ch(ncharm), thcdch(ncharm),  &
                            cck   (ncharm), hck   (ncharm), ock   (ncharm),  &
                            sck (ncharm), nck(ncharm),                       &
                            gamma (ncharm), delta (ncharm), kappa(ncharm),   &
                            zeta(ncharm),                                    &
                            rhock (ncharm), pcick (ncharm),                  &
                            cpashc(ncharm),                                  &
                            h0ashc(ncharm),                                  &
                            h02ch (ncharm),                                  &
                            cp2wat(ncharm),                                  &
                            crepn1(2,ncharm),crepn2(2,ncharm)


  !        cp2ch        --> Cp du Charbon
  !        xashch(ch)   --> Taux de cendre (kg/kg)
  !        xwatch(ch)   --> Taux d'humidite (kg/kg)

  real(c_double), pointer, save :: cp2ch(:), xashch(:), xwatch(:)

  !      - Parametres cinetiques pour la devolatilisation
  !         (Modele de Kobayashi)
  !        iy1ch(ch)    --> Indicateur : 0 si MVl = {CH4;CO}
  !                                      1 si MVl = {CHz;CO}
  !        y1ch(ch)     --> Coefficient stoechiometrique (adim)
  !                         calcule si IY1CH = 0 ; donne si IY1CH = 1
  !        a1ch(ch)     --> Facteur pre-exponetielle (1/s)
  !        e1ch(ch)     --> Energie d'activation (J/mol)
  !        iy2ch(ch)    --> Indicateur : 0 si MVL = {C2H4;CO}
  !                                      1 si MVL = {CxHy;CO}
  !        y2ch(ch)     --> Coefficient stoechiometrique (adim)
  !                         calcule si IY2CH = 0 ; donne si IY2CH = 1
  !        a2ch(ch)     --> Constante preexponetielle (1/s)
  !        e2ch(ch)     --> Energie d'activation (J/mol)

  !        - Parametres cinetiques pour la combustion heterogene du coke avec O2
  !           (Modele a sphere retrecissante)
  !        ahetch(ch)   --> Constante pre-exponentielle (kg/m2/s/atm)
  !        ehetch(ch)   --> Energie d'activation (kcal/mol)
  !        iochet(ch)   --> Ordre de la reaction 0.5 si = 0 1 si = 1

  !        - Parametres cinetiques pour la combustion heterogene du coke avec CO2
  !           (Modele a sphere retrecissante)
  !        ahetc2(ch)   --> Constante pre-exponentielle (kg/m2/s/atm)
  !        ehetc2(ch)   --> Energie d'activation (kcal/mol)
  !        ioetc2(ch)   --> Ordre de la reaction 0.5 si = 0 1 si = 1

  !        - Parametres cinetiques pour la combustion heterogene du coke avec H2O
  !           (Modele a sphere retrecissante)
  !        ahetwt(ch)   --> Constante pre-exponentielle (kg/m2/s/atm)
  !        ehetwt(ch)   --> Energie d'activation (kcal/mol)
  !        ioetwt(ch)   --> Ordre de la reaction 0.5 si = 0 1 si = 1

  integer, save ::          iy1ch (ncharm), iy2ch (ncharm)
  integer, save ::          iochet (ncharm) , ioetc2(ncharm), ioetwt(ncharm)
  double precision, save :: y1ch  (ncharm), a1ch  (ncharm), e1ch  (ncharm),  &
                            y2ch  (ncharm), a2ch  (ncharm), e2ch  (ncharm),  &
                            ahetch(ncharm), ehetch(ncharm),                  &
                            ahetc2(ncharm), ehetc2(ncharm),                  &
                            ahetwt(ncharm), ehetwt(ncharm)

  !      - Enthalpie du charbon reactif, coke et cendres
  !     ich(ch)      --> Pointeur dans le tableau ehsoli pour
  !                         le Charbon Reactif
  !     ick(ch)      --> Pointeur dans le tableau ehsoli pour le Coke
  !     iash(ch)     --> Pointeur dans le tableau ehsoli pour les cendres
  !     iwat(ch)     --> Pointeur dans le tableau ehsoli pour l'humidite
  !     nsolid       --> Nb constituants solides (Ch.Reactif, Coke, Ash)
  !     nsolim       --> Nb maximal de constituants solides
  !     ehsoli(s,it) --> Enthalpie massique (J/kg) du constituant solide
  !                         no S a la temperature T(it)
  !     wmols(s)     --> Masse molaire du constituant solide
  !     eh0sol(s)    --- Enthalpie de formation (J/kg) du constituant solide
  !                      no S

  integer    nsolim
  parameter( nsolim = 4*ncharm )

  integer, save ::          nsolid

  integer(c_int), pointer, save :: ich(:), ick(:), iash(:), iwat(:)
  real(c_double), pointer, save :: ehsoli(:,:), wmols(:), eh0sol(:)

  ! By class (deduced quantities)

  ! Number of classes
  integer(c_int), pointer, save :: nclacp

  !      - Proprietes
  !        ichcor(cl)  --> = ich si la classe consideree appartient
  !                        au charbon ich(1, 2, ...)
  !        diam20(cl)  --> Diametre initial (m)
  !        dia2mn(cl)  --> Diametre minimum (m)
  !        rho20(cl)   --> Masse volumique initiale (kg/m3)
  !        rho2mn(cl)  --> Masse volumique minimale (kg/m3)
  !        xmp0(cl)    --> Masse initiale de la particule (m)
  !        xmash(cl)   --> Masse de cendres de la particule (m)

  integer(c_int), pointer, save :: ichcor(:)
  real(c_double), pointer, save :: diam20(:), dia2mn(:),                  &
                                   rho20 (:), rho2mn(:),                  &
                                   xmp0(:),   xmash (:)

  !--> Donnees relatives a la combustion des especes gazeuses

  !        ichx1c(ch)  --> Pointeur CHx1  pour ehgaze et wmole
  !        ichx2c(ch)  --> Pointeur CHx2  pour ehgaze et wmole
  !        ichx1       --> Pointeur CHx1m pour ehgaze et wmole
  !        ichx2       --> Pointeur CHx2m pour ehgaze et wmole
  !        ico         --> Pointeur CO    pour ehgaze et wmole
  !        io2         --> Pointeur O2    pour ehgaze et wmole
  !        ico2        --> Pointeur CO2   pour ehgaze et wmole
  !        ih2o        --> Pointeur H2O   pour ehgaze et wmole
  !        in2         --> Pointeur N2    pour ehgaze et wmole
  !        chx1(ch)    --> Composition de l'hydrocarbure relatif
  !                        au MVl : CH(X1)
  !        chx2(ch)    --> Composition de l'hydrocarbure relatif
  !                        au MVL : CH(X2)
  !        a1(ch),     --> Coefficients stoechiometriques molaires pour
  !        b1(ch)          la reaction de devolatilisation a basses T
  !        c1(ch)
  !        d1(ch)
  !        e1(ch)
  !        f1(ch)
  !        a2(ch),     --> Coefficients stoechiometriques molaires pour
  !        b2(ch)          la reaction de devolatilisation a basses T
  !        c2(ch)
  !        d2(ch)
  !        e2(ch)
  !        f2(ch)

  integer, save ::          ichx1c(ncharm), ichx2c(ncharm),                  &
                            ichx1, ichx2

  integer(c_int), pointer, save :: ico, ico2, ih2o, io2, in2

  double precision, save :: chx1(ncharm), chx2(ncharm),                      &
                            a1(ncharm), b1(ncharm),c1(ncharm),d1(ncharm),    &
                            e1(ncharm), f1(ncharm),                          &
                            a2(ncharm), b2(ncharm),c2(ncharm),d2(ncharm),    &
                            e2(ncharm), f2(ncharm)

  !--> Pointeurs dans le tableau tbmcr

  integer, save :: if1mc(ncharm) , if2mc(ncharm)
  integer, save :: ix1mc ,ix2mc, ichx1f1, ichx2f2
  integer, save :: icof1, icof2, ih2of1 , ih2of2
  integer, save :: ih2sf1, ih2sf2 , ihcnf1 , ihcnf2

  ! Complement Table

  double precision, save :: thc(npot)
  integer, save ::          npoc

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

  !> Pointer to Np*age(particles)
  integer, save :: inagecp(nclcpm)

  !>
  integer, save :: iv_p_x(nclcpm)

  !>
  integer, save :: iv_p_y(nclcpm)

  !>
  integer, save :: iv_p_z(nclcpm)

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

    subroutine cs_f_cpincl_get_pointers(p_ico, p_ico2, p_ih2o, p_io2, p_in2)    &
      bind(C, name='cs_f_cpincl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ico, p_ico2, p_ih2o, p_io2, p_in2
    end subroutine cs_f_cpincl_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_coal_get_pointers(p_ncharb, p_nclacp,               &
                                             p_nclpch, p_idrift,               &
                                             p_ich, p_ick, p_iash, p_iwat,     &
                                             p_ehsoli, p_wmols, p_eh0sol,      &
                                             p_ichcor, p_cp2ch, p_xashch,      &
                                             p_xwatch, p_diam20, p_dia2mn,     &
                                             p_rho20, p_rho2mn,                &
                                             p_xmp0, p_xmash)                  &
      bind(C, name='cs_f_cpincl_coal_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ncharb, p_nclacp,                          &
                                  p_nclpch, p_idrift,                          &
                                  p_ich, p_ick, p_iash, p_iwat,                &
                                  p_ehsoli, p_wmols, p_eh0sol,                 &
                                  p_ichcor, p_cp2ch, p_xashch,                 &
                                  p_xwatch, p_diam20, p_dia2mn,                &
                                  p_rho20, p_rho2mn,                           &
                                  p_xmp0, p_xmash
    end subroutine cs_f_cpincl_coal_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine cp_models_init() &
    bind(C, name='cs_f_cp_models_init')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ico, p_ico2, p_ih2o, p_io2, p_in2

    call cs_f_cpincl_get_pointers(p_ico, p_ico2, p_ih2o, p_io2, p_in2)

    call c_f_pointer(p_ico, ico)
    call c_f_pointer(p_io2, io2)
    call c_f_pointer(p_in2, in2)
    call c_f_pointer(p_ico2, ico2)
    call c_f_pointer(p_ih2o, ih2o)

  end subroutine cp_models_init

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
                   p_ehsoli, p_wmols, p_eh0sol,                     &
                   p_ichcor, p_cp2ch, p_xashch,                     &
                   p_xwatch, p_diam20, p_dia2mn, p_rho20, p_rho2mn, &
                   p_xmp0, p_xmash

    call cs_f_cpincl_coal_get_pointers(p_ncharb, p_nclacp,               &
                                       p_nclpch, p_idrift,               &
                                       p_ich, p_ick, p_iash, p_iwat,     &
                                       p_ehsoli, p_wmols, p_eh0sol,      &
                                       p_ichcor, p_cp2ch, p_xashch,      &
                                       p_xwatch, p_diam20, p_dia2mn,     &
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

    call c_f_pointer(p_cp2ch, cp2ch, [ncharm])
    call c_f_pointer(p_xashch, xashch, [ncharm])
    call c_f_pointer(p_xwatch, xwatch, [ncharm])

    call c_f_pointer(p_diam20, diam20, [nclcpm])
    call c_f_pointer(p_dia2mn, dia2mn, [nclcpm])
    call c_f_pointer(p_rho20,  rho20,  [nclcpm])
    call c_f_pointer(p_rho2mn, rho2mn, [nclcpm])
    call c_f_pointer(p_xmp0,   xmp0,   [nclcpm])
    call c_f_pointer(p_xmash,  xmash,  [nclcpm])

  end subroutine cp_model_map_coal

  !=============================================================================

end module cpincl

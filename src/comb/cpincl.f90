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
  integer, save :: ncharb

  ! By coal (given quantities)

  ! Granulometric distribution

  !> Number of classes per coal
  integer, save :: nclpch(ncharm)

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
  !        xashch(ch)   --> Taux de cendre (kg/kg)
  !        cpashc(ch)   --> Cp des cendres (J/kg/K)
  !        h0ashc(ch)   --> Enthalpie de formation des cendres (J/kg)
  !        h02ch        --> H0 du Charbon
  !        cpch         --> Cp du Charbon
  !        xwatch(ch)   --> Taux d'humidite (kg/kg)
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
                          h02ch (ncharm), cp2ch (ncharm),                  &
                          xwatch(ncharm), cp2wat(ncharm),                  &
                          crepn1(2,ncharm),crepn2(2,ncharm)

real(c_double), pointer, save :: xashch(:)

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

  integer, save ::          nsolid, ich(ncharm), ick(ncharm), iash(ncharm)
  integer, save ::          iwat(ncharm)
  double precision, save :: ehsoli(nsolim,npot), wmols(nsolim)
  double precision, save :: eh0sol(nsolim)

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
                            ichx1, ichx2, ico, io2, in2

  integer(c_int), pointer, save :: ico2, ih2o

  double precision, save :: chx1(ncharm), chx2(ncharm),                      &
                            a1(ncharm), b1(ncharm),c1(ncharm),d1(ncharm),    &
                            e1(ncharm), f1(ncharm),                          &
                            a2(ncharm), b2(ncharm),c2(ncharm),d2(ncharm),    &
                            e2(ncharm), f2(ncharm)

  !--> Donnees complementaires relatives au calcul de rho
  !    sur les facettes de bord

  !       ientat(ient) --> Indicateur air par type de facette d'entree
  !       ientcp(ient) --> Indicateur Cp  par type de facette d'entree
  !       timpat(ient) --> Temperature en K pour l'air relative
  !                         a l'entree ient
  !       x20(ient,    --> Fraction massique dans le melange de charbon
  !           icla   )     de la classe icla relative a l'entree ient

  integer, save ::          ientat(nozppm), ientcp(nozppm)
  double precision, save :: timpat(nozppm), x20(nozppm,nclcpm)

  !--> Pointeurs dans le tableau tbmcr

  integer, save :: if1mc(ncharm) , if2mc(ncharm)
  integer, save :: ix1mc ,ix2mc, ichx1f1, ichx2f2
  integer, save :: icof1, icof2, ih2of1 , ih2of2
  integer, save :: ih2sf1, ih2sf2 , ihcnf1 , ihcnf2

  !--> Grandeurs fournies par l'utilisateur en conditions aux limites
  !      permettant de calculer automatiquement la vitesse, la turbulence,
  !      l'enthalpie d'entree.

  !    Pour les entrees uniquement, ient etant le numero de zone frontiere

  !       qimpat(ient)           --> Debit       air          en kg/s
  !       timpat(ient)           --> Temperature air          en K
  !       qimpcp(ient,icha)      --> Debit       charbon icha en kg/s
  !       timpcp(ient,icha)      --> Temperature charbon icha en K
  !       distch(ient,icha,icla) --> Distribution en %masse de la classe icla
  !                                  pour le charbon icha

  double precision, save ::  qimpat(nozppm)
  double precision, save ::  qimpcp(nozppm,ncharm), timpcp(nozppm,ncharm)
  double precision, save ::  distch(nozppm,ncharm,ncpcmx)

  ! Complement Table

  double precision, save :: thc(npot)
  integer, save ::          npoc

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_cpincl_get_pointers(p_ico2, p_ih2o, p_nclacp, p_ichcor,    &
                                        p_xashch, p_diam20, p_dia2mn,          &
                                        p_rho20, p_rho2mn,                     &
                                        p_xmp0, p_xmash)                       &
      bind(C, name='cs_f_cpincl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ico2, p_ih2o, p_nclacp, p_ichcor,          &
                                  p_xashch, p_diam20, p_dia2mn,                &
                                  p_rho20, p_rho2mn,                           &
                                  p_xmp0, p_xmash
    end subroutine cs_f_cpincl_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine cp_models_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ico2, p_ih2o, p_nclacp, p_ichcor, p_xashch,             &
                   p_diam20, p_dia2mn, p_rho20, p_rho2mn, p_xmp0, p_xmash

    call cs_f_cpincl_get_pointers(p_ico2, p_ih2o, p_nclacp, p_ichcor,        &
                                  p_xashch, p_diam20, p_dia2mn,              &
                                  p_rho20, p_rho2mn,     &
                                  p_xmp0, p_xmash)

    call c_f_pointer(p_ico2, ico2)
    call c_f_pointer(p_ih2o, ih2o)

    call c_f_pointer(p_nclacp, nclacp)

    call c_f_pointer(p_ichcor, ichcor, [nclcpm])

    call c_f_pointer(p_xashch, xashch, [ncharm])

    call c_f_pointer(p_diam20, diam20, [nclcpm])
    call c_f_pointer(p_dia2mn, dia2mn, [nclcpm])
    call c_f_pointer(p_rho20,  rho20,  [nclcpm])
    call c_f_pointer(p_rho2mn, rho2mn, [nclcpm])
    call c_f_pointer(p_xmp0,   xmp0,   [nclcpm])
    call c_f_pointer(p_xmash,  xmash,  [nclcpm])

  end subroutine cp_models_init

  !=============================================================================

end module cpincl

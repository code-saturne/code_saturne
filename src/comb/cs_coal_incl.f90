!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_coal_incl.f90
!> Module for coal combustion

module cs_coal_incl

  !=============================================================================

  use ppppar

  implicit none

  ! Combustion du coke par H2O

  integer, save :: ihth2o , ighh2o(nclcpm), ipci(ncharm)

  ! Modele de NOx
  ! qpr : % d'azote libere pendant la devol.% de MV libere pendant la devol.
  ! fn : concentration en azote sur pur

  double precision, save :: qpr(ncharm), fn(ncharm), xashsec(ncharm)

  ! Nouveau modele NOx
  ! ==================
  ! ychxle   : Fraction massique de CHx1(icha) dans le MV legeres.
  ! ychxlo   : Fraction massique de CHx1(icha) dans le MV lourdes.
  ! yhcnle   : Fraction massique de HCN dans le MV legeres.
  ! yhcnlo   : Fraction massique de HCN dans les MV lourdes.
  ! ynh3le   : Fraction massique de NH3 dans le MV legeres.
  ! ynh3lo   : Fraction massique de NH3 dans le MV lourdes.
  ! ychxmvlo : Fraction massique de CHx2 dans les MV lourdes.
  ! repnle/lo: Pourcentage de l'azote total du charbon du char1/char2.
  ! repnck   : Pourcentage du HCN libere lors de la combustion heterogene.
  ! ycoch1   : Fraction massique de CO dans les produits de la combustion
  !            heterogene du Char 1.
  ! yhcnc1   : Fraction massique de HCN dans les produits de la combustion
  !            heterogene du Char 1.
  ! ynoch1   : Fraction massique de NO dans les produits de la combustion
  !            heterogene du Char 1.
  ! ycoch2   : Fraction massique de CO dans les produits de la combustion
  !            heterogene du Char 2.
  ! yhcnc2   : Fraction massique de HCN dans les produits de la combustion
  !            heterogene du Char 2.
  ! ynoch2   : Fraction massique de NO dans les produits de la combustion
  !            heterogene du Char 2.
  ! nnch     : Fraction massique d'azote dans le charbon.
  ! nnckle   : Fraction massique d'azote dans le Char 1.
  ! nhckle   : Fraction massique d'hydrogen dans le Char 1.
  ! ncckle   : Fraction massique de carbone dans le Char 1.
  ! nncklo   : Fraction massique d'azote dans le Char 2.
  ! nhcklo   : Fraction massique d'hydrogen dans le Char 2.
  ! nccklo   : Fraction massique de carbone dans le Char 2.
  ! wchx1c   : masse molaire de CHx1 du charbon i
  ! wchx2c   : masse molaire de Chx2 du charbon i
  ! wmchx1   : masse molaire des CHx1
  ! wmchx2   : masse molaire des CHx2

  double precision, save :: repnck(ncharm), repnle(ncharm), repnlo(ncharm)
  double precision, save :: ychxle(ncharm), ychxlo(ncharm)
  double precision, save :: yhcnle(ncharm), yhcnlo(ncharm), ynh3le(ncharm),    &
                            ynh3lo(ncharm)
  double precision, save :: ycoch1(ncharm), yhcnc1(ncharm), ynoch1(ncharm)
  double precision, save :: ycoch2(ncharm), yhcnc2(ncharm), ynoch2(ncharm)
  double precision, save :: nnch(ncharm),   nnckle(ncharm), nhckle(ncharm),    &
                            ncckle(ncharm)
  double precision, save :: nncklo(ncharm), nhcklo(ncharm), nccklo(ncharm)
  double precision, save :: wchx1c(ncharm), wchx2c(ncharm), wmchx1,wmchx2

  !=============================================================================

end module cs_coal_incl

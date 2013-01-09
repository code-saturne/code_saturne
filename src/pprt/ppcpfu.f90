!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file ppcpfu.f90
!> Module for specific physics common variables
!> between combustion of pulverized coal and heavy fuel

module ppcpfu

  !===========================================================================

  use ppppar
  use ppthch

  ! XSI         --> XSI = 3,76 pour de l'air

  double precision, save ::  xsi

  ! nb de moles de I dans J

  double precision, save :: ao2f3,acof3,an2f3,ah2of3
  double precision, save :: ao2f4,an2f4,ah2of4,aco2f4
  double precision, save :: ah2of5
  double precision, save :: ao2f6,an2f6,ah2of6,aco2f6
  double precision, save :: ao2f7,an2f7,ah2of7,aco2f7

  ! nb de moles de I dans J : nouvelle version

  double precision, save :: af3(ngazgm),af4(ngazgm),af5(ngazgm),af6(ngazgm)
  double precision, save :: af7(ngazgm),af8(ngazgm),af9(ngazgm)

  ! prise en compte H2  , H2S , SO2 , HCN , NH3
  integer, save ::         ihy , ih2s , iso2  , ihcn , inh3

  ! Equation sur YCO2

  integer, save ::         ieqco2 , iyco2

  ! Combustion heterogene avec le  CO2

  integer, save ::         ihtco2

  ! Equation sur NOX :
  ! ================

  !   IEQNOX = 0 pas de NOx
  !          = 1 calcul du NOx

  integer, save ::         ieqnox

  !   Scalaires supplementaires : fraction massique de H2, HCN et NO
  !                               temperature air

  integer, save ::         iyhcn , iyno , itaire, ihox

  !   Propce supplementaires :

  !         Conversion HCN en NO       : EXP(-E1/RT)
  !         Conversion HCN en NO       : EXP(-E2/RT)
  !         NO thermique (Zel'dovitch) : EXP(-E3/RT)

  integer, save ::         ighcn1 , ighcn2 , ignoth

  !   Temperature moyenne d'entree
  !   Taux de vapeur moyen

  double precision, save :: taire

  !--> DONNEES RELATIVES AUX OXYDANTS

  !       NOXYD        --> Nombre d'Oxydants (Maximum 3)

  integer, save :: noxyd

  !       OXYO2       --> composition des oxydants en O2
  !       OXYN2       --> composition des oxydants en N2
  !       OXYH2O      --> composition des oxydants en H2O
  !       OXYCO2      --> composition des oxydants en CO2

  double precision, save :: oxyo2(3),oxyn2(3),oxyh2o(3),oxyco2(3)

  !--> Conditions aux limites

  integer, save :: inmoxy(nozppm)

  !=============================================================================

end module ppcpfu

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

  !=============================================================================

  use ppppar
  use ppthch

  implicit none

  !===========================================================================

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
  !
  !   ieqnox = 0     pas de NOx
  !          = 1     calcul du NOx
  !
  !   NOx Model Features:
  !   imdnox = 0 : - HCN is the only intermediate nitrogen species
  !                  liberated during the devolatilisation process.
  !                - HCN is the only intermediate nitrogen species
  !                  liberated during char combustion.
  !                - Constant ratio of the nitrogen mass liberated
  !                  during devolatilisation and the nitrogen mass
  !                  remaining in char.
  !          = 1 (only with iccoal = 0(1)):
  !                - HCN and NH3 are the intermediate nitrogen
  !                  species liberated during the devolatilisation
  !                  process.
  !                - HCN and NO are the intermediate nitrogen species
  !                  liberated during char combustion.
  !                - Temperature dependent ratios of the nitrogen
  !                  mass liberated during devolatilisation and the
  !                  nitrogen mass remaining in char.
  !                - Activation of Reburning kinetics is possible.
  !
  !   irb    = 0     Pas de reburning
  !          = 1     Modele de Chen et al.
  !          = 2     Modele de Dimitriou et al.
  !
  integer, save ::         ieqnox, imdnox, irb
  !
  !   Scalaires supplementaires : fraction massique de H2, HCN et NO
  !                               temperature air

  integer, save ::         iyhcn , iyno , itaire, ihox
  !
  integer, save ::         iynh3
  !
  !
  !   Propriétés supplementaires :

  !         Conversion HCN en NO       : EXP(-E1/RT)
  !         Conversion HCN en NO       : EXP(-E2/RT)
  !         NO thermique (Zel'dovitch) : EXP(-E3/RT)

  integer, save ::         ighcn1 , ighcn2 , ignoth
  !
  !   Propriétés supplementaires :

  !         Conversion NH3 en NO       : EXP(-E4/RT)
  !         Conversion NH3 en NO       : EXP(-E5/RT)
  !
  integer, save ::         ignh31 , ignh32
  !
  !   Affichage des termes source:
  !
  !         ifhcnd: Liberation de HCN au cours de la devolatilisation
  !         ifhcnc: Liberation de HCN au cours de la combustion heterogene
  !         ifnh3d: Liberation de NH3 au cours de la devolatilisation
  !         ifnh3c: Liberation de NH3 au cours de la combustion heterogene
  !         ifnohc: Formation de NO selon la reaction HCN + O2 -> NO + ...
  !         ifnonh: Formation de NO selon la reaction NH3 + O2 -> NO + ...
  !         ifnoch: Liberation de NO au cours de la combustion heterogene
  !         ifnoth: Formation de NO thermique
  !         ifhcnr: Formation de NO selon le mecanisme du "Reburning"
  !         icnohc: Consommation de NO selon la reaction HCN + NO -> Produits
  !         icnonh: Consommation de NO selon la reaction NH3 + NO -> Produits
  !         icnorb: Consommation de NO selon le mecanisme du "Reburning"
  !
  integer, save ::         ifhcnd, ifhcnc, ifnh3d, ifnh3c
  integer, save ::         ifnohc, ifnonh, ifnoch, ifnoth, ifhcnr
  integer, save ::         icnohc, icnonh, icnorb
  !
  !
  !   Temperature moyenne d'entree
  !   Taux de vapeur moyen

  double precision, save :: taire
  !
  !   Reburning
  !
  !   Tableau de la temperature utilise pour determiner la cinetique du
  !   "Reburning".
  double precision, save :: teno(npot)
  !
  !   Tableau des constantes cinetiques (Model de Dimitriou)
  double precision, save :: ka(4,npot), kb(4,npot), kc(4,npot), chi2(npot)
  !
  !   Constante cinetique (Model de Chen)
  integer, save          :: igrb
  !
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

  ! The gas enthalpy is transported (x1 h1 to be correct)
  integer, save :: ihgas

  !=============================================================================

end module ppcpfu

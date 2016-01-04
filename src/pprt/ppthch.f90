!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

  !> number of tabulation points
  integer, save ::           npo
  !> number of elementary gas components
  integer, save ::           ngaze
  !> number of global species
  integer, save ::           ngazg
  !> number of atomic species
  integer, save ::           nato

  !> number of global reactions in gas phase
  integer, save ::           nrgaz

  !> temperature (in K)
  double precision, save ::  th(npot)

  !> engaze(ij) is the massic enthalpy (J/kg) of the i-th elementary gas component
  !> at temperature  th(j)
  double precision, save ::  ehgaze(ngazem,npot)

  !> engazg(ij) is the massic enthalpy (J/kg) of the i-th global secies
  !> at temperature  th(j)
  double precision, save ::  ehgazg(ngazgm,npot)

  !> cpgazg(ij) is the massic calorific capacity (J/kg/K) of the i-th global secies
  !> at temperature  th(j)
  double precision, save ::  cpgazg(ngazgm,npot)

  !> molar mass of an elementary gas component
  double precision, save ::  wmole(ngazem)

  !> molar mass of a global species
  double precision, save ::  wmolg(ngazgm)

  !> molar mass of atoms
  double precision, save ::  wmolat(natom)

  !> Stoichiometry in reaction global species.  Negative for the reactants,
  !> and positive for the products
  double precision, save ::  stoeg(ngazgm,nrgazm)

  !> Mixing rate at the stoichiometry
  double precision, save ::  fs(nrgazm)

  !> Absorption coefficient of global species
  double precision, save ::  ckabsg(ngazgm)

  !> Absorption coefficient of gas mixture
  double precision, save ::  ckabs1

  !> molecular diffusivity for the enthalpy (\f$kg.m^{-1}.s^{-1}\f$)
  !> for gas or coal combustion (the code then automatically sets
  !> \ref optcal::visls0 "visls0" to \ref diftl0 for the scalar
  !> representing the enthalpy).
  !>
  !> Always useful for gas or coal combustion.
  double precision, save ::  diftl0

  !> Molar coefficient of CO2
  double precision, save ::  xco2
  !> Molar coefficient of H2O
  double precision, save ::  xh2o

  !=============================================================================

  !> \}

end module ppthch

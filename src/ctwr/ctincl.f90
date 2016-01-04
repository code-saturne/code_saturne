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

! Module for cooling towers

module ctincl

  !=============================================================================

  implicit none

  !=============================================================================

  ! ichrze : activate postprocessing of cooling tower exchange zones

  integer, save :: ichrze

  ! iaeeri : activation de l'ecart impose
  ! iaeerp : frequence de modification de la temperature

  integer, save :: iaeeri, iaeerp, nbzsup, nbzinf

  ! vaeeri : ecart de refrigeration a imposer
  ! paseri : pas de temperature pour le calcul de la pente de ecartref(teau)
  ! aetemn : minimum de la temperature d'eau refroidie moyenne ponderee
  ! aetemx : maximum de la temperature d'eau chaude moyenne ponderee

  double precision, save :: vaeeri, paseri, aetemn, aetemx, inbaei, &
                            lizsup(100), lizinf(100)

  !=============================================================================

end module ctincl

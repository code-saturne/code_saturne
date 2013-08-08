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

!> \file cplsat.f90
!> Module for code/code coupling

module cplsat

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !  couplage code / code - gestion des parametres principaux

  ! nbrcpl : nombre de couplage Code_Saturne / Code_Saturne
  ! ifaccp : indicateur de couplage face/face uniquement
  ! imobil : indicateur de maillage mobile pour les turbomachines

  integer, save :: nbrcpl, ifaccp, imobil

  ! nbcpmx : nombre de couplage max admissible

  integer   nbcpmx
  parameter(nbcpmx=10)

  ! iturcp(nbcpmx) : modele de turbulence de l'instance distante
  ! imajcp(nbcpmx) : indice de mise a jour de la localisation du couplage
  ! icormx(nbcpmx) : indice de presence de calcul en repere relatif
  ! nvarcp(nbcpmx) : nombre de variables a envoyer/recevoir
  ! nvarto(nbcpmx) : taille des tableaux d'echange

  integer, save :: iturcp(nbcpmx), imajcp(nbcpmx), icormx(nbcpmx)
  integer, save :: nvarcp(nbcpmx), nvarto(nbcpmx)

  !> Absolute time value after the mesh starts to rotate (if it does),
  !> for previous calculation
  double precision, save :: ttpmob

  !> Current absolute time after the mesh starts to rotate (if it does).
  !> In case of restart, this is equal to ttpmob + additional computed time.
  double precision, save :: ttcmob

  !=============================================================================

end module cplsat

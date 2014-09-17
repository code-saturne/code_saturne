!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file lagdim.f90
!> Module for Lagrangian dimensions

module lagdim

  !=============================================================================

  implicit none

  !===========================================================================

  !         Trois modules complementaires
  !                            lagran qui porte les non dimensions
  !                            lagdim qui porte les dimensions variables
  !                            lagpar qui porte les parametres

  !===========================================================================
  ! 1. Connectivite

  !     LONGUEUR DU TABLEAU DU CONNECTIVITE CELLULES -> FACES
  !     (calcule dans le sous-programme LAGINI)

  integer, save :: lndnod

  !============================================================================
  ! 2. Classes et particules

  !     NBPMAX : NOMBRE MAXIMAL DE PARTICULES AUTORISE DANS LE DOMAINE
  !              AU COUR DU CALCUL (UTILE SI INJECTION INSTATIONNAIRE)

  integer, save :: nbpmax

  !=============================================================================
  ! 3. Dimensions des tableaux particulaires

  !     nvp1         : number of real particle variables
  !     ntersl       : nombre de termes sources pour couplage retour
  !     nvlsta       : nombre de variables statistiques
  !     nvisbr       : nombre de variables a enregistrer sur les frontieres

  integer, save ::  nvp1, ntersl, nvlsta, nvisbr

  !=============================================================================

end module lagdim

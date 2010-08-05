!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

! Module for cooling towers

module ctincl

  !=============================================================================

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

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

! Module for specific physics common variables
! between combustion of pulverized coal and heavy fuel

module ppcpfu

  !===========================================================================

  ! XSI         --> XSI = 3,76 pour de l'air

  double precision, save ::  xsi

  ! nb de moles de I dans J

  double precision, save :: ao2f3,acof3,an2f3,ah2of3
  double precision, save :: ao2f4,an2f4,ah2of4,aco2f4
  double precision, save :: ah2of5
  double precision, save :: ao2f6,an2f6,ah2of6,aco2f6
  double precision, save :: ao2f7,an2f7,ah2of7,aco2f7

  ! Equation sur YCO2

  integer, save ::         ieqco2 , iyco2

  ! Combustion heterogene avec le  CO2

  integer, save ::         ihtco2

  ! Equation sur NOX :
  ! ================

  !   IEQNOX = 0 pas de NOx
  !          = 1 calcul du NOx

  integer, save ::         ieqnox

  !   Scalaires supplementaires : fraction massique de HCN et NO
  !                               temperature air

  integer, save ::         iyhcn , iyno , itaire

  !   Propce supplementaires :

  !         Conversion HCN en NO       : EXP(-E1/RT)
  !         Conversion HCN en NO       : EXP(-E2/RT)
  !         NO thermique (Zel'dovitch) : EXP(-E3/RT)

  integer, save ::         ighcn1 , ighcn2 , ignoth

  !   Temperature moyenne d'entree
  !   Taux de vapeur moyen

  double precision, save :: taire

  !=============================================================================

end module ppcpfu

!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

! Module for multigrid parameters

module albase

  !=============================================================================

  !  Methode ale
  !  iale   : utilisation de la methode ALE
  !         = 0 sans methode ALE
  !         = 1 avec methode ALE
  !  iimpal : pointeur sur impale, indicateur de deplacement impose
  !  ixyzn0 : pointeur sur xyzno0, position initiale du maillage
  !  idepal : pointeur sur depale, deplacement du maillage
  !  iialty : pointeur sur ialtyb, type de bord
  !  nalinf : nombre d'iterations d'initialisation du fluide
  !  nalimx : nombre maximal d'iterations d'implicitation du deplacement
  !           des structures
  !  iortvm : type de viscosite de maillage
  !         = 0 isotrope
  !         = 1 orthotrope
  !  epalim : precision relative d'implicitation du deplacement des
  !           structures
  !  italin : iteration d'initialisation de l'ALE
  !         = 0 non
  !         = 1 oui

  integer, save :: iale  , iimpal, ixyzn0, idepal, iialty, nalinf
  integer, save :: nalimx, iortvm, italin

  double precision, save :: epalim

  !=============================================================================

end module albase

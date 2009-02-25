!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

!                              albase.h
!===============================================================================

!  METHODE ALE
!  IALE   : UTILISATION DE LA METHODE ALE
!         = 0 SANS METHODE ALE
!         = 1 AVEC METHODE ALE
!  IIMPAL : POINTEUR SUR IMPALE, INDICATEUR DE DEPLACEMENT IMPOSE
!  IXYZN0 : POINTEUR SUR XYZNO0, POSITION INITIALE DU MAILLAGE
!  IDEPAL : POINTEUR SUR DEPALE, DEPLACEMENT DU MAILLAGE
!  IIALTY : POINTEUR SUR IALTYB, TYPE DE BORD
!  NALINF : NOMBRE D'ITERATIONS D'INITIALISATION DU FLUIDE
!  NALIMX : NOMBRE MAXIMAL D'ITERATIONS D'IMPLICITATION DU DEPLACEMENT
!           DES STRUCTURES
!  IORTVM : TYPE DE VISCOSITE DE MAILLAGE
!         = 0 ISOTROPE
!         = 1 ORTHOTROPE
!  EPALIM : PRECISION RELATIVE D'IMPLICITATION DU DEPLACEMENT DES
!           STRUCTURES
!  ITALIN : ITERATION D'INITIALISATION DE l'ALE
!         = 0 NON
!         = 1 OUI

integer           iale  , iimpal, ixyzn0, idepal, iialty, nalinf
integer           nalimx, iortvm, italin
common / icmale / iale  , iimpal, ixyzn0, idepal, iialty, nalinf, &
                  nalimx, iortvm, italin

double precision epalim
common / rcmale / epalim
! FIN


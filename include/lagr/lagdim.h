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

!                              lagdim.h
!===============================================================================

!===============================================================================

!     Include pour le module Lagrangien : dimensions

!         Trois fichiers complementaires
!                            lagran.h qui porte les non dimensions
!                            lagdim.h qui porte les dimensions variables
!                            lagpar.h qui porte les parametres

!===============================================================================
! 1. Connectivite

!     LONGUEUR DU TABLEAU DU CONNECTIVITE CELLULES -> FACES
!     (calcule dans le sous-programme LAGINI)

integer           lndnod
common / ilagd1 / lndnod

!===============================================================================
! 2. Classes et particules

!     NBPMAX : NOMBRE MAXIMAL DE PARTICULES AUTORISE DANS LE DOMAINE
!              AU COUR DU CALCUL (UTILE SI INJECTION INSTATIONNAIRE)

integer           nbpmax
common / ilagd2 / nbpmax

!===============================================================================
! 3. Dimensions des tableaux particulaires

!     NVP          : NOMBRE DE VARIABLES SUR LES PARTICULES

!     NVP1         : NOMBRE DE VARIABLES SUR LES PARTICULES
!                     EN ENLEVANT POSITION, VITESSE PARTICULE
!                     ET VITESSE FLUIDE

!     NVEP         : NOMBRE D'INFO SUR LES PARTICULES (REELS)

!     NIVEP        : NOMBRE D'INFO SUR LES PARTICULES (ENTIERS)

!     NTERSL       : NOMBRE DE TERMES SOURCES POUR COUPLAGE RETOUR

!     NVLSTA       : NOMBRE DE VARIABLES STATISTIQUES

!     NVISBR       : NOMBRE DE VARIABLES A ENREGISTRER SUR LES FRONTIERES


integer           nvp    , nvp1   , nvep   , nivep  ,             &
                  ntersl , nvlsta , nvisbr
common / ilagd3 / nvp    , nvp1   , nvep   , nivep  ,             &
                  ntersl , nvlsta , nvisbr

! FIN


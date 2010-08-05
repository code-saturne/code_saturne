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

!                              alaste.h
!===============================================================================

!  METHODE ALE - MOUVEMENT DE STRUCTURES EN COUPLAGE AVEC CODE_ASTER

! NTCAST : NUMERO D'ITERATION DE COUPLAGE AVEC CODE_ASTER
! NBASTE : NOMBRE DE STRUCTURES MOBILES
! NBFAST : NOMBRE DE FACES COUPLEES
! NBNAST : NOMBRE DE NOEUDS COUPLES
! ISYNCP : INDICATEUR D'IMPRESSION DES RESULTATS DES DEUX CODES
!          AUX MEMES INSTANTS (SORTIE ENSIGHT POUR ASTER)
! ASDDLF : BLOCAGE DES DDL DE FORCE
! ASDDLC : BLOCAGE DES DDL CINEMATIQUES

integer           ntcast
integer           nbaste, nbfast, nbnast
integer           iforas
integer           isyncp
integer           asddlf(3,nastmx), asddlc(3,nastmx)

common / iaster / ntcast, nbaste, nbfast, nbnast,                 &
                  iforas,                                         &
                  isyncp,                                         &
                  asddlf, asddlc

! FIN


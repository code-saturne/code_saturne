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

!                              vector.h
!===============================================================================

! Longueur de registre LREGIS
!           VR number    VR lenght
!               4       4096
!               8       2048
!              16       1024
!              32        512
!              64        256
!             128        128
!             256         64

! IVECTI,IVECTB             INDICATEUR (1/0) DE VECTORISATION
!                           FORCEE OU NON (FACES INTERN ET DE BRD

integer   lregis
parameter(lregis=1024)
integer           ivecti , ivectb
common / ivecto / ivecti , ivectb

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

!                              parall.h
!===============================================================================


! GESTION DU PARALLELLISME
!   IRANGP : RANG DU PROCESSUS
!     = -1 EN MODE SEQUENTIEL
!     =  R (0 < R < NB_PROCESSUS) EN EXECUTION PARALLELLE
!   NRANGP : NOMBRE DE PROCESSUS (=1 SI SEQUENTIEL)

integer           irangp, nrangp
common / iparal / irangp, nrangp


! DIMENSIONS GLOBALES (I.E. INDEPENDANTES DE DECOUPAGE PARALLELE)
!   NCELGB : NOMBRE DE CELLULES GLOBAL
!   NFACGB : NOMBRE DE FACES INTERNES GLOBAL
!   NFBRGB : NOMBRE DE FACES DE BORD GLOBAL
!   NSOMGB : NOMBRE DE SOMMETS GLOBAL

integer           ncelgb, nfacgb, nfbrgb, nsomgb
common / igeogb / ncelgb, nfacgb, nfbrgb, nsomgb

! FIN

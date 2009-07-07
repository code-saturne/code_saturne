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

!                              parall.h
!===============================================================================


! Gestion du parallellisme

!   irangp : rang du processus
!     = -1 en mode sequentiel
!     =  r (0 < r < nb_processus) en execution parallelle distribuee
!   nrangp : nombre de processus (=1 si sequentiel)
!   nthrdi : nombre de sous ensembles independants de faces internes max
!            dans un groupe (>= 1 avec OpenMP, 1 sinon)
!   nthrdb : nombre de sous ensembles independants de faces de bord max
!            dans un groupe (>= 1 avec OpenMP, 1 sinon)
!   ngrpi  : nombre de groupes de faces internes (> 1 avec OpenMP, 1 sinon)
!   ngrpb  : nombre de groupes de faces de bord (> 1 avec OpenMP, 1 sinon)
!   iompli : bornes par thread pour les faces internes
!   iomplb : bornes par thread pour les faces de bord
!            pour le groupe (couleur) j et le thread i, boucles
!            de iompl.(1, j, i) à iompl.(2, j, i)

integer           irangp, nrangp, nthrdi, nthrdb, ngrpi, ngrpb,           &
                  iompli(2, nthrd1, nthrd2), iomplb(2, nthrd1, nthrd2)
common / iparal / irangp, nrangp, nthrdi, nthrdb, ngrpi, ngrpb,           &
                  iompli, iomplb


! Dimensions globales (i.e. independantes de decoupage parallele)
!   ncelgb : nombre de cellules global
!   nfacgb : nombre de faces internes global
!   nfbrgb : nombre de faces de bord global
!   nsomgb : nombre de sommets global

integer           ncelgb, nfacgb, nfbrgb, nsomgb
common / igeogb / ncelgb, nfacgb, nfbrgb, nsomgb

! FIN

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

! Module for ALE with Code_Aster coupling

module alaste

  !=============================================================================

  ! Nombre de structures max en ALE et Couplage Code_Aster

  integer nastmx
  parameter (nastmx=200)

  !  Methode ALE - mouvement de structures en couplage avec Code_Aster

  ! ntcast : numero d'iteration de couplage avec Code_Aster
  ! nbaste : nombre de structures mobiles
  ! nbfast : nombre de faces couplees
  ! nbnast : nombre de noeuds couples
  ! asddlf : blocage des ddl de force
  ! asddlc : blocage des ddl cinematiques

  integer, save ::  ntcast
  integer, save ::  nbaste, nbfast, nbnast
  integer, save ::  iforas
  integer, save ::  asddlf(3,nastmx), asddlc(3,nastmx)

  !=============================================================================

end module alaste



!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

! Module for multigrid parameters

module mltgrd

  !=============================================================================

  use paramx

  !=============================================================================

  ! Multigrid
  ! -----------
  !   ncegrm : maximum number of cells on coarsest mesh
  !   ngrmax : maximum number of grids

  !   mltmmn : mean number of cells under which merging should take place
  !   mltmgl : global number of cells under which merging should take place
  !   mltmmr : number of active ranks under which no merging is done
  !   mltmst : number of ranks over which merging takes place (stride)
  !
  !   nagmx0 : parametre construction de  maillage automatique
  !   iagmx0 : parametre construction de  maillage automatique
  !   ncpmgr : si > 0, active le post traitement de l'agglomeration, en
  !            projetant les numeros de cellules grossieres sur le
  !            maillage fin (modulo ncpmgr(ivar))

  integer, save :: ncegrm, ngrmax
  integer, save :: mltmmn, mltmgl, mltmst, mltmmr
  integer, save :: nagmx0(nvarmx), iagmx0(nvarmx), ncpmgr(nvarmx)

  !=============================================================================

end module mltgrd

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

  integer, save :: iale  , nalinf
  integer, save :: nalimx, iortvm, italin

  double precision, save :: epalim

  !  impale : indicateur de deplacement impose
  !  xyzno0 : position initiale du maillage
  !  depale : deplacement du maillage
  !  ialtyb : type de bord

  integer, allocatable, dimension(:) :: impale, ialtyb

  double precision, allocatable, dimension(:,:) :: xyzno0, depale

contains

  !=============================================================================

  subroutine init_ale ( ncelet , ncel , nfac , nfabor , nnod )

    use cplsat

    ! Arguments

    integer, intent(in) :: ncelet, ncel, nfac, nfabor, nnod

    if (iale.eq.1.or.imobil.eq.1) then
      allocate(xyzno0(3,nnod))
    endif

    if (iale.eq.1) then
      allocate(impale(nnod))
      allocate(ialtyb(nfabor))
      allocate(depale(nnod,3))
    endif

  end subroutine init_ale

  !=============================================================================

  subroutine finalize_ale

    use cplsat

    if (iale.eq.1.or.imobil.eq.1) then
      deallocate(xyzno0)
    endif

    if (iale.eq.1) then
      deallocate(impale)
      deallocate(depale)
      deallocate(ialtyb)
    endif

  end subroutine finalize_ale

  !=============================================================================

end module albase

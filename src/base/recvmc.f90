!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine recvmc &
!================

 ( rom    , flumas , flumab ,                                     &
   ux     , uy     , uz     ,                                     &
   bx     , by     , bz     )

!===============================================================================
! FONCTION :
! ----------

! RECONSTRUCTION DE LA VITESSE A PARTIR DU FLUX DE MASSE
!     PAR MOINDRES CARRES (VITESSE CONSTANTE PAR ELEMENT)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! rom(ncelet       ! tr ! <-- ! masse volumique aux cellules                   !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! ux   uy          ! tr ! --> ! vitesse reconstruite                           !
! uz   (ncelet)    ! tr !     !                                                !
! bx,y,z(ncelet)   ! tr ! --- ! tableau de travail                             !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use parall
use mesh

!===============================================================================

implicit none

! Arguments

double precision rom(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision ux  (ncelet), uy  (ncelet), uz  (ncelet)
double precision bx(ncelet),   by(ncelet),   bz(ncelet)

! Local variables

integer          ig, it, ii, jj, iel, ifac
integer          irel, idim1, idim2
double precision a11, a22, a33, a12, a13, a23, unsdet
double precision cocg11, cocg12, cocg13, cocg21, cocg22, cocg23
double precision cocg31, cocg32, cocg33
double precision smbx, smby, smbz, unsrho
double precision vecfac, pfacx, pfacy, pfacz

double precision, allocatable, dimension(:,:,:) :: cocg

!===============================================================================

!===============================================================================
! 1. Compute matrix
!===============================================================================

! Allocate a temporary array
allocate(cocg(3,3,ncelet))

! Initialization

!$omp parallel do private(ii, jj)
do iel = 1, ncelet
  do ii = 1, 3
    do jj = 1, 3
      cocg(jj,ii,iel) = 0.d0
    enddo
  enddo
enddo

! Contribution from interior faces

do ig = 1, ngrpi
  !$omp parallel do private(ifac, ii, jj, idim1, idim2, vecfac)
  do it = 1, nthrdi
    do ifac = iompli(1,ig,it), iompli(2,ig,it)
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      do idim1 = 1, 3
        do idim2 = idim1, 3
          vecfac = surfac(idim1,ifac)*surfac(idim2,ifac)
          cocg(idim2,idim1,ii) = cocg(idim2,idim1,ii) + vecfac
          cocg(idim2,idim1,jj) = cocg(idim2,idim1,jj) + vecfac
        enddo
      enddo
    enddo
  enddo
enddo

! Contribution from boundary faces

do ig = 1, ngrpb
  !$omp parallel do private(ifac, ii, idim1, idim2)
  do it = 1, nthrdb
    do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
      ii = ifabor(ifac)
      do idim1 = 1, 3
        do idim2 = idim1, 3
          cocg(idim2,idim1,ii) =   cocg(idim2,idim1,ii) &
                                 + surfbo(idim1,ifac)*surfbo(idim2,ifac)
        enddo
      enddo
    enddo
  enddo
enddo

! Symetrization

!$omp parallel do
do iel = 1, ncel
  cocg(1,2,iel) = cocg(2,1,iel)
  cocg(1,3,iel) = cocg(3,1,iel)
  cocg(2,3,iel) = cocg(3,2,iel)
enddo

!===============================================================================
! 2. Invert matrix
!===============================================================================


!$omp parallel do private(cocg11, cocg12, cocg13, cocg21, cocg22, cocg23, &
!$omp                     cocg31, cocg32, cocg33, &
!$omp                     a11, a12, a13, a22, a23, a33, unsdet)
do iel = 1, ncel

  cocg11 = cocg(1,1,iel)
  cocg12 = cocg(2,1,iel)
  cocg13 = cocg(3,1,iel)
  cocg21 = cocg(1,2,iel)
  cocg22 = cocg(2,2,iel)
  cocg23 = cocg(3,2,iel)
  cocg31 = cocg(1,3,iel)
  cocg32 = cocg(2,3,iel)
  cocg33 = cocg(3,3,iel)

  a11=cocg22*cocg33-cocg32*cocg23
  a12=cocg32*cocg13-cocg12*cocg33
  a13=cocg12*cocg23-cocg22*cocg13
  a22=cocg11*cocg33-cocg31*cocg13
  a23=cocg21*cocg13-cocg11*cocg23
  a33=cocg11*cocg22-cocg21*cocg12

  unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

  cocg(1,1,iel) = a11 *unsdet
  cocg(2,1,iel) = a12 *unsdet
  cocg(3,1,iel) = a13 *unsdet
  cocg(1,2,iel) = a12 *unsdet
  cocg(2,2,iel) = a22 *unsdet
  cocg(3,2,iel) = a23 *unsdet
  cocg(1,3,iel) = a13 *unsdet
  cocg(2,3,iel) = a23 *unsdet
  cocg(3,3,iel) = a33 *unsdet

enddo


!===============================================================================
! 3. Compute RHS
!===============================================================================

!$omp parallel do
do iel = 1, ncelet
  bx(iel) = 0.d0
  by(iel) = 0.d0
  bz(iel) = 0.d0
enddo

! Contribution from interior faces

do ig = 1, ngrpi
  !$omp parallel do private(ifac, ii, jj, pfacx, pfacy, pfacz)
  do it = 1, nthrdi
    do ifac = iompli(1,ig,it), iompli(2,ig,it)
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      pfacx = flumas(ifac)*surfac(1,ifac)
      pfacy = flumas(ifac)*surfac(2,ifac)
      pfacz = flumas(ifac)*surfac(3,ifac)
      bx(ii) = bx(ii) + pfacx
      by(ii) = by(ii) + pfacy
      bz(ii) = bz(ii) + pfacz
      bx(jj) = bx(jj) + pfacx
      by(jj) = by(jj) + pfacy
      bz(jj) = bz(jj) + pfacz
    enddo
  enddo
enddo

! Contribution from boundary faces

do ig = 1, ngrpb
  !$omp parallel do private(ifac, ii)
  do it = 1, nthrdb
    do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
      ii = ifabor(ifac)
      bx(ii) = bx(ii) + flumab(ifac)*surfbo(1,ifac)
      by(ii) = by(ii) + flumab(ifac)*surfbo(2,ifac)
      bz(ii) = bz(ii) + flumab(ifac)*surfbo(3,ifac)
    enddo
  enddo
enddo

!===============================================================================
! 4. Resolution
!===============================================================================

!$omp parallel do private(unsrho, smbx, smby, smbz)
do iel = 1, ncel
  unsrho = 1.d0/rom(iel)
  smbx = bx(iel)
  smby = by(iel)
  smbz = bz(iel)
  ux(iel) = (cocg(1,1,iel)*smbx+cocg(2,1,iel)*smby+cocg(3,1,iel)*smbz)*unsrho
  uy(iel) = (cocg(1,2,iel)*smbx+cocg(2,2,iel)*smby+cocg(3,2,iel)*smbz)*unsrho
  uz(iel) = (cocg(1,3,iel)*smbx+cocg(2,3,iel)*smby+cocg(3,3,iel)*smbz)*unsrho
enddo

! Free memory
deallocate(cocg)

!----
! End
!----

return

end subroutine

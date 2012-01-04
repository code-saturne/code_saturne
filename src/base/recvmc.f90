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

subroutine recvmc &
!================

 ( nvar   , nscal  ,                                              &
   rom    , flumas , flumab ,                                     &
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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! rom(ncelet       ! tr ! <-- ! masse volumique aux cellules                   !
! flumas(nfac)     ! tr ! <-- ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! <-- ! flux de masse aux faces de bord                !
! ux   uy          ! tr ! --> ! vitesse reconstruite                           !
! uz   (ncelet     ! tr !     !                                                !
! bx,y,z(ncelet    ! tr ! --- ! tableau de travail                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use parall
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal


double precision rom(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision ux  (ncelet), uy  (ncelet), uz  (ncelet)
double precision bx(ncelet),   by(ncelet),   bz(ncelet)

! Local variables

integer          lbloc
parameter       (lbloc = 1024)

integer          ii, jj, iel, ifac
integer          ibloc, nbloc, irel, idim1, idim2
double precision aa(lbloc,3,3)
double precision a11, a22, a33, a12, a13, a23, unsdet
double precision cocg11, cocg12, cocg13, cocg21, cocg22, cocg23
double precision cocg31, cocg32, cocg33
double precision smbx, smby, smbz, unsrho
double precision vecfac, pfacx, pfacy, pfacz

double precision, allocatable, dimension(:,:,:) :: cocg

!===============================================================================


!===============================================================================
! 1. CALCUL DE LA MATRICE
!===============================================================================

! Allocate a temporary array
allocate(cocg(ncelet,3,3))

!   INITIALISATION

do ii = 1, 3
  do jj = 1, 3
    do iel = 1, ncelet
      cocg(iel,ii,jj) = 0.d0
    enddo
  enddo
enddo

!   ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

do idim1 = 1, 3
  do idim2 = idim1, 3

    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      vecfac = surfac(idim1,ifac)*surfac(idim2,ifac)
      cocg(ii,idim1,idim2) = cocg(ii,idim1,idim2) + vecfac
      cocg(jj,idim1,idim2) = cocg(jj,idim1,idim2) + vecfac
    enddo

    do ifac = 1, nfabor
      ii = ifabor(ifac)
      cocg(ii,idim1,idim2) = cocg(ii,idim1,idim2)               &
                          + surfbo(idim1,ifac)*surfbo(idim2,ifac)
    enddo

  enddo
enddo


!   SYMETRISATION

do iel = 1, ncel
  cocg(iel,2,1) = cocg(iel,1,2)
  cocg(iel,3,1) = cocg(iel,1,3)
  cocg(iel,3,2) = cocg(iel,2,3)
enddo

!===============================================================================
! 2. INVERSION DE LA MATRICE
!===============================================================================


nbloc = ncel/lbloc
if (nbloc.gt.0) then
  do ibloc = 1, nbloc
    do ii = 1, lbloc
      iel = (ibloc-1)*lbloc+ii

      cocg11 = cocg(iel,1,1)
      cocg12 = cocg(iel,1,2)
      cocg13 = cocg(iel,1,3)
      cocg21 = cocg(iel,2,1)
      cocg22 = cocg(iel,2,2)
      cocg23 = cocg(iel,2,3)
      cocg31 = cocg(iel,3,1)
      cocg32 = cocg(iel,3,2)
      cocg33 = cocg(iel,3,3)

      a11=cocg22*cocg33-cocg32*cocg23
      a12=cocg32*cocg13-cocg12*cocg33
      a13=cocg12*cocg23-cocg22*cocg13
      a22=cocg11*cocg33-cocg31*cocg13
      a23=cocg21*cocg13-cocg11*cocg23
      a33=cocg11*cocg22-cocg21*cocg12

      unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

      aa(ii,1,1) = a11 *unsdet
      aa(ii,1,2) = a12 *unsdet
      aa(ii,1,3) = a13 *unsdet
      aa(ii,2,2) = a22 *unsdet
      aa(ii,2,3) = a23 *unsdet
      aa(ii,3,3) = a33 *unsdet

    enddo

    do ii = 1, lbloc
      iel = (ibloc-1)*lbloc+ii
      cocg(iel,1,1) = aa(ii,1,1)
      cocg(iel,1,2) = aa(ii,1,2)
      cocg(iel,1,3) = aa(ii,1,3)
      cocg(iel,2,2) = aa(ii,2,2)
      cocg(iel,2,3) = aa(ii,2,3)
      cocg(iel,3,3) = aa(ii,3,3)
    enddo

  enddo

endif

irel = mod(ncel,lbloc)
if (irel.gt.0) then
  ibloc = nbloc + 1
  do ii = 1, irel
    iel = (ibloc-1)*lbloc+ii

    cocg11 = cocg(iel,1,1)
    cocg12 = cocg(iel,1,2)
    cocg13 = cocg(iel,1,3)
    cocg21 = cocg(iel,2,1)
    cocg22 = cocg(iel,2,2)
    cocg23 = cocg(iel,2,3)
    cocg31 = cocg(iel,3,1)
    cocg32 = cocg(iel,3,2)
    cocg33 = cocg(iel,3,3)

    a11=cocg22*cocg33-cocg32*cocg23
    a12=cocg32*cocg13-cocg12*cocg33
    a13=cocg12*cocg23-cocg22*cocg13
    a22=cocg11*cocg33-cocg31*cocg13
    a23=cocg21*cocg13-cocg11*cocg23
    a33=cocg11*cocg22-cocg21*cocg12

    unsdet = 1.d0/(cocg11*a11+cocg21*a12+cocg31*a13)

    aa(ii,1,1) = a11 *unsdet
    aa(ii,1,2) = a12 *unsdet
    aa(ii,1,3) = a13 *unsdet
    aa(ii,2,2) = a22 *unsdet
    aa(ii,2,3) = a23 *unsdet
    aa(ii,3,3) = a33 *unsdet

  enddo

  do ii = 1, irel
    iel = (ibloc-1)*lbloc+ii
    cocg(iel,1,1) = aa(ii,1,1)
    cocg(iel,1,2) = aa(ii,1,2)
    cocg(iel,1,3) = aa(ii,1,3)
    cocg(iel,2,2) = aa(ii,2,2)
    cocg(iel,2,3) = aa(ii,2,3)
    cocg(iel,3,3) = aa(ii,3,3)
  enddo
endif


!         MATRICE SYMETRIQUE

do iel = 1, ncel
  cocg(iel,2,1) = cocg(iel,1,2)
  cocg(iel,3,1) = cocg(iel,1,3)
  cocg(iel,3,2) = cocg(iel,2,3)
enddo


!===============================================================================
! 3. CALCUL DU SECOND MEMBRE
!===============================================================================

do iel = 1, ncelet
  bx(iel) = 0.d0
  by(iel) = 0.d0
  bz(iel) = 0.d0
enddo


!     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES

do ifac = 1,nfac
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


!     ASSEMBLAGE A PARTIR DES FACETTES DE BORD

do ifac = 1,nfabor
  ii = ifabor(ifac)
  bx(ii) = bx(ii) + flumab(ifac)*surfbo(1,ifac)
  by(ii) = by(ii) + flumab(ifac)*surfbo(2,ifac)
  bz(ii) = bz(ii) + flumab(ifac)*surfbo(3,ifac)
enddo

!===============================================================================
! 4. RESOLUTION
!===============================================================================


do iel = 1, ncel
  unsrho = 1.d0/rom(iel)
  smbx = bx(iel)
  smby = by(iel)
  smbz = bz(iel)
  ux  (iel) = (cocg(iel,1,1)*smbx+cocg(iel,1,2)*smby              &
              +cocg(iel,1,3)*smbz)*unsrho
  uy  (iel) = (cocg(iel,2,1)*smbx+cocg(iel,2,2)*smby              &
              +cocg(iel,2,3)*smbz)*unsrho
  uz  (iel) = (cocg(iel,3,1)*smbx+cocg(iel,3,2)*smby              &
              +cocg(iel,3,3)*smbz)*unsrho
enddo

! Free memory
deallocate(cocg)

!----
! FIN
!----

return

end subroutine

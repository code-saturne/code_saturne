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

subroutine promav &
!================

 ( ncelet , ncel   , nfac   , isym   , iinvpe ,                   &
   ifacel , da     , xa     , vx     , vy     )

!===============================================================================
! FONCTION :
! ----------

! PRODUIT MATRICE VECTEUR Y = (DA+ XA).X

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! isym             ! e  ! <-- ! indicateur = 1 matrice symetrique              !
!                  !    !     !            = 2 matrice non symetrique          !
! iinvpe           ! e  ! <-- ! indicateur pour annuler les increment          !
!                  !    !     ! en periodicite de rotation (=2) ou             !
!                  !    !     ! pour les echanger normalement de               !
!                  !    !     ! maniere scalaire (=1)                          !
! ifacel(2,nfac    ! te ! <-- ! no des elts voisins d'une face intern          !
! da(ncelet        ! tr ! <-- ! diagonale de la matrice                        !
! xa(nfac,isym)    ! tr ! <-- ! extra diagonale de la matrice                  !
! vx(ncelet        ! tr ! <-- ! vecteur a multiplier                           !
! vy(ncelet        ! tr ! --> ! vecteur resultat                               !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfac   , isym , iinvpe
integer          ifacel(2,nfac)
double precision da(ncelet),xa(nfac,isym),vx(ncelet),vy(ncelet)

! Local variables

integer ifac,ii,jj,iel,idimte,itenso

!===============================================================================

!     1 - PRODUIT MATRICE/VECTEUR PAR LA DIAGONALE
!     --------------------------------------------

!     ESSAYER AUSSI AVEC LES BLAS :
!     CALL DNDOT(NCEL,1,VY,1,1,DA,1,1,VX,1,1)

do iel = 1, ncel
 vy(iel) = da(iel)*vx(iel)
enddo
if(ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    vy(iel) = 0.d0
  enddo
endif


!     2 - PRODUIT MATRICE/VECTEUR TERMES X-TRADIAGONAUX
!     -------------------------------------------------

! ---> TRAITEMENT DU PARALLELISME

if(irangp.ge.0) call parcom (vx)
                !==========

! --> TRAITEMENT DE LA PERIODICITE
if(iperio.eq.1) then
  if(iinvpe.eq.1) then
    idimte = 0
    itenso = 0
    call percom                                                   &
    !==========
  ( idimte , itenso ,                                             &
    vx     , vx     , vx    ,                                     &
    vx     , vx     , vx    ,                                     &
    vx     , vx     , vx    )
   elseif(iinvpe.eq.2) then
      idimte = 0
      itenso = 11
      call percom                                                 &
      !==========
  ( idimte , itenso ,                                             &
    vx     , vx     , vx    ,                                     &
    vx     , vx     , vx    ,                                     &
    vx     , vx     , vx    )

!        Utile a Codits (produit avec une variable non en increment)
    elseif(iinvpe.eq.3) then
      idimte = 0
      itenso = 1
      call percom                                                 &
      !==========
  ( idimte , itenso ,                                             &
    vx     , vx     , vx    ,                                     &
    vx     , vx     , vx    ,                                     &
    vx     , vx     , vx    )
    endif
endif


if( isym.eq.1 ) then

  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    vy(ii) = vy(ii) +xa(ifac,1)*vx(jj)
    vy(jj) = vy(jj) +xa(ifac,1)*vx(ii)
  enddo

else

  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    vy(ii) = vy(ii) +xa(ifac,1)*vx(jj)
    vy(jj) = vy(jj) +xa(ifac,2)*vx(ii)
  enddo

endif

!--------
! FORMATS
!--------


!----
! FIN
!----

return

end subroutine

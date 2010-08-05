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

subroutine vectds &
!================

 ( ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   ifacel , ifabor , ia     ,                                     &
   surfac , surfbo , pond   ,                                     &
   vectx  , vecty  , vectz  ,                                     &
   valf   , valb   , ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL A LA FACE DE (VECT)ij . S
!  A PARTIR DU VECTEUR VECTX, VECTY, VECTZ
! UTILISE POUR LE CALCUL DU TERME DE DIFFUSION DE Rij ET Epsilon
!  EN Rij-Epsilon LRR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! ifacel(2,nfac    ! te ! <-- ! no des elts voisins d'une face intern          !
! ifabor(nfabor    ! te ! <-- ! no de l'elt voisin d'une face de bord          !
! ia(*)            ! ia ! --- ! main integer work array                        !
! surfac(3,nfac    ! tr ! <-- ! surf vectorielle des surfaces interne          !
! surfbo(          ! tr ! <-- ! surf vectorielle des surfaces                  !
!  (3,nfabor       !    !     !                               de bord          !
! pond(nfac        ! tr ! <-- ! ponderation pour interpolation faces           !
! vectx (ncelet    ! tr ! <-- ! composante x du vecteur   entre                !
! vecty (ncelet    ! tr ! <-- ! composante y du vecteur   entre                !
! vectz (ncelet    ! tr ! <-- ! composante z du vecteur   entre                !
! valf (nfac)      ! tr ! --> ! vect*surface      aux faces internes           !
! valb (nfabor     ! tr ! --> ! vect*surface      aux faces de bord            !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks

include "paramx.f90"
include "pointe.f90"
include "period.f90"
include "parall.f90"

!===============================================================================

! Arguments

integer          ndim, ncelet, ncel, nfac, nfabor
integer          ifacel(2,nfac), ifabor(nfabor)
integer          ia(*)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision pond(nfac)
double precision vectx(ncelet), vecty(ncelet), vectz(ncelet)
double precision valf(nfac), valb(nfabor)
double precision ra(*)

! Local variables

integer          ifac, iel1, iel2
double precision valfx, valfy, valfz

!===============================================================================

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(vectx, vecty, vectz)
  !==========
endif


do ifac = 1 , nfac

  iel1 = ifacel(1,ifac)
  iel2 = ifacel(2,ifac)

  valfx =       pond(ifac)  * vectx(iel1) +                       &
          (1.d0-pond(ifac)) * vectx(iel2)
  valfy =       pond(ifac)  * vecty(iel1) +                       &
          (1.d0-pond(ifac)) * vecty(iel2)
  valfz =       pond(ifac)  * vectz(iel1) +                       &
          (1.d0-pond(ifac)) * vectz(iel2)

  valf(ifac) = valfx*surfac(1,ifac) +                             &
               valfy*surfac(2,ifac) +                             &
               valfz*surfac(3,ifac)
 enddo

 do ifac = 1 , nfabor

!     On met VALB a zero, ce qui revient a negliger la partie
!       extradiagonale du tenseur de diffusion au bord.
!MO          IEL1 = IFABOR(IFAC)
!MOC
!MO          VALB(IFAC) = VECTX(IEL1)*SURFBO(1,IFAC) +
!MO     &                 VECTY(IEL1)*SURFBO(2,IFAC) +
!MO     &                 VECTZ(IEL1)*SURFBO(3,IFAC)
    valb(ifac) = 0.d0

 enddo

 return
 end

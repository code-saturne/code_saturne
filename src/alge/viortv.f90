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

subroutine viortv &
!================

 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DE LA VITESSE DE DIFFUSION "ORTHOTROPE"
! VISCF,B = VISCOSITE*SURFACE/DISTANCE, HOMOGENE A UN DEBIT EN KG/S

!         =
! (NX**2*VISC11_MOY_FACE
! +NY**2*VISC22_MOY_FACE+NZ**2*VISC33_MOY_FACE)*SURFACE/DISTANCE

! LA VISCOSITE EST DONNE PAR W1, W2, W3

! RQE : A PRIORI, PAS BESOIN DE TECHNIQUE DE RECONSTRUCTION
!  ( A AMELIORER SI NECESSAIRE )

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! imvisf           ! e  ! <-- ! methode de calcul de la visc face              !
!                  !    !     !  = 0 arithmetique                              !
!                  !    !     !  = 1 harmonique                                !
! w1(3,ncelet)     ! tr ! <-- ! valeurs de la viscosite                        !
! viscf(nfac)      ! tr ! --> ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --> ! visc*surface/dist aux faces de bord            !
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
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          imvisf


double precision w1(3,ncelet)
double precision viscf(nfac), viscb(nfabor)

! Local variables

integer          ifac, ii, jj
double precision viscxi, viscxj, viscyi, viscyj, visczi, visczj
double precision sx2, sy2, sz2, distbf, pnd, surfn

!===============================================================================

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call syndin(w1)
  !==========
endif


if( imvisf.eq.0 ) then

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    surfn = surfan(ifac)

    viscxi = w1(1,ii)
    viscxj = w1(1,jj)
    viscyi = w1(2,ii)
    viscyj = w1(2,jj)
    visczi = w1(3,ii)
    visczj = w1(3,jj)

    sx2    = surfac(1,ifac)**2
    sy2    = surfac(2,ifac)**2
    sz2    = surfac(3,ifac)**2

    viscf(ifac) = 0.5d0*(                                         &
       (viscxi+viscxj)*sx2                                        &
     + (viscyi+viscyj)*sy2                                        &
     + (visczi+visczj)*sz2 ) / (surfn*dist(ifac))

  enddo

else

  do ifac = 1,nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    surfn = surfan(ifac)
    pnd  = pond(ifac)

    viscxi = w1(1,ii)
    viscxj = w1(1,jj)
    viscyi = w1(2,ii)
    viscyj = w1(2,jj)
    visczi = w1(3,ii)
    visczj = w1(3,jj)

    sx2    = surfac(1,ifac)**2
    sy2    = surfac(2,ifac)**2
    sz2    = surfac(3,ifac)**2

    viscf(ifac) =                                                 &
      ( viscxi*viscxj*sx2                                         &
              /(pnd*viscxi+(1.d0-pnd)*viscxj)                     &
      + viscyi*viscyj*sy2                                         &
              /(pnd*viscyi+(1.d0-pnd)*viscyj)                     &
      + visczi*visczj*sz2                                         &
              /(pnd*visczi+(1.d0-pnd)*visczj)                     &
       ) /(surfn*dist(ifac))
  enddo

endif

do ifac=1,nfabor

  ii = ifabor(ifac)

  surfn = surfbn(ifac)
  distbf = distb(ifac)

  viscxi = w1(1,ii)
  viscyi = w1(2,ii)
  visczi = w1(3,ii)

  sx2    = surfbo(1,ifac)**2
  sy2    = surfbo(2,ifac)**2
  sz2    = surfbo(3,ifac)**2

  viscb(ifac) =                                                   &
    (viscxi*sx2+viscyi*sy2+visczi*sz2)/(surfn*distbf)

enddo

!----
! FIN
!----

return

end subroutine

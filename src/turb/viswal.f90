!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine viswal

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE "TURBULENTE" POUR
!          UN MODELE LES WALE

! VISCT = ROM * (CWALE*L)**2 *
!   [(Sijd.Sijd)**(3/2)] / [(Sij.Sij)**(5/2) + (Sijd.Sijd)**(5/4)]
!
! avec
!       Sij = 0.5*[DUi/Dxj + DUj/Dxi]
! et
!       Sijd = 0.5*[DUi/Dxk.DUk/Dxj + DUj/Dxk.DUk/Dxi] - 1/3*Delta_ij.Gkk**2

! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use numvar
use optcal
use dimens, only: nvar
use cstphy
use entsor
use parall
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iel, inc
integer          ipcvst, ipcvis, iprev
integer          i, j, k

double precision coef, deux, delta, tiers
double precision sij, sijd, s, sd, sinv
double precision xfil, xa  , xb  , radeux, con
double precision dudx(ndim,ndim), kdelta(ndim,ndim)

double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_val_s(iprpfl(iviscl), viscl)
call field_get_val_s(iprpfl(ivisct), visct)
call field_get_val_s(icrom, crom)

! --- Pour le calcul de la viscosite de sous-maille
xfil   = xlesfl
xa     = ales
xb     = bles
deux   = 2.d0
radeux = sqrt(deux)
tiers  = 1.d0/3.d0


!===============================================================================
! 2.  CALCUL DU GRADIENT DE VITESSE
!       W1 = DU/DX, W2 = DU/DY, W3 = DU/DZ
!       W4 = DV/DX, W5 = DV/DY, W6 = DV/DZ
!       W7 = DW/DX, W8 = DW/DY, W9 = DW/DZ
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

! Kronecker delta Dij

kdelta(1,1) = 1
kdelta(1,2) = 0
kdelta(1,3) = 0
kdelta(2,1) = 0
kdelta(2,2) = 1
kdelta(2,3) = 0
kdelta(3,1) = 0
kdelta(3,2) = 0
kdelta(3,3) = 1

coef = cwale**2 * radeux

do iel = 1, ncel

  ! Dudx is interleaved, but not gradv...
  ! gradv(iel, xyz, uvw)
  dudx(1,1) = gradv(1,1,iel)
  dudx(1,2) = gradv(2,1,iel)
  dudx(1,3) = gradv(3,1,iel)
  dudx(2,1) = gradv(1,2,iel)
  dudx(2,2) = gradv(2,2,iel)
  dudx(2,3) = gradv(3,2,iel)
  dudx(3,1) = gradv(1,3,iel)
  dudx(3,2) = gradv(2,3,iel)
  dudx(3,3) = gradv(3,3,iel)

  s  = 0.d0
  sd = 0.d0

  do i = 1, ndim
    do j = 1, ndim

      ! Sij = 0.5 * (dUi/dXj + dUj/dXi)

      sij = 0.5d0*(dudx(i,j)+dudx(j,i))

      s = s + sij**2

      do k = 1, ndim

!  traceless symmetric part of the square of the velocity gradient tensor
!    Sijd = 0.5 * ( dUi/dXk dUk/dXj + dUj/dXk dUk/dXi) - 1/3 Dij dUk/dXk dUk/dXk

        sijd = 0.5d0*(dudx(i,k)*dudx(k,j)+ dudx(j,k)*dudx(k,i)) &
              -tiers*kdelta(i,j)*dudx(k,k)**2

        sd = sd + sijd**2

      enddo
    enddo
  enddo

!===============================================================================
! 3.  CALCUL DE LA VISCOSITE TURBULENTE
!===============================================================================

  ! Turbulent inverse time scale =
  !   (Sijd Sijd)^3/2 / [ (Sij Sij)^5/2 + (Sijd Sijd)^5/4 ]

  sinv = (s**2.5d0 + sd**1.25d0)
  if (sinv.gt.0.d0) then
    con = sd**1.5d0 / sinv
  else
    con = 0.d0
  endif

  delta = xfil* (xa*volume(iel))**xb
  delta = coef * delta**2

  visct(iel) = crom(iel) * delta * con

enddo

! Free memory
deallocate(gradv)

!----
! FORMAT
!----

!----
! FIN
!----

return

end subroutine

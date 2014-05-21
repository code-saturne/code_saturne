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

subroutine vislmg

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE TURBULENTE POUR
!          UN MODELE DE LONGUEUR DE MELANGE SIMPLE

! VISCT = ROM * (XKAPPA * L) **2 * SQRT ( 2 * Sij.Sij )
!       Sij = (DUi/Dxj + DUj/Dxi)/2

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
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iel, inc
integer          ipcvst, iprev

double precision coef, deux

double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: visct

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

call field_get_val_s(iprpfl(ivisct), visct)
call field_get_val_s(icrom, crom)

!===============================================================================
! 2.  CALCUL DES GRADIENTS DE VITESSE ET DE
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

do iel = 1, ncel
  visct(iel) = &
      gradv(1, 1, iel)**2 + gradv(2, 2, iel)**2 + gradv(3, 3, iel)**2  &
    + 0.5d0*( (gradv(2, 1, iel) + gradv(1, 2, iel))**2                 &
            + (gradv(3, 1, iel) + gradv(1, 3, iel))**2                 &
            + (gradv(3, 2, iel) + gradv(2, 3, iel))**2 )
enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  CALCUL DE LA VISCOSITE (DYNAMIQUE)
!===============================================================================

deux = 2.d0
coef = (xkappa*xlomlg)**2 * sqrt(deux)

do iel = 1, ncel
  visct(iel) = crom(iel) * coef * sqrt(visct(iel))
enddo

!----
! FORMAT
!----


!----
! FIN
!----

return
end subroutine

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine laggra &
!================

 ( iprev, gradpr , gradvf )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    1)   Calcul de (- (GRADIENT DE PRESSION) / ROM )

!    2)   Calcul du gradient de Vitesse

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iprev            ! i  ! --> ! use value from previous time step              !
! gradpr(3,ncelet) ! ra ! --> ! gradient de pression                           !
! gradvf           ! ra ! --> ! gradient de vitesse fluide                     !
!   (3,3,ncelet)   !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use pointe
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

integer          iprev
double precision gradpr(3,ncelet) , gradvf(3,3,ncelet)

! Local variables

integer          inc , iccocg
integer          iel
double precision unsrho

double precision, dimension(:), pointer :: cromf

!===============================================================================

!===============================================================================
! 0. Parameters for gradient computation
!===============================================================================

inc     = 1
iccocg  = 1

!===============================================================================
! 1. Compute pressure gradient / rho
!===============================================================================
! FIXME for iphydr = 1 and 2
call field_gradient_scalar(ivarfl(ipr), iprev, imrgra, inc,          &
                           iccocg,                                   &
                           gradpr)

! Pointeur sur la masse volumique en fonction de l'ecoulement

if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then
  call field_get_val_s(iprpfl(ipproc(irom1)), cromf)
else
  call field_get_val_s(icrom, cromf)
endif

! Calcul de -Grad P / Rom

! Warning, in standard calculation, the computed pressure is
! the hydrostatic pressure and not the real one
if (iphydr.eq.0.and.ippmod(icompf).lt.0) then
  do iel = 1, ncel
    gradpr(1,iel) = gradpr(1,iel) + ro0*gx
    gradpr(2,iel) = gradpr(2,iel) + ro0*gy
    gradpr(3,iel) = gradpr(3,iel) + ro0*gz
  enddo
endif

do iel = 1, ncel
  unsrho = 1.d0 / cromf(iel)
  gradpr(1,iel) = -gradpr(1,iel) * unsrho
  gradpr(2,iel) = -gradpr(2,iel) * unsrho
  gradpr(3,iel) = -gradpr(3,iel) * unsrho
enddo

!===============================================================================
! 2. Compute velocity gradient
!===============================================================================

if (modcpl.gt.0 .and. iplas.ge.modcpl) then

  call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,          &
                             gradvf)

endif

!----
! End
!----

return

end subroutine

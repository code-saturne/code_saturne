!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine tspdcv &
 ( ncepdp , icepdc ,                                              &
   vela   ,                                                       &
   ckupdc , trav   )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES TERMES DE PERTE DE CHARGE
!        POUR LE BILAN EXPLICITE
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! trav(ncelet,3    ! tr ! <-- ! tableau des second membres                     !
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
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          ncepdp
integer          icepdc(ncepdp)

double precision ckupdc(6,ncepdp)
double precision trav(3,ncelet)
double precision vela  (3  ,ncelet)

! Local variables

integer          iel   , ielpdc
double precision romvom, vit1  , vit2  , vit3
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision, dimension(:), pointer :: crom, croma

!===============================================================================


call field_get_val_s(icrom, crom)

if (iroext.gt.0.and.isno2t.gt.0) then
  call field_get_val_prev_s(icrom, croma)
endif

!     La matrice est toujours "implicite"

do ielpdc = 1, ncepdp

  iel    = icepdc(ielpdc)
  romvom =-crom(iel)*cell_f_vol(iel)
  cpdc11 = ckupdc(1,ielpdc)
  cpdc22 = ckupdc(2,ielpdc)
  cpdc33 = ckupdc(3,ielpdc)
  cpdc12 = ckupdc(4,ielpdc)
  cpdc23 = ckupdc(5,ielpdc)
  cpdc13 = ckupdc(6,ielpdc)
  vit1   = vela(1,iel)
  vit2   = vela(2,iel)
  vit3   = vela(3,iel)

  trav(1,iel) = trav(1,iel) + romvom*(cpdc11*vit1 + cpdc12*vit2 + cpdc13*vit3)
  trav(2,iel) = trav(2,iel) + romvom*(cpdc12*vit1 + cpdc22*vit2 + cpdc23*vit3)
  trav(3,iel) = trav(3,iel) + romvom*(cpdc13*vit1 + cpdc23*vit2 + cpdc33*vit3)

enddo

return
end subroutine

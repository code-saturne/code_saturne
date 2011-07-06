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

subroutine tspdcv &
!================

 ( nvar   , nscal  , ncepdp ,                                     &
   icepdc ,                                                       &
   rtpa   , vela   ,                                              &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  , ckupdc , trav   )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp
integer          icepdc(ncepdp)

double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision trav(3,ncelet)
double precision vela  (3  ,ncelet)

! Local variables

integer          iel   , ielpdc
integer          ipcrom, ipcroo
double precision romvom, vit1  , vit2  , vit3
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23

!===============================================================================


ipcrom = ipproc(irom)

ipcroo = ipcrom
if(iroext.gt.0.and.isno2t.gt.0) then
  ipcroo = ipproc(iroma )
endif

!     La matrice est toujours "implicite"

do ielpdc = 1, ncepdp

  iel    = icepdc(ielpdc)
  romvom =-propce(iel,ipcrom)*volume(iel)
  cpdc11 = ckupdc(ielpdc,1)
  cpdc22 = ckupdc(ielpdc,2)
  cpdc33 = ckupdc(ielpdc,3)
  cpdc12 = ckupdc(ielpdc,4)
  cpdc13 = ckupdc(ielpdc,5)
  cpdc23 = ckupdc(ielpdc,6)
  vit1   = vela(1,iel)
  vit2   = vela(2,iel)
  vit3   = vela(3,iel)

  trav(1,iel) = trav(1,iel) + romvom*(cpdc11*vit1 + cpdc12*vit2 + cpdc13*vit3)
  trav(2,iel) = trav(2,iel) + romvom*(cpdc12*vit1 + cpdc22*vit2 + cpdc23*vit3)
  trav(3,iel) = trav(3,iel) + romvom*(cpdc13*vit1 + cpdc23*vit2 + cpdc33*vit3)

enddo

return
end subroutine

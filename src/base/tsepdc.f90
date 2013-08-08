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

subroutine tsepdc &
!================

 ( ncepdp , idiaex , icepdc ,                                     &
   rtpa   , propce ,                                              &
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
! idiaex           ! e  ! <-- ! indicateur de traitement de la                 !
!                  !    !     ! diagonale (=1) ou extradiagonale (=2)          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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

integer          ncepdp
integer          idiaex

integer          icepdc(ncepdp)

double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6)
double precision trav(ncelet,3)

! Local variables

integer          iel   , ielpdc
integer          ipcrom, ipcroo
double precision romvom, vit1  , vit2  , vit3
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23

!===============================================================================


ipcrom = ipproc(irom  )

ipcroo = ipcrom
if(iroext.gt.0.and.isno2t.gt.0) then
  ipcroo = ipproc(iroma )
endif

!     La diagonale est toujours "implicite"

if(idiaex.eq.1) then

  do ielpdc = 1, ncepdp

    iel    = icepdc(ielpdc)
    romvom =-propce(iel,ipcrom)*volume(iel)
    cpdc11 = ckupdc(ielpdc,1)
    cpdc22 = ckupdc(ielpdc,2)
    cpdc33 = ckupdc(ielpdc,3)
    vit1   = rtpa(iel,iu)
    vit2   = rtpa(iel,iv)
    vit3   = rtpa(iel,iw)

    trav(iel,1) = trav(iel,1) +                                   &
         romvom * ( cpdc11*vit1                             )
    trav(iel,2) = trav(iel,2) +                                   &
         romvom * (               cpdc22*vit2               )
    trav(iel,3) = trav(iel,3) +                                   &
         romvom * (                             cpdc33*vit3 )

  enddo

endif

!     L'extradiagonale est explicite mais peut etre extrapolee

if(idiaex.eq.2) then

  do ielpdc = 1, ncepdp

    iel    = icepdc(ielpdc)
    romvom =-propce(iel,ipcroo)*volume(iel)
    cpdc12 = ckupdc(ielpdc,4)
    cpdc13 = ckupdc(ielpdc,5)
    cpdc23 = ckupdc(ielpdc,6)
    vit1   = rtpa(iel,iu)
    vit2   = rtpa(iel,iv)
    vit3   = rtpa(iel,iw)

    trav(iel,1) = trav(iel,1) +                                   &
         romvom * (               cpdc12*vit2 + cpdc13*vit3 )
    trav(iel,2) = trav(iel,2) +                                   &
         romvom * ( cpdc12*vit1               + cpdc23*vit3 )
    trav(iel,3) = trav(iel,3) +                                   &
         romvom * ( cpdc13*vit1 + cpdc23*vit2               )

  enddo

endif


return

end subroutine

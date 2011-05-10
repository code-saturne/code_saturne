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

subroutine tsepdc &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp ,                            &
   idiaex ,                                                       &
   icepdc ,                                                       &
   ia     ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! idiaex           ! e  ! <-- ! indicateur de traitement de la                 !
!                  !    !     ! diagonale (=1) ou extradiagonale (=2)          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! trav(ncelet,3    ! tr ! <-- ! tableau des second membres                     !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp
integer          idiaex

integer          icepdc(ncepdp)
integer          ia(*)

double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision trav(ncelet,3)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel   , ielpdc
integer          iuiph , iviph , iwiph , ipcrom, ipcroo
double precision romvom, vit1  , vit2  , vit3
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23

!===============================================================================

idebia = idbia0
idebra = idbra0

iuiph  = iu
iviph  = iv
iwiph  = iw
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
    vit1   = rtpa(iel,iuiph)
    vit2   = rtpa(iel,iviph)
    vit3   = rtpa(iel,iwiph)

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
    vit1   = rtpa(iel,iuiph)
    vit2   = rtpa(iel,iviph)
    vit3   = rtpa(iel,iwiph)

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

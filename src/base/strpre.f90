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

subroutine strpre &
!================

 ( itrale , italim , ineefl ,                                     &
   impale ,                                                       &
   flmalf , flmalb , xprale , cofale )

!===============================================================================
! FONCTION :
! ----------

! PREDICTION DU DEPLACEMENT DES STRUCTURES EN ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itrale           ! e  ! <-- ! numero d'iteration pour l'ale                  !
! italim           ! e  ! <-- ! numero d'iteration couplage implicite          !
! ineefl           ! e  ! <-- ! indicateur de sauvegarde des flux              !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! flmalf(nfac)     ! tr ! --> ! sauvegarde du flux de masse faces int          !
! flmalb(nfabor    ! tr ! --> ! sauvegarde du flux de masse faces brd          !
! cofale           ! tr ! --> ! sauvegarde des cl de p et u                    !
!    (nfabor,8)    !    !     !                                                !
! xprale(ncelet    ! tr ! --> ! sauvegarde de la pression, si nterup           !
!                  !    !     !    est >1                                      !
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
use optcal
use numvar
use pointe
use albase, only: nalimx
use alstru
use alaste
use parall
use period
use entsor
use albase, only: fdiale
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          itrale , italim , ineefl

integer          impale(nnod)

double precision flmalf(nfac), flmalb(nfabor), xprale(ncelet)
double precision cofale(nfabor,11)

! Local variables

integer          istr, ii, ifac, inod, iel
integer          iflmas, iflmab

double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:), pointer :: disale
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: coefap, coefbp

double precision, dimension(:), pointer :: cvara_pr

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

call field_get_val_v(fdiale, disale)

call field_get_val_prev_s(ivarfl(ipr), cvara_pr)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

call field_get_coefa_s(ivarfl(ipr), coefap)
call field_get_coefb_s(ivarfl(ipr), coefbp)

!===============================================================================
! 2. PREDICTION DU DEPLACEMENT DES STRUCTURES
!===============================================================================

! 2.1 STRUCTURES INTERNES :
! -----------------------

! Lors de la phase d'initialisation de l'ALE (ITRALE=0), XSTP contient
!    - la valeur du deplacement initial des structures si l'utilisateur
!        les a touchees (suite ou pas)
!    - 0 si debut de calcul avec les structures
!    - le deplacement utilise pour le calcul precedent si suite sans
!        modification utilisateur
!   Il faut cependant transferer sa valeur dans XSTR (qui est utilise
!   par Newmark)
! Lors des iterations suivantes (ITRALE>0) on utilise les schemas standard
!   de calcul de XSTP

if (nbstru.gt.0) then

  if (itrale.eq.0) then

    do istr = 1, nbstru
      do ii = 1, ndim
        xstr(ii,istr) = xstp(ii,istr)
      enddo
    enddo

  else

! 2.1.1 : SCHEMA DE COUPLAGE EXPLICITE
!---------------------------------------------
    if (nalimx.eq.1) then
      do istr = 1, nbstru
        do ii = 1, 3
          xstp(ii,istr) = xstr(ii,istr)                           &
             + aexxst*dtstr(istr)*xpstr(ii,istr)                  &
             + bexxst*dtstr(istr)*(xpstr(ii,istr)-xpsta(ii,istr))
        enddo
      enddo

! 2.1.2 : SCHEMA DE COUPLAGE IMPLICITE
!---------------------------------------------
    else
      do istr = 1, nbstru
        do ii = 1, 3
          xstp(ii,istr)  = xstr(ii,istr)
        enddo
      enddo
    endif

  endif

  do ifac = 1, nfabor
    istr = idfstr(ifac)
    if (istr.gt.0) then
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        impale(inod) = 1
        disale(1,inod) = xstp(1,istr)
        disale(2,inod) = xstp(2,istr)
        disale(3,inod) = xstp(3,istr)
      enddo
    endif
  enddo

endif

! 2.2 STRUCTURES EXTERNES (COUPLAGE CODE_ASTER) :
! -----------------------

if (nbaste.gt.0) then

  do ifac = 1, nfabor
    istr = idfstr(ifac)
    if (istr.lt.0) then
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        impale(inod) = 1
      enddo
    endif
  enddo

! Si ITRALE = 0, on ne fait rien pour l'instant, mais il faudrait
! prevoir la reception de deplacement initiaux venant de Code_Aster
  if (itrale.gt.0) then

    ntcast = ntcast + 1

    ! Reception des deplacements predits et remplissage de disale

    call astcin(ntcast, disale)
    !==========

  endif

endif

!===============================================================================
! 3. DEPLACEMENT AU PAS DE TEMPS PRECEDENT ET SAUVEGARDE FLUX ET DE P
!===============================================================================

if (italim.eq.1) then
  do istr = 1, nbstru
    do ii = 1, 3
      xsta(ii,istr)   = xstr(ii,istr)
      xpsta(ii,istr)  = xpstr(ii,istr)
      xppsta(ii,istr) = xppstr(ii,istr)
    enddo
  enddo
  if (ineefl.eq.1) then
    do ifac = 1, nfac
      flmalf(ifac) = imasfl(ifac)
    enddo
    do ifac = 1, nfabor
      flmalb(ifac) = bmasfl(ifac)
      cofale(ifac,1)  = coefap(ifac)
      cofale(ifac,2)  = coefau(1, ifac)
      cofale(ifac,3)  = coefau(2, ifac)
      cofale(ifac,4)  = coefau(3, ifac)
      cofale(ifac,5)  = coefbp(ifac)
      ! the coefficient B is supposed to be symmetric
      cofale(ifac,6)  = coefbu(1, 1, ifac)
      cofale(ifac,7)  = coefbu(2, 2, ifac)
      cofale(ifac,8)  = coefbu(3, 3, ifac)
      cofale(ifac,9)  = coefbu(1, 2, ifac)
      cofale(ifac,10) = coefbu(2, 3, ifac)
      cofale(ifac,11) = coefbu(1, 3, ifac)
    enddo
    if (nterup.gt.1) then
      do iel = 1, ncelet
        xprale(iel) = cvara_pr(iel)
      enddo
    endif
  endif
endif

!--------
! Formats
!--------

!----
! End
!----

end subroutine

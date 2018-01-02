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

subroutine vor2cl &
 ( itypfb ,                                                       &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    TRANSFERT DES VORTEX DANS LES TABLEAUX RCDOCL
!    AVEC CHANGEMENT DE REPERE EVENTUEL
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itypfb           ! ia ! --> ! boundary face types                            !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!                  !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/turb_schmidt)*gradt    !
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
use cstphy
use cstnum
use dimens, only: nvar
use entsor
use parall
use period
use vorinc
use mesh

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iel, ii, ient
double precision xu, xv, xw

integer          ipass
data             ipass /0/
save             ipass
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


ipass = ipass + 1
if(irangp.ge.0.and.ipass.eq.1) then
  do ii = 1, nnent
    call parbcr(0,3,dir1(1,ii))
    call parbcr(0,3,dir2(1,ii))
    call parbcr(0,3,dir3(1,ii))
  enddo
endif

! on envoie la vitesse calcule par le processeur 0
! a tous les autres processeurs

if(irangp.ge.0) then
  do ient = 1, nnent
    call parbcr(0,icvmax,uvort(1,ient))
    call parbcr(0,icvmax,vvort(1,ient))
    call parbcr(0,icvmax,wvort(1,ient))
  enddo
endif

do ii = 1, nnent
  icvor2(ii) = 0
enddo

do ifac = 1, nfabor

  iel = ifabor(ifac)
  ient = irepvo(ifac)

  if(ient.ne.0) then

    icvor2(ient) = icvor2(ient) + 1

    itypfb(ifac) = ientre

    ii = ifacgl(icvor2(ient),ient)

    xu = uvort(ii,ient)
    xv = vvort(ii,ient)
    xw = wvort(ii,ient)

    rcodcl(ifac,iu,1) = xu*dir3(1,ient) + xv*dir1(1,ient) + xw*dir2(1,ient)
    rcodcl(ifac,iv,1) = xu*dir3(2,ient) + xv*dir1(2,ient) + xw*dir2(2,ient)
    rcodcl(ifac,iw,1) = xu*dir3(3,ient) + xv*dir1(3,ient) + xw*dir2(3,ient)

  endif

enddo

! ---
! FIN
! ---

return
end subroutine

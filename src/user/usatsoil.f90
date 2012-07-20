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

subroutine usatsoil &
     !==================
     ( iappel )

!===============================================================================
! Purpose:
! -------

!     User subroutine.

!     Data Entry for the atmospheric ground model .


! Introduction:
!=============

! Define the different values which can be taken by iappel:
!--------------------------------------------------------

! iappel = 1 (only one call on initialization):
!            Computation of the cells number where we impose a
!            Ground Model

! iappel = 2 (only one call on initialization):
!            users may defined the ground face composition
!            Warning : be coherent with the dimension of the array pourcent_sol
!            It's also possible to modified the tab_sol array of the ground
!            type constants
!
!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use mesh

!===============================================================================

implicit none

! Arguments
!-------------------------------------------------------------------
integer          iappel

! Local variables
!-------------------------------------------------------------------
integer          ifac , ifbt1d , ilelt , nlelt , isol

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

ifbt1d = 0
allocate(lstelt(nfabor))

!===============================================================================
! APPEL 1.  INITIALISATIONS
!===============================================================================
if (iappel.eq.1) then
!On precise la couleur du sol
    CALL GETFBR('75',NLELT,LSTELT)
    do ilelt = 1, nlelt
      ifbt1d =ifbt1d + 1
    enddo
    nfmodsol = ifbt1d

    allocate(indsol(nfmodsol))

    do ilelt = 1, nlelt
      ifac = lstelt(ilelt)
      indsol(ilelt) = ifac
    enddo
!On precise le nombre sol utilise pour le modele
! 5 dans le cas bati, 7 dans le cas bati dense/mixte/diffus
    nbrsol   = 5
!On renseigne la teneur en eau des deux reservoirs
!(necessaire pour l'initialisation) arguments rajoutes dans
! atincl.h
    w1ini = 0.d0
    w2ini = 0.0d0
endif


if (iappel.eq.2) then

!Modification pour cas Wangara, dans ce cas la on a Csol(mineral=4) = 1.7e-5
! ainsi que zoth = 1.2e-3
    tab_sol(4)%csol = 1.7e-5
    tab_sol(4)%rugthe = 0.0012

!
!Initialisation of the pourcent_sol array
    do ifac = 1,nfmodsol
        do isol = 1,nbrsol
            pourcent_sol(ifac,isol)=0
        enddo
        pourcent_sol(ifac,4) = 100
    enddo
endif

!===============================================================================

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine usatsoil

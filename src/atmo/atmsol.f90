!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!> \file atmsol.f90
!> \brief    build constants and variables to describe the ground model
!
!> \brief    build constants and variables to describe ground model
!>-     NB : soil model structures defined in module atsoil.f90

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     iappel        Calling number: 1 allocation, 2 fillign arrays
!______________________________________________________________________________

subroutine atmsol &
     ( iappel )

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

integer          iappel

! Local variables

integer          error

!===============================================================================

if (iatsoil.ge.0) then
  ! First call, define nfmodsol
  if (iappel.eq.1) then

    call usatsoil(iappel)

    ! Allocation of table of values, TODO move to zone
    allocate(tab_sol(nbrsol),stat = error)
    call solcat( error )

    if (error /= 0) then
      write(nfecra,*) "Allocation error of atmodsol::tab_sol"
      call csexit(1)
    endif

    ! We continue only if nfmodsol > 0

    if (nfmodsol.gt.0) then

      ! We allocate the percentage of soil for all faces, its definition is done th the second
      ! call of usatsoil
      allocate(pourcent_sol(nfmodsol,nbrsol),stat = error)

      if (error /= 0) then
        write(nfecra,*) "Allocation error of atmodsol::pourcent_sol"
        call csexit(1)
      endif

      ! We allocate the structure use for the solving with soil constants for
      ! each soil face and the 3 variables
      allocate(solution_sol(nfmodsol),stat = error)

      if (error /= 0) then
        write(nfecra,*) "Allocation error of atmodsol::solution_sol"
        call csexit(1)
      endif

    endif ! nfmodsol > 0

  endif ! End of the first call

  if (iappel.eq.2.and.nfmodsol.gt.0) then
    call usatsoil(iappel)

    call solmoy(error)
    if (error /= 0) then
      write(nfecra,*) "Allocation error of atmodsol::solmoy"
      call csexit(1)
    endif

    !Initialisation des variables Temps , Tempp , Total Water W1 et W2
    call soliva()

  endif ! End of second call

endif ! iatsoil > 0

!----
! End
!----

return
end subroutine atmsol

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
subroutine atmsol
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
!===============================================================================

! Local variables

integer          iappel
integer          error

!===============================================================================

if (iatsoil.ge.0) then
  ! Premier appel: definition de nfmodsol
  iappel = 1
  call usatsoil(iappel)
  !============

  ! On fabrique une table de valeur des constantes utilisees dans le
  ! modele sol
  allocate(tab_sol(nbrsol),stat = error)
  call solcat( error )

  if (error /= 0) then
    write(nfecra,*) "Allocation error of atmodsol::tab_sol"
    call csexit(1)
  endif

  ! On continue seulement si nfmodsol > 0

  if (nfmodsol.gt.0) then

    ! On alloue le pourcentage de presence de sol pour chaque face de bord
    ! sa definition se fera dans l'appel2 de usatsoil
    allocate(pourcent_sol(nfmodsol,nbrsol),stat = error)

    if (error /= 0) then
      write(nfecra,*) "Allocation error of atmodsol::pourcent_sol"
      call csexit(1)
    endif

    iappel = 2
    call usatsoil(iappel)
    !============

    ! On definit une structure dediee a la resolution du probleme,
    ! avec presence des constantes  propre a chaque face ainsi que
    ! des 3 variables que l'on traitera
    allocate(solution_sol(nfmodsol),stat = error)

    if (error /= 0) then
      write(nfecra,*) "Allocation error of atmodsol::solution_sol"
      call csexit(1)
    endif

    call solmoy( error )
    if (error /= 0) then
      write(nfecra,*) "Allocation error of atmodsol::solmoy"
      call csexit(1)
    endif

    !Initialisation des variables Temps , Tempp , Total Water W1 et W2
    call soliva()

  endif ! nfmodsol > 0

endif ! iatsoil > 0

!----
! End
!----

return
end subroutine atmsol

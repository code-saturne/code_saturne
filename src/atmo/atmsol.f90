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
!> \brief    build constants and variables to describe ground model
!>-     NB : soil model structures defined in module atsoil.f90

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine atmsol( )

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

! Local variables

integer          error, n_g_soil_elts
integer, dimension(:), pointer :: elt_ids

!===============================================================================

! Get the number of element in the soil zone
call atmo_get_soil_zone(nfmodsol, nbrsol, elt_ids)

! Global number for all ranks
n_g_soil_elts = nfmodsol

if (irangp.ge.0) then
  call parcpt(n_g_soil_elts)
endif

! There are some soil faces on some ranks
! Note: we can use soil categories without soil model
! (which solve Temperature and humidity)
if (n_g_soil_elts.gt.0) then
  ! Second pass, print and check soil categories parameters
  call solcat(2)

  call solmoy(error)
  if (error /= 0) then
    write(nfecra,*) "Allocation error of atmodsol::solmoy"
    call csexit(1)
  endif

  ! Initialization of soil variables
  ! Only if soil is activated
  if (iatsoil.eq.1) then
    call soliva()
  endif

endif ! End of second call

!----
! End
!----

return
end subroutine atmsol

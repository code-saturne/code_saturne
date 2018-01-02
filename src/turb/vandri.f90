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


!===============================================================================
! Function:
! ---------
!> \file vandri.f90
!> \brief Imposition of an amortization of Van Driest type for the LES.
!>        \f$ \nu_T \f$ is absorbed by \f$ (1-\exp(\dfrac{-y^+}{d^+}))^2 \f$
!>        where \f$ d^+ \f$ is set at 26.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
! modename        name          role
!______________________________________________________________________________!
!> \param[in]     visvdr        dynamic viscosity in edge cells after
!>                               driest velocity amortization
!> \param[in]     yplusc        \f$ y^+\f$ value in cells
!______________________________________________________________________________!

subroutine vandri &
 (visvdr, yplusc)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use parall
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision visvdr(ncelet)
double precision yplusc(ncelet)

! Local variables

integer          iel
double precision yplus
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct

!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

do iel = 1, ncel
  yplus = yplusc(iel)
  visct(iel) = visct(iel)*(1.0d0-exp(-yplus/cdries))**2
enddo

! For the wall cells we add the turbulent viscosity which was absorbed
! in clptur and which has served to calculate the boundary conditions
do iel = 1, ncel
  if (visvdr(iel).gt.-900.d0) visct(iel) = visvdr(iel)
enddo

!----
! End
!----

return
end subroutine

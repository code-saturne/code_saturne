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

!===============================================================================
! Function:
! ---------

!> \file cffana.f90
!>
!> \brief Computes the analytical flux at the boundary for Euler and Energy
!>
!> The Euler equations used to compute the flux are:
!> \f{eqnarray*}{
!>    \der{\rho}{t} + \divs \left(\rho\vect{u}\right) &=&0
!> \\ \der{\rho \vect{u}}{t} + \divv \left(\vect{u}\otimes\rho\vect{u}\right)&=&0
!> \\ \der{E}{t} + \divs \left(\rho\vect{u} E\right) &=&0
!> \f}
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     imodif        indicator of what is computed
!> \param[in,out] bc_en         dirichlet value for the total energy
!> \param[in,out] bc_pr         dirichlet value for the pressure
!> \param[in,out] bc_vel        dirichlet value for the velocity
!_______________________________________________________________________________


subroutine cffana &
 ( imodif , bc_en  , bc_pr  , bc_vel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use parall
use pointe
use entsor
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          imodif

double precision bc_en(nfabor), bc_pr(nfabor), bc_vel(3,nfabor)

! Local variables

integer          iel    , ifac
integer          ien
double precision und    , rund
double precision, dimension(:,:), pointer :: cofacv
double precision, dimension(:), pointer :: coface
double precision, dimension(:), pointer :: brom

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

ien = isca(ienerg)

call field_get_coefac_v(ivarfl(iu), cofacv)
call field_get_coefac_s(ivarfl(ien), coface)

call field_get_val_s(ibrom, brom)

ifac  = imodif
iel   = ifabor(ifac)

!===============================================================================
! 1. COMPUTE VALUES NEEDED FOR ANALYTICAL FLUX
!===============================================================================

und   = (bc_vel(1,ifac)*surfbo(1,ifac)                          &
       + bc_vel(2,ifac)*surfbo(2,ifac)                          &
       + bc_vel(3,ifac)*surfbo(3,ifac))/surfbn(ifac)
rund  = brom(ifac)*und

!===============================================================================
! 2. CONVECTIVE ANALYTICAL FLUX
!===============================================================================

! Tag the faces where an analytical flux is computed
! The tag will be used in bilsc2 to retrieve the faces where an analytical flux
! has to be imposed
icvfli(ifac) = 1

! Momentum flux (the centered pressure contribution is directly taken into account
! in the pressure BC)
cofacv(1,ifac) = suffbn(ifac) * rund * bc_vel(1,ifac)

cofacv(2,ifac) = suffbn(ifac) * rund * bc_vel(2,ifac)

cofacv(3,ifac) = suffbn(ifac) * rund * bc_vel(3,ifac)

! Total energy flux
coface(ifac) = suffbn(ifac) * (rund * bc_en(ifac) +  &
                               und  * bc_pr(ifac))

return

end subroutine

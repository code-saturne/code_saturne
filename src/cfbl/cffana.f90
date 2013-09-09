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
!> \param[in]     nvar          total number of variables
!> \param[in]     imodif        indicator of what is computed
!> \param[in,out] propfb        physical properties at boundary face centers
!> \param[in,out] bval          dirichlet value for all variables
!_______________________________________________________________________________


subroutine cffana &
 ( nvar   , imodif , propfb , bval   )

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

integer          nvar
integer          imodif

double precision propfb(nfabor,*)
double precision bval(nfabor,nvar)

! Local variables

integer          iel    , ifac
integer          ien
integer          iflmab , ipcrom , ipbrom
double precision und    , rund
double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:,:), pointer :: cofacv
double precision, dimension(:), pointer :: coface

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

ipcrom = ipproc(irom)
ipbrom = ipprob(irom)
ien = isca(ienerg)

call field_get_key_int(ivarfl(ien), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

call field_get_coefac_v(ivarfl(iu), cofacv)
call field_get_coefac_s(ivarfl(ien), coface)

ifac  = imodif
iel   = ifabor(ifac)

!===============================================================================
! 1. COMPUTE VALUES NEEDED FOR ANALYTICAL FLUX
!===============================================================================

und   = (bval(ifac,iu)*surfbo(1,ifac)                          &
       + bval(ifac,iv)*surfbo(2,ifac)                          &
       + bval(ifac,iw)*surfbo(3,ifac))/surfbn(ifac)
rund  = propfb(ifac,ipbrom)*und

!===============================================================================
! 2. CONVECTIVE ANALYTICAL FLUX
!===============================================================================

! Tag the faces where an analytical flux is computed
! The tag will be used in bilsc2 to retrieve the faces where an analytical flux
! has to be imposed
icvfli(ifac) = 1

! Mass flux
bmasfl(ifac) = rund * surfbn(ifac)

! Momentum flux (the centered pressure contribution is directly taken into account
! in the pressure BC)
cofacv(1,ifac) = surfbn(ifac) * rund * bval(ifac,iu)

cofacv(2,ifac) = surfbn(ifac) * rund * bval(ifac,iv)

cofacv(3,ifac) = surfbn(ifac) * rund * bval(ifac,iw)

! Total energy flux
coface(ifac) = surfbn(ifac) * (rund * bval(ifac,ien) +              &
                               und  * bval(ifac,ipr))

return

end subroutine

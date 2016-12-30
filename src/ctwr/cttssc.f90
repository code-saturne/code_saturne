!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine cttssc &
!================

 ( iscal  , ncesmp , icetsm ,                                          &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iscal            ! i  ! <-- ! scalar number                                  !
! ncesmp           ! i  ! <-- ! number of cells with mass source terms         !
! icetsm           ! i  ! <-- ! index of cells with mass source term           !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use period
use ppppar
use ppthch
use ppincl
use ctincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal, ncesmp
integer          icetsm(ncesmp)

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

double precision, dimension (:), allocatable :: ct_smbrs, ct_rovsdt

integer evap_model       ! evaporation model
integer i, iel

!===============================================================================

allocate(ct_smbrs(1:ncesmp),ct_rovsdt(1:ncesmp))

!===============================================================================
! 1. Compute the phase change source terms per unit volume in the packing zone
!===============================================================================

do i = 1, ncesmp
  ct_smbrs(i) = 0d0
  ct_rovsdt(i) = 0d0
end do

call cs_ctwr_source_term(ivarfl(isca(iscal)),            &
                         p0,                             &
                         molmass_rat,                    &
                         ct_smbrs, ct_rovsdt);

!===============================================================================
! 2. Compute the scalar source terms
!===============================================================================

do i = 1, ncesmp
  iel = icetsm(i)
  smbrs(iel) = smbrs(iel) + ct_smbrs(i)*cell_f_vol(iel)
  rovsdt(iel) = rovsdt(iel) + ct_rovsdt(i)*cell_f_vol(iel)
end do

deallocate(ct_smbrs,ct_rovsdt)

!----
! End
!----

return
end subroutine

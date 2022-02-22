!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_extra_operations-parallel_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \ref cs_user_extra_operations_examples and
!> \ref cs_user_extra_operations-nusselt_calculation.c for examples.
!>
!> This is an example of cs_user_extra_operations.f90 which
!> performs parallel operations

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel    , ifac
integer          ii     , nbr    , irangv
integer          itab(3)

double precision rrr
double precision xyz(3)

double precision, dimension(:,:), pointer :: cvar_vel

!< [loc_var_dec]

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), cvar_vel)

!===============================================================================
! Example use of parallel utility functions for several operations
!===============================================================================

! This example demonstrates the parallel utility functions that may be used
! to simplify parallel operations.

! CAUTION: these routines modify their input

! ----------------------------------------------

! Sum of an integer counter 'ii', here the number of cells

!< [example_1]
! local value
ii = ncel
! global sum
if (irangp.ge.0) then
  call parcpt(ii)
endif
! print the global sum
write(nfecra,5020)ii
 5020 format(' cs_user_extra_operations: total number of cells = ', i10)
!< [example_1]

! Maximum of an integer counter 'ii', here the number of cells

!< [example_2]
! local value
ii = ncel
! global maximum
if (irangp.ge.0) then
  call parcmx(ii)
endif
! print the global maximum value
write(nfecra,5010)ii
 5010 format(' cs_user_extra_operations: max. number of cells per process = ',&
            i10)
!< [example_2]

! Sum of a real 'rrr', here the volume

!< [example_3]
! local value
rrr = 0.d0
do iel = 1, ncel
  rrr = rrr + volume(iel)
enddo
! global sum
if (irangp.ge.0) then
  call parsom(rrr)
endif
! print the global sum
write(nfecra,5030)rrr
 5030 format(' cs_user_extra_operations: total domain volume = ', e14.5)
!< [example_3]

! Maximum of a real 'rrr', here the volume

!< [example_4]
! local value
rrr = 0.d0
do iel = 1, ncel
  if (volume(iel).gt.rrr) rrr = volume(iel)
enddo
! global maximum
if (irangp.ge.0) then
  call parmax(rrr)
endif
! print the global maximum
write(nfecra,5040)rrr
 5040 format(' cs_user_extra_operations: max volume per process = ', e14.5)
!< [example_4]

! Minimum of a real 'rrr', here the volume

!< [example_5]
! local value
rrr = grand
do iel = 1, ncel
  if (volume(iel).lt.rrr) rrr = volume(iel)
enddo
! global minimum
if (irangp.ge.0) then
  call parmin(rrr)
endif
! print the global minimum
write(nfecra,5050)rrr
 5050 format(' cs_user_extra_operations: min volume per process = ', e14.5)
!< [example_5]

! Maximum of a real and associated real values;
! here the volume and its location (3 coordinates)

!< [example_6]
nbr = 3
rrr  = -1.d0
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  if (rrr.lt.volume(iel)) then
    rrr = volume(iel)
    xyz(1) = xyzcen(1,iel)
    xyz(2) = xyzcen(2,iel)
    xyz(3) = xyzcen(3,iel)
  endif
enddo
! global maximum and associated location
if (irangp.ge.0) then
  call parmxl(nbr, rrr, xyz)
endif
! print the global maximum and its associated values
write(nfecra,5060) rrr, xyz(1), xyz(2), xyz(3)
 5060 format(' Cs_user_extra_operations: Max. volume =      ', e14.5, /,  &
             '         Location (x,y,z) = ', 3e14.5)
!< [example_6]

! Minimum of a real and associated real values;
! here the volume and its location (3 coordinates)

!< [example_7]
nbr = 3
rrr  = 1.d+30
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  if (rrr.gt.volume(iel)) then
    rrr = volume(iel)
    xyz(1) = xyzcen(1,iel)
    xyz(2) = xyzcen(2,iel)
    xyz(3) = xyzcen(3,iel)
  endif
enddo
! global minimum and associated location
if (irangp.ge.0) then
  call parmnl(nbr,rrr,xyz)
endif
! print the global minimum and its associated values
write(nfecra,5070) rrr, xyz(1), xyz(2), xyz(3)
 5070 format(' Cs_user_extra_operations: Min. volume =      ', e14.5, /,  &
             '         Location (x,y,z) = ', 3e14.5)
!< [example_7]

! Sum of an array of integers;
! here, the number of cells, faces, and boundary faces

! local values; note that to avoid counting interior faces on
! parallel boundaries twice, we check if 'ifacel(1,ifac) .le. ncel',
! as on a parallel boundary, this is always true for one domain
! and false for the other.

!< [example_8]
nbr = 3
itab(1) = ncel
itab(2) = 0
itab(3) = nfabor
do ifac = 1, nfac
  if (ifacel(1, ifac).le.ncel) itab(2) = itab(2) + 1
enddo
! global sum
if (irangp.ge.0) then
  call parism(nbr, itab)
endif
! print the global sums
write(nfecra,5080) itab(1), itab(2), itab(3)
 5080 format(' cs_user_extra_operations: Number of cells =         ', i10, /,  &
             '         Number of interior faces = ', i10, /,  &
             '         Number of boundary faces = ', i10)
!< [example_8]

! Maxima from an array of integers;
! here, the number of cells, faces, and boundary faces

!< [example_9]
! local values
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
! global maxima
if (irangp.ge.0) then
  call parimx(nbr, itab)
endif
! print the global maxima
write(nfecra,5090) itab(1), itab(2), itab(3)
 5090 format(' cs_user_extra_operations: Max. number of cells per rank = ',  &
             i10, /,  &
             '         Max. number of interior faces per rank = ', i10, /,   &
             '         Max. number of boundary faces per rank = ', i10)
!< [example_9]

! Minima from an array of integers;
! here, the number of cells, faces, and boundary faces

!< [example_10]
! local values
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
! global minima
if (irangp.ge.0) then
  call parimn(nbr, itab)
endif
! print the global minima
write(nfecra,5100) itab(1), itab(2), itab(3)
 5100 format(' cs_user_extra_operations: Min. number of cells per rank = ',  &
             i10, /,  &
             '         Min. number of interior faces per rank = ', i10, /,   &
             '         Min. number of boundary faces per rank = ', i10)
!< [example_10]

! Sum of an array of reals;
! here, the 3 velocity components (so as to compute a mean for example)

! local values
!< [example_11]
nbr = 3
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  xyz(1) = xyz(1) + cvar_vel(1,iel)
  xyz(2) = xyz(2) + cvar_vel(2,iel)
  xyz(3) = xyz(3) + cvar_vel(3,iel)
enddo
! global sum
if (irangp.ge.0) then
  call parrsm(nbr, xyz)
endif
! print the global sums
write(nfecra,5110) xyz(1), xyz(2), xyz(3)
 5110 format(' cs_user_extra_operations: Sum of U on the domain = ', e14.5, /, &
             '         Sum of V on the domain = ', e14.5, /,   &
             '         Sum of V on the domain = ', e14.5)
!< [example_11]

! Maximum of an array of reals;
! here, the 3 velocity components

! local values
!< [example_12]
nbr = 3
xyz(1) = cvar_vel(1,1)
xyz(2) = cvar_vel(2,1)
xyz(3) = cvar_vel(3,1)
do iel = 1, ncel
  xyz(1) = max(xyz(1),cvar_vel(1,iel))
  xyz(2) = max(xyz(2),cvar_vel(2,iel))
  xyz(3) = max(xyz(3),cvar_vel(3,iel))
enddo
! global maximum
if (irangp.ge.0) then
  call parrmx(nbr, xyz)
endif
! print the global maxima
write(nfecra,5120) xyz(1), xyz(2), xyz(3)
 5120 format(' cs_user_extra_operations: Maximum of U on the domain = ',  &
             e14.5, /,   &
             '         Maximum of V on the domain = ', e14.5, /,   &
             '         Maximum of V on the domain = ', e14.5)
!< [example_12]

! Maximum of an array of reals;
! here, the 3 velocity components

!< [example_13]
! local values
nbr = 3
xyz(1) = cvar_vel(1,1)
xyz(2) = cvar_vel(2,1)
xyz(3) = cvar_vel(3,1)
do iel = 1, ncel
  xyz(1) = min(xyz(1),cvar_vel(1,iel))
  xyz(2) = min(xyz(2),cvar_vel(2,iel))
  xyz(3) = min(xyz(3),cvar_vel(3,iel))
enddo
! global minimum
if (irangp.ge.0) then
  call parrmn(nbr, xyz)
endif
! print the global maxima
write(nfecra,5130) xyz(1), xyz(2), xyz(3)
 5130 format(' cs_user_extra_operations: Minimum of U on the domain = ', &
             e14.5, /,   &
             '         Minimum of V on the domain = ', e14.5, /,   &
             '         Minimum of V on the domain = ', e14.5)
!< [example_13]

! Broadcast an array of local integers to other ranks;
! in this example, we use the number of cells, interior faces, and boundary
! faces from process rank 0 (irangv).

!< [example_14]
! local values
irangv = 0
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
! broadcast from rank irangv to all others
if (irangp.ge.0) then
  call parbci(irangv, nbr, itab)
endif
! print values broadcast and received from rank 'irangv'
write(nfecra,5140) irangv, itab(1), itab(2), itab(3)
 5140 format(' cs_user_extra_operations: On rank ', i10 , /,   &
             '         Number of cells          = ', i10, /,   &
             '         Number of interior faces = ', i10, /,   &
             '         Number of boundary faces = ', i10)
!< [example_14]

! Broadcast an array of local reals to other ranks;
! in this example, we use 3 velocity values from process rank 0 (irangv).

!< [example_15]
! local values
irangv = 0
nbr = 3
xyz(1) = cvar_vel(1,1)
xyz(2) = cvar_vel(2,1)
xyz(3) = cvar_vel(3,1)
! broadcast from rank irangv to all others
if (irangp.ge.0) then
  call parbcr(irangv, nbr, xyz)
endif
! print values broadcast and received from rank 'irangv'
write(nfecra,5150) irangv, xyz(1), xyz(2), xyz(3)
 5150 format(' cs_user_extra_operations: On rank ', i10 , /,     &
             '         Velocity U in first cell = ', e14.5, /,   &
             '         Velocity V in first cell = ', e14.5, /,   &
             '         Velocity W in first cell = ', e14.5)
!< [example_15]

return
end subroutine cs_f_user_extra_operations

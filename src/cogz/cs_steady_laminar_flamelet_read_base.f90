!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file cs_steady_laminar_flamelet_read_base.f90
!>
!> \brief Specific physic subroutine: gas combustion diffusion flames
!>
!> Read data file Turbulent_Library for steady laminar flamelet model
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine cs_steady_laminar_flamelet_read_base () &
  bind(C, name='cs_f_steady_laminar_flamelet_read_base')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppthch
use coincl
use ppincl
use radiat

!=========================================================

implicit none

call init_steady_laminar_flamelet_library

! Lecture de librairie flamelettes turbulentes
call read_flamelet_library

! Lecture de librairie flamelettes radiatives
if (iirayo .eq. 1) call read_radiation_library

return

end subroutine

!=====================================================================
!   read flamelet library
!=====================================================================

subroutine read_flamelet_library()

use entsor
use coincl

integer ::  ii, jj, mm, nn

character(len=64) :: ficfpp
integer(c_int) :: c_name_max, c_name_len
type(c_ptr) :: c_name_p
character(kind=c_char, len=1), dimension(:), pointer :: c_name

!===============================================================================
! Interfaces
!===============================================================================

  interface

    subroutine cs_f_combustion_gas_get_data_file_name  &
      (f_name_max, f_name, f_name_len)                 &
      bind(C, name='cs_f_combustion_gas_get_data_file_name')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_name_max
      type(c_ptr), intent(out)    :: f_name
      integer(c_int), intent(out) :: f_name_len
    end subroutine cs_f_combustion_gas_get_data_file_name

  end interface

!========================================================

! Get thermochemistry data file name

c_name_max = len(ficfpp)
call cs_f_combustion_gas_get_data_file_name(c_name_max, c_name_p, c_name_len)
call c_f_pointer(c_name_p, c_name, [c_name_len])

do ii = 1, c_name_len
  ficfpp(ii:ii) = c_name(ii)
enddo
do ii = c_name_len + 1, c_name_max
  ficfpp(ii:ii) = ' '
enddo

open(unit = impfpp, file = ficfpp, status='old')

! skip the 2 head lines

read(impfpp,*)
read(impfpp,*)

do mm = 1, nki
  do nn = 1, nxr
    do jj = 1, nzvar
      read(impfpp,*) ! skip the 1 seperation line in .Dat
      do ii = 1, nzm
        read(impfpp,*) flamelet_library(:,nn,mm,jj,ii)
      enddo
    enddo
  enddo
enddo

close(impfpp)

end subroutine

!=====================================================================
!   read radiation library
!=====================================================================

subroutine read_radiation_library()

use entsor
use coincl

implicit none

integer ::  ii, jj, mm, nn,ig
integer :: imprad

imprad = -1
open(newunit=imprad, file=ficrad, status='old')

do mm = 1, nki
  do nn = 1, nxr
    do jj = 1, nzvar
      do ii = 1, nzm
        do ig=1, nwsgg
          read(imprad,*) radiation_library(:,ig,nn,mm,jj,ii)
        enddo
      enddo
    enddo
  enddo
enddo

close(imprad)

end subroutine



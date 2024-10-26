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

subroutine csopli                    &
 (infecr, isuppr, ierror)            &
  bind(C, name='cs_f_open_run_log')

!===============================================================================
! Purpose:
! -------

!    Open log files using Fortran IO.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! infecr           ! i  ! <-- ! value to assign to nfecra                      !
! isuppr           ! i  ! <-- ! supress output if ~                            !
! ierror           ! i  ! --> ! error code                                     !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor

!===============================================================================

use iso_c_binding
implicit none

! Arguments

integer(c_int) :: infecr, isuppr, ierror, i

! Local variables

character(kind=c_char, len=1), dimension(64) :: c_path
character(len=64) :: name, t_name

!===============================================================================

  interface

    subroutine cs_f_base_log_name(lmax, path)     &
      bind(C, name='cs_f_base_log_name')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: lmax
      character(kind=c_char, len=1), dimension(*), intent(out) :: path
    end subroutine cs_f_base_log_name

  end interface

!===============================================================================

ierror = 0

nfecra = infecr

if (nfecra .eq. 6) return

call cs_f_base_log_name(64, c_path)
t_name = " "
do i = 1, 64
  if (c_path(i) == c_null_char) exit
  t_name(i:i) = c_path(i)
enddo
name = trim(t_name)

if (isuppr .eq. 0) then
  open(file=name, unit=nfecra, form='formatted', status='old',   &
       position='append', action='write', err=900)
else
  open(file=name, unit=nfecra, form='formatted', status='unknown', err=900)
endif

goto 950

900 ierror = 1

950 continue

return
end subroutine

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

!> \file entsor.f90
!> \brief Module for input/output

module entsor

  !=============================================================================
  use, intrinsic :: iso_c_binding
  use paramx

  implicit none

  !=============================================================================

  !> \defgroup entsor Module for input/output

  !> \addtogroup entsor
  !> \{

  !> standard output
  integer, save :: nfecra

  !> field key for output label (\ref label field keyword).
  integer, save :: keylbl = -1

  !> field key for logging (\ref log field keyword).
  integer, save :: keylog = -1

  !> <a name="keyvis"></a>
  !> field key for postprocessing output (\ref post_vis field keyword).
  integer, save :: keyvis = -1

  !> \}

  !> \defgroup userfile Additional user files

  !> \addtogroup userfile
  !> \{

  !> name of the thermochemical data file for combustion.
  !>
  !> Useful in case of gas combustion.
  character(len=64), save :: ficrad

  !> logical unit of the thermochemical data file.
  !> Useful in case of gas combustion;
  integer, save :: impfpp = 25

  !> \}

  !=============================================================================

contains

  !=============================================================================

  !> \brief Flush Fortran log

  subroutine flush_nfecra() &
    bind(C, name='cs_f_flush_logs')
    flush(nfecra)
  end subroutine flush_nfecra

  !=============================================================================

  !> \brief Open log files using Fortran IO.

  !> \param[in]   infecr   value to assign to nfecra
  !> \param[in]   isuppr   supress output if ~
  !> \param[out]  ierror   error code

  subroutine csopli(infecr, isuppr, ierror) &
    bind(C, name='cs_f_open_run_log')

    use iso_c_binding
    implicit none

    ! Arguments

    integer(c_int) :: infecr, isuppr, ierror, i

    ! Local variables

    character(kind=c_char, len=1), dimension(64) :: c_path
    character(len=64) :: name, t_name

    interface
      subroutine cs_f_base_log_name(lmax, path)     &
        bind(C, name='cs_f_base_log_name')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), value :: lmax
        character(kind=c_char, len=1), dimension(*), intent(out) :: path
      end subroutine cs_f_base_log_name
    end interface

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
  end subroutine csopli

  !=============================================================================

  !> \brief Close log files using Fortran IO.

  subroutine csclli() bind(C, name='cs_f_close_run_log')

    implicit none

    ! If output has been redirected, it uses unit nfecra = 9 instead of 6

    if (nfecra.ne.6) then
      close(nfecra)
      nfecra = 6
    endif

  end subroutine csclli

  !=============================================================================

  !> \brief Log a character string

  !> \param[in]   str      character string
  !> \param[in]   l        string length

  subroutine csprnt(str, l) &
     bind(C, name='cs_f_print')

    implicit none

    ! Arguments

    character      :: str(*)
    integer(c_int) :: l

    ! Local variables

    character     chloc*16384
    integer       ii

    l = min(l, 16384 - 1)

    do ii = 1, l
      chloc(ii:ii) = str(ii)
    enddo

    write(nfecra, '(a)', advance='no') chloc(1:l)

    return

  end subroutine csprnt

  !=============================================================================

end module entsor

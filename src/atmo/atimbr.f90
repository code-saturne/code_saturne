!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!==============================================================================
! Function:
! ---------

!> \file atimbr.f90
!>
!> \brief Atmospheric Imbrication module.
!>  This module contains the data structure and subroutines to
!>  perform atmospheric imbrication or nesting of a CFD domain within
!>  a large scale meteorological field.
!>  Starting from a set of large  scale meteological profiles (in the
!>  format of meteo files) an interpolation is performed for each
!>  boundary face both spatially and temporally (using Cressman method)
!>
!>  This imbrication reads additional meteo files.
!______________________________________________________________________________

module atimbr
!> \defgroup at_imbrication

use entsor
use atincl
implicit none

!> \addtogroup at_imbrication
!> \{
! --------------------------------------------------------------
!> activation flag
! --------------------------------------------------------------
logical imbrication_flag
logical imbrication_verbose
save imbrication_verbose

! --------------------------------------------------------------
!> Flags for activating the cressman interpolation for the boundary
!> conditions
! --------------------------------------------------------------
logical cressman_u
logical cressman_v
logical cressman_tke
logical cressman_eps
logical cressman_theta
logical cressman_qw
logical cressman_nc

! --------------------------------------------------------------
!> numerical parameters for the cressman interpolation formulas
! --------------------------------------------------------------
double precision :: horizontal_influence_radius
double precision :: vertical_influence_radius

! --------------------------------------------------------------
!> Parameter for "meteo" files
! --------------------------------------------------------------
integer line_len
parameter(line_len = 132)
character(line_len) :: imbrication_files_list
character(line_len), dimension(:), allocatable :: imbrication_files
integer number_of_files
character*(3) skip_chars
parameter (skip_chars = "/#!")
!> Profile dimension variable
integer thermal_profile_dim
integer dynamical_profile_dim
!> Time sections per files
integer sections_per_file
data thermal_profile_dim /-1/
data dynamical_profile_dim /-1/
data sections_per_file /-1/

! --------------------------------------------------------------
!> read data from "meteo" files
! --------------------------------------------------------------
!> Time variables
integer, dimension(:,:), allocatable :: years
integer, dimension(:,:), allocatable :: ordinals
integer, dimension(:,:), allocatable :: hours
integer, dimension(:,:), allocatable :: minutes
double precision, dimension(:,:),allocatable :: seconds
!> Positions
double precision, dimension(:,:), allocatable :: xpos
double precision, dimension(:,:), allocatable :: ypos
double precision, dimension(:,:), allocatable :: ground_pressure
!> Vertical grid for temperature and humidity variables
double precision, dimension(:,:,:), allocatable :: zt
double precision, dimension(:,:,:), allocatable :: tempC
double precision, dimension(:,:,:), allocatable :: qw
double precision, dimension(:,:,:), allocatable :: Nc
!> Vertical grid for wind variables
double precision, dimension(:,:,:), allocatable :: zd
double precision, dimension(:,:,:), allocatable :: u
double precision, dimension(:,:,:), allocatable :: v
double precision, dimension(:,:,:), allocatable :: tke
double precision, dimension(:,:,:), allocatable :: eps

! --------------------------------------------------------------
!> derived data
! --------------------------------------------------------------
double precision, dimension(:,:), allocatable,target :: times
double precision, dimension(:,:,:), allocatable :: pressure
double precision, dimension(:,:,:), allocatable :: theta
double precision, dimension(:,:,:), allocatable :: density

! --------------------------------------------------------------
!> time interpolated profiles
! --------------------------------------------------------------
double precision, dimension(:,:), allocatable :: ti_zt
double precision, dimension(:,:), allocatable :: ti_tempC
double precision, dimension(:,:), allocatable :: ti_qw
double precision, dimension(:,:), allocatable :: ti_Nc
double precision, dimension(:,:), allocatable :: ti_zd
double precision, dimension(:,:), allocatable :: ti_u
double precision, dimension(:,:), allocatable :: ti_v
double precision, dimension(:,:), allocatable :: ti_tke
double precision, dimension(:,:), allocatable :: ti_eps
double precision, dimension(:,:), allocatable :: ti_pressure
double precision, dimension(:,:), allocatable :: ti_theta
double precision, dimension(:,:), allocatable :: ti_density

! --------------------------------------------------------------
!> additional variables
! --------------------------------------------------------------
double precision, dimension(:,:,:), allocatable :: coordinates_th
double precision, dimension(:,:,:), allocatable :: influence_param_th
double precision, dimension(:,:,:), allocatable :: coordinates_dyn
double precision, dimension(:,:,:), allocatable :: influence_param_dyn
integer id_u
integer id_v
integer id_tke
integer id_eps
integer id_theta
integer id_qw
integer id_nc

! --------------------------------------------------------------
!> 1D array of times at which profiles are given
! --------------------------------------------------------------
double precision, dimension(:), pointer :: times_sequence=>null()
!> \}
contains

! ----------------------------------------------------------------
!> \brief Allocate variables adapted to the number of files
!>  and time step to be considered
! ----------------------------------------------------------------
subroutine allocate_all()
implicit none
logical error
data error /.false./
save error

! --- section of explicit allocations -----------------
allocate(years(sections_per_file,number_of_files))
allocate(ordinals(sections_per_file,number_of_files))
allocate(hours(sections_per_file,number_of_files))
allocate(minutes(sections_per_file,number_of_files))
allocate(seconds(sections_per_file,number_of_files))
allocate(xpos(sections_per_file,number_of_files))
allocate(ypos(sections_per_file,number_of_files))
allocate(ground_pressure(sections_per_file,number_of_files))
allocate(zt(thermal_profile_dim,sections_per_file,number_of_files))

! --------------------------------------------------------------
!tempC and qw are always allocated for dry and humid atmosphere
! --------------------------------------------------------------
if (ippmod(iatmos).ge.0) then
  allocate(tempC(thermal_profile_dim,sections_per_file,number_of_files))
  allocate(qw(thermal_profile_dim,sections_per_file,number_of_files))
endif

! --------------------------------------------------------------
! Nc only used in humid atmosphere
! --------------------------------------------------------------
if (ippmod(iatmos).ge.2) then
  allocate(Nc(thermal_profile_dim,sections_per_file,number_of_files))
endif

allocate(zd(dynamical_profile_dim,sections_per_file,number_of_files))
allocate(u(dynamical_profile_dim,sections_per_file,number_of_files))
allocate(v(dynamical_profile_dim,sections_per_file,number_of_files))
allocate(tke(dynamical_profile_dim,sections_per_file,number_of_files))
allocate(eps(dynamical_profile_dim,sections_per_file,number_of_files))
end subroutine allocate_all

! ----------------------------------------------------------------
!> \brief  Final step for deallocation
! ----------------------------------------------------------------
subroutine finalize_imbrication()
implicit none
! --- section of explicit deallocations -----------------
deallocate(imbrication_files)
deallocate(years)
deallocate(ordinals)
deallocate(hours)
deallocate(minutes)
deallocate(seconds)
deallocate(xpos)
deallocate(ypos)
deallocate(ground_pressure)
deallocate(zt)
if (ippmod(iatmos).ge.0) then
  deallocate(tempC)
  deallocate(qw)
endif
if (ippmod(iatmos).ge.2) then
  deallocate(Nc)
endif
deallocate(zd)
deallocate(u)
deallocate(v)
deallocate(tke)
deallocate(eps)
deallocate(times)
deallocate(pressure)
deallocate(theta)
deallocate(density)
deallocate(ti_zt)
if (ippmod(iatmos).ge.0) then
  deallocate(ti_tempC)
  deallocate(ti_qw)
endif
if (ippmod(iatmos).ge.2) then
  deallocate(ti_Nc)
endif
deallocate(ti_zd)
deallocate(ti_u)
deallocate(ti_v)
deallocate(ti_tke)
deallocate(ti_eps)
deallocate(ti_pressure)
deallocate(ti_theta)
deallocate(ti_density)
deallocate(coordinates_th)
deallocate(influence_param_th)
deallocate(coordinates_dyn)
deallocate(influence_param_dyn)
end subroutine finalize_imbrication

! ---------------------------------------------------------------------------
!> \brief Time interpolation of all profiles -
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   the_time    current time
!-------------------------------------------------------------------------------
subroutine interpolate_all_profiles(the_time)
implicit none
double precision the_time
integer lb,ub
integer i
double precision :: time_was
save time_was
logical first_call
data first_call /.true./
save first_call
if (first_call) then
  time_was = the_time
  first_call = .false.
else
  if (time_was.eq.the_time) then
    if(imbrication_verbose) &
         write(nfecra,*) &
              "interpolate_all_profiles:the_time==time_was==", the_time
    return
  endif
endif
if (imbrication_verbose) &
     write(nfecra,*) "interpolate_all_profiles:the_time==", the_time
if (imbrication_verbose) &
     write(nfecra,*) "interpolate_all_profiles:time_was==", time_was
time_was = the_time

! -----------
!  -   zt    -
!  -----------
if (allocated(zt)) then
  lb = lbound(zt(:,1,1),1)
  ub = ubound(zt(:,1,1),1)
  if (allocated(ti_zt)) then
    if( lb.ne.lbound(ti_zt,1) .or. ub.ne.ubound(ti_zt,1)) then
      deallocate(ti_zt)
    endif
  endif
  if (.not.allocated(ti_zt)) then
    allocate(ti_zt(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),zt(:,:,i),ti_zt(:,i))
  enddo
endif

!  -----------
!  -   tempC    -
!  -----------
if (allocated(tempC)) then
  lb = lbound(tempC(:,1,1),1)
  ub = ubound(tempC(:,1,1),1)
  if (allocated(ti_tempC))then
    if (lb.ne.lbound(ti_tempC,1) .or. ub.ne.ubound(ti_tempC,1)) then
      deallocate(ti_tempC)
    endif
  endif
  if (.not.allocated(ti_tempC))then
    allocate(ti_tempC(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),tempC(:,:,i),ti_tempC(:,i))
  enddo
endif

!  -----------
!  -   qw    -
!  -----------
if (allocated(qw)) then
  lb = lbound(qw(:,1,1),1)
  ub = ubound(qw(:,1,1),1)
  if (allocated(ti_qw)) then
    if (lb.ne.lbound(ti_qw,1) .or. ub.ne.ubound(ti_qw,1)) then
      deallocate(ti_qw)
    endif
  endif
  if (.not.allocated(ti_qw)) then
    allocate(ti_qw(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),qw(:,:,i),ti_qw(:,i))
  enddo
endif

!  -----------
!  -   Nc    -
!  -----------
if (allocated(Nc)) then
  lb = lbound(Nc(:,1,1),1)
  ub = ubound(Nc(:,1,1),1)
  if (allocated(ti_Nc)) then
    if (lb.ne.lbound(ti_Nc,1) .or. ub.ne.ubound(ti_Nc,1)) then
      deallocate(ti_Nc)
    endif
  endif
  if (.not.allocated(ti_Nc)) then
    allocate(ti_Nc(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),Nc(:,:,i),ti_Nc(:,i))
  enddo
endif

!  -----------
!  -   zd    -
!  -----------
if (allocated(zd)) then
  lb = lbound(zd(:,1,1),1)
  ub = ubound(zd(:,1,1),1)
  if (allocated(ti_zd)) then
    if (lb.ne.lbound(ti_zd,1) .or. ub.ne.ubound(ti_zd,1)) then
      deallocate(ti_zd)
    endif
  endif
  if (.not.allocated(ti_zd)) then
    allocate(ti_zd(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),zd(:,:,i),ti_zd(:,i))
  enddo
endif

!  -----------
!  -   u    -
!  -----------
if (allocated(u)) then
  lb = lbound(u(:,1,1),1)
  ub = ubound(u(:,1,1),1)
  if (allocated(ti_u)) then
    if (lb.ne.lbound(ti_u,1) .or. ub.ne.ubound(ti_u,1)) then
      deallocate(ti_u)
    endif
  endif
  if (.not.allocated(ti_u)) then
    allocate(ti_u(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),u(:,:,i),ti_u(:,i))
  enddo
endif

!  -----------
!  -   v    -
!  -----------
if (allocated(v)) then
  lb = lbound(v(:,1,1),1)
  ub = ubound(v(:,1,1),1)
  if (allocated(ti_v)) then
    if (lb.ne.lbound(ti_v,1) .or. ub.ne.ubound(ti_v,1)) then
      deallocate(ti_v)
    endif
  endif
  if (.not.allocated(ti_v)) then
    allocate(ti_v(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),v(:,:,i),ti_v(:,i))
  enddo
endif

!  -----------
!  -   tke    -
!  -----------
if (allocated(tke)) then
  lb = lbound(tke(:,1,1),1)
  ub = ubound(tke(:,1,1),1)
  if (allocated(ti_tke)) then
    if( lb.ne.lbound(ti_tke,1) .or. ub.ne.ubound(ti_tke,1))then
      deallocate(ti_tke)
    endif
  endif
  if (.not.allocated(ti_tke)) then
    allocate(ti_tke(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),tke(:,:,i),ti_tke(:,i))
  enddo
endif

!  -----------
!  -   eps    -
!  -----------
if (allocated(eps)) then
  lb = lbound(eps(:,1,1),1)
  ub = ubound(eps(:,1,1),1)
  if (allocated(ti_eps)) then
    if (lb.ne.lbound(ti_eps,1) .or. ub.ne.ubound(ti_eps,1)) then
      deallocate(ti_eps)
    endif
  endif
  if (.not.allocated(ti_eps)) then
    allocate(ti_eps(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),eps(:,:,i),ti_eps(:,i))
  enddo
endif

!  -----------
!  -   pressure    -
!  -----------
if (allocated(pressure)) then
  lb = lbound(pressure(:,1,1),1)
  ub = ubound(pressure(:,1,1),1)
  if (allocated(ti_pressure)) then
    if (lb.ne.lbound(ti_pressure,1) .or. ub.ne.ubound(ti_pressure,1)) then
      deallocate(ti_pressure)
    endif
  endif
  if (.not.allocated(ti_pressure)) then
    allocate(ti_pressure(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),pressure(:,:,i),      &
                            ti_pressure(:,i))
  enddo
endif

!  -----------
!  -   theta     -
!  -----------
if (allocated(theta)) then
  lb = lbound(theta (:,1,1),1)
  ub = ubound(theta (:,1,1),1)
  if (allocated(ti_theta)) then
    if (lb.ne.lbound(ti_theta ,1) .or. ub.ne.ubound(ti_theta ,1)) then
      deallocate(ti_theta)
    endif
  endif
  if (.not.allocated(ti_theta)) then
    allocate(ti_theta (lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),theta (:,:,i),ti_theta (:,i))
  enddo
endif

!  -----------
!  -   density    -
!  -----------
if (allocated(density)) then
  lb = lbound(density(:,1,1),1)
  ub = ubound(density(:,1,1),1)
  if (allocated(ti_density)) then
    if (lb.ne.lbound(ti_density,1) .or. ub.ne.ubound(ti_density,1)) then
      deallocate(ti_density)
    endif
  endif
  if (.not.allocated(ti_density)) then
    allocate(ti_density(lb:ub,number_of_files))
  endif
  do i = 1, number_of_files, 1
    call time_interpolation(the_time,times(:,1),density(:,:,i),       &
                            ti_density(:,i))
  enddo
endif

end subroutine interpolate_all_profiles

! ----------------------------------------------------------------
!> \brief Print the interpolated profiles for checking purposes
! ----------------------------------------------------------------
subroutine dump_interpolated_profiles
implicit none
integer lb, ub
integer i, j

!TODO: this prints could be merged with a function
! ---------------------------
!  zt
! ---------------------------
if (allocated(ti_zt)) then
  lb = lbound(ti_zt,1)
  ub = ubound(ti_zt,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "zt=", ti_zt(j,i)
    enddo
  enddo
endif

! ---------------------------
!  tempC
! ---------------------------
if (allocated(ti_tempC)) then
  lb = lbound(ti_tempC,1)
  ub = ubound(ti_tempC,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "tempC=", ti_tempC(j,i)
    enddo
  enddo
endif

! ---------------------------
!  qw
! ---------------------------
if (allocated(ti_qw)) then
  lb = lbound(ti_qw,1)
  ub = ubound(ti_qw,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "qw=", ti_qw(j,i)
    enddo
  enddo
endif

! ---------------------------
!  Nc
! ---------------------------
if (allocated(ti_Nc)) then
  lb = lbound(ti_Nc,1)
  ub = ubound(ti_Nc,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "Nc=", ti_Nc(j,i)
    enddo
  enddo
endif

! ---------------------------
!  zd
! ---------------------------
if (allocated(ti_zd)) then
  lb = lbound(ti_zd,1)
  ub = ubound(ti_zd,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "zd=", ti_zd(j,i)
    enddo
  enddo
endif

! ---------------------------
!  u
! ---------------------------
if (allocated(ti_u)) then
  lb = lbound(ti_u,1)
  ub = ubound(ti_u,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "u=", ti_u(j,i)
    enddo
  enddo
endif

! ---------------------------
!  v
! ---------------------------
if (allocated(ti_v)) then
  lb = lbound(ti_v,1)
  ub = ubound(ti_v,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "v=", ti_v(j,i)
    enddo
  enddo
endif

! ---------------------------
!  tke
! ---------------------------
if (allocated(ti_tke)) then
  lb = lbound(ti_tke,1)
  ub = ubound(ti_tke,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "tke=", ti_tke(j,i)
    enddo
  enddo
endif

! ---------------------------
!  eps
! ---------------------------
if (allocated(ti_eps)) then
  lb = lbound(ti_eps,1)
  ub = ubound(ti_eps,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "eps=", ti_eps(j,i)
    enddo
  enddo
endif

! ---------------------------
!  pressure
! ---------------------------
if (allocated(ti_pressure)) then
  lb = lbound(ti_pressure,1)
  ub = ubound(ti_pressure,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "pressure=", ti_pressure(j,i)
    enddo
  enddo
endif

! ---------------------------
!  theta
! ---------------------------
if (allocated(ti_theta)) then
  lb = lbound(ti_theta ,1)
  ub = ubound(ti_theta ,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "theta=", ti_theta(j,i)
    enddo
  enddo
endif

! ---------------------------
!  density
! ---------------------------
if (allocated(ti_density)) then
  lb = lbound(ti_density,1)
  ub = ubound(ti_density,1)
  do i = 1, number_of_files, 1
    write(nfecra,*) imbrication_files(i)
    do j = lb, ub, 1
      write(nfecra,*) "j=", j, "density=", ti_density(j,i)
    enddo
  enddo
endif

end subroutine dump_interpolated_profiles


! ----------------------------------------------------------------
!> \brief  Converts a (year,ordinal) date to julian calendar date
!> for calculating time shifts
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   year
!> \param[in]   ordinal     number of the day in the year
!>                          e.g 1st january has ordinal 1
!>                               31 december 365 or 366
! ----------------------------------------------------------------
integer function yo2j(year,ordinal)
!-------------------------------------------------------------------------------
!
implicit none
integer year,ordinal
! I am puzzled by the (1-14)/12 : why not write -1?

! yo2j= ordinal + ((1461 * (year + 4800 + (1 - 14) / 12)) / 4 +        &
!                  (367 * (1 - 2 - 12 * ((1 - 14) / 12))) / 12 -       &
!                  (3 * ((year + 4900 + (1 - 14) / 12) / 100)) / 4     &
!                   + 1 - 32075) - 1

yo2j= ordinal + ((1461 * (year + 4800)) / 4     &
      -30 - (3 * ((year + 4900) / 100)) / 4     &
      + 1 - 32075) - 1
end function yo2j


! ----------------------------------------------------------------
!> \brief Reads a file having in each significative line a file name
!>  it returns then as 'the_list' the list of lines read
!>  a line is significative if it's first char is not / or # or !
!>  The following 3 lines give an example from which one must remove the
!>  first two characters.
!>- /list of files
!>- profile_one.txt
!>- profile_two.txt
!>
!> Beware that blank lines or lines starting with blanks+ comment_char are
!> NOT ignored.
! ----------------------------------------------------------------
! Arguments
! ----------------------------------------------------------------
!> \param[in]   a_file      the file with list of file names
!> \param[out]  the_list    the list of file names
! ----------------------------------------------------------------
subroutine read_files_list(a_file,the_list)

implicit none
! ----------------------------------------------------------------
! declarations
! ----------------------------------------------------------------
character(line_len) :: a_file
character(line_len), dimension(:), allocatable :: the_list
integer unilog
parameter (unilog = 10)
character(line_len) :: current_line
integer l_iostat, l_stat
integer counter

! ----------------------------------------------------------------
! executables
! ----------------------------------------------------------------
open(unilog,file=imbrication_files_list,status='old', &
     form='formatted',iostat=l_iostat)

counter = 0
do while (.true.)
  call find_next_line(unilog,current_line,a_file,l_iostat)
  if (l_iostat.ne.0) exit
  counter = counter + 1
enddo
number_of_files = counter

allocate(the_list(number_of_files),stat=l_stat)
! find_next closes the file when reaching the end
open(unilog,file=imbrication_files_list,status='old',&
     form='formatted',iostat=l_iostat)
counter = 0
do while (.true.)
  call find_next_line(unilog,current_line,a_file,l_iostat)
  if(l_iostat.ne.0) exit
  counter = counter + 1
the_list(counter) = current_line
enddo
end subroutine read_files_list


! ----------------------------------------------------------------
!> \brief  Find next validated line
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   unilog          logical unit number of the reading file
!> \param[out]  current_line    the characters of the line
!> \param[in]   meteo_file      name of the 'meteo' file
!> \param[out]  l_iostat        logical status of I/O following a read statement
!-------------------------------------------------------------------------------
subroutine find_next_line(unilog,current_line,meteo_file,l_iostat)
implicit none
integer unilog
character(line_len) :: current_line
character(line_len) :: meteo_file
integer l_iostat
integer first,last

next_line:do while (.true.)
  read (unilog,'(A132)',iostat=l_iostat) current_line
  if (l_iostat.gt.0) then
    call bounds(meteo_file,len(meteo_file),first,last)
    write(nfecra,*)"unexpected read error (1) on file ",meteo_file(first:last)
    write(nfecra,*)"connected logical unit :",unilog
    call bounds(current_line,len(current_line),first,last)
    write(nfecra,*)"current_line is (was?):>",current_line(1:last),"<"
    stop
  endif
  if (l_iostat.lt.0) then
    close(unilog)
    return
  endif
  call bounds(current_line,len(current_line),first,last)
  if (first.le.last .and.                                  &
      verify(current_line(first:first),skip_chars).ne.0) exit
enddo next_line
end subroutine find_next_line


! ----------------------------------------------------------------
!> \brief Reads a meteo_file for Code_Saturne Atmospheric Physics option
!>
!> They contain an arbitrary number (>=1) of sections having the following
!> structure.
!> Comment lines start with a slash / as first character
!>
!>- yyyy,dd,hh,mm,ss
!>- xpos,ypos
!>- ground pressure
!>- nt (thermal profile dimension)
!>- nt lines of
!>- zt,tempC,qw(kg/kg),Ndrops(1/cm3)
!>- nd (thermal profile dimension)
!>- nd lines of
!>- zd,u,v,k,eps
!>
!> WARNINGS:
!> Beware that all dimensions nt,nd must be the same as the first one.
!>
!> Beware that blank lines or lines starting with blanks+ comment_char are
!> NOT ignored.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   meteo_file      "meteo" file name
!-------------------------------------------------------------------------------
subroutine read_meteo_file(meteo_file)
!
character(line_len) :: meteo_file
integer unilog
parameter (unilog = 10)
character(line_len) :: current_line
integer l_iostat
integer first, last
!
integer nt, nd, ns
!
integer count_lines
integer section_length
integer section_count
integer file_count
data file_count /0/
save file_count
integer i

! ----------------------------------------------------------------
! executables
! ----------------------------------------------------------------
file_count = file_count + 1

open(unilog,file=meteo_file,status='old',form='formatted',iostat=l_iostat)
count_lines = 0

! --------------------------------------------------------------
! first loop on file for counting purposes
! --------------------------------------------------------------
init_loop: do while (.true.)
  call find_next_line(unilog,current_line,meteo_file,l_iostat)
  if(l_iostat.lt.0)exit
  count_lines = count_lines + 1
  call bounds(current_line,len(current_line),first,last)
  if (count_lines.eq.4) then
    ! grab the first nt of the current file
    read(current_line,*,iostat=l_iostat) nt
    if(l_iostat.ne.0)then
      write(nfecra,*)"unexpected read error (2) on line: ",current_line
      write(nfecra,*)"found in file : ",meteo_file
      write(nfecra,*)"check and correct your data files"
      write(nfecra,*)"all calculations will be stopped"
      stop
    endif
    if (thermal_profile_dim.lt.0) then
      thermal_profile_dim = nt
    endif
    if (thermal_profile_dim.ne.nt)then
      write(nfecra,*) &
            "all files should provide the same nt for thermal profiles"
      write(nfecra,*) "nt read on file >",meteo_file,"< is ",nt
      write(nfecra,*) "global nt read on first file is ",thermal_profile_dim
      stop
    endif
  endif
  if (count_lines.ge.4) then
    if (count_lines.eq.4+nt+1) then
      read(current_line,*,iostat=l_iostat) nd
      if(l_iostat.ne.0) then
        write(nfecra,*) "while reading nd (dynamical profiles dim)"
        write(nfecra,*) "unexpected read error (3) on line: ",current_line
        write(nfecra,*) "found in file : ",meteo_file
        stop
      endif
      if (dynamical_profile_dim.lt.0) then
        dynamical_profile_dim = nd
      endif
      if (dynamical_profile_dim.ne.nd) then
        write(nfecra,*) &
              "all files should provide the same nd for dynamic profiles"
        call bounds(meteo_file,len(meteo_file),first,last)
        write(nfecra,*) "nd read on file >",meteo_file(first:last),"< is ",nd
        write(nfecra,*) "global nd read on first file is ",dynamical_profile_dim
        stop
      endif
    endif
  endif
enddo init_loop
if (modulo(count_lines,nt+nd+5).ne.0) then
  call bounds(meteo_file,len(meteo_file),first,last)
  write(nfecra,*) "read_meteo_file encountered an error on file >", &
       meteo_file(first:last),"<"
  write(nfecra,*) &
       "all the sections date,pos,pressure+profiles should have the same length"
  write(nfecra,*) &
       "this length should be multiple of ",                        &
       "'5+nt(thermal dim)+nd(dynamical dim)':",                    &
       5+nt+nd
  write(nfecra,*) "but one found ",count_lines," lines of data"
  write(nfecra,*) &
       "probable error cause : ",                                   &
       "the actual and expected length of the profiles differ"
  write(nfecra,*) "all thermal profiles should have ",nt," lines"
  write(nfecra,*) "all dynamical profiles should have ",nd," lines"
  stop
else
  ns = count_lines/(nt+nd+5)
  if (sections_per_file.le.0) then
    sections_per_file = ns
  endif
  if (ns.ne.sections_per_file) then
    call bounds(meteo_file,len(meteo_file),first,last)
    write(nfecra,*)"read_meteo_file encountered an error on file >",&
         meteo_file(first:last),"<"
    write(nfecra,*) &
         " all the files  should contain the same number ",         &
         "of sections date,pos,pressure+profiles"
    write(nfecra,*)"number of sections of current file is ",ns
    write(nfecra,*)"number of sections given by first file is ",    &
         sections_per_file
    stop
  endif
endif

! here the end of file is reached and nt,nd are consistent with global values
! as well as the number of sections in the file
! --------------------------------------------------------------
! allocation of arrays for reading data (only the first time)
! --------------------------------------------------------------
if (file_count.eq.1) then
  call allocate_all()
endif

! -------------------------------------------------------------------
! second loop on file : one can safely read data to load the profiles
! -------------------------------------------------------------------
call bounds(meteo_file,len(meteo_file),first,last)
if (imbrication_verbose) write(nfecra,*) "reopening the file'", &
     meteo_file(first:last),"'"
open(unilog,file=meteo_file,status='old',form='formatted',iostat=l_iostat)
if (l_iostat.ne.0) then
  call bounds(meteo_file,len(meteo_file),first,last)
  write(nfecra,*) "read_meteo_file could not open file '", &
       meteo_file(first:last),"'"
  stop
endif

section_length = nt + nd + 5
section_count = 1
count_lines = 0
read_loop: do while(.true.)
  if (imbrication_verbose) write(nfecra,*) "section count=",section_count
  if (imbrication_verbose) write(nfecra,*) "file count =",file_count
  call find_next_line(unilog,current_line,meteo_file,l_iostat)
  if(l_iostat.lt.0)exit
  count_lines = count_lines + 1
  if (modulo(count_lines,section_length).eq.1) then
    if (imbrication_verbose) write(nfecra,*) &
         "reading years,ord,hour,minute,seconds in:",current_line
    read(current_line,*,iostat=l_iostat)     &
         years(section_count,file_count),   &
         ordinals(section_count,file_count),&
         hours(section_count,file_count),   &
         minutes(section_count,file_count), &
         seconds(section_count,file_count)
    if(l_iostat.ne.0)then
      write(nfecra,*) "while reading years,ord,hour,minute,seconds"
      write(nfecra,*) "unexpected read error (4) on line: ",current_line
      write(nfecra,*) "found in file : ",meteo_file
      stop
    endif
  endif
  if (modulo(count_lines,section_length).eq.2) then
    call bounds(current_line,len(current_line),first,last)
    if (imbrication_verbose) write(nfecra,*) "reading xpos ypos in:", &
         current_line(1:last)
    read(current_line,*,iostat=l_iostat)     &
         xpos(section_count,file_count),    &
         ypos(section_count,file_count)
    if (l_iostat.ne.0) then
      write(nfecra,*) "while reading xpos ypos "
      write(nfecra,*) "unexpected read error (5) on line: ",current_line(1:last)
      write(nfecra,*) "found in file : ",meteo_file
      stop
    endif
  endif
  if (modulo(count_lines,section_length).eq.3) then
    call bounds(current_line,len(current_line),first,last)
    if (imbrication_verbose) write(nfecra,*)"reading ground pressure in:", &
         current_line(1:last)
    read(current_line,*,iostat=l_iostat)     &
         ground_pressure(section_count,file_count)
    if (l_iostat.ne.0) then
      write(nfecra,*) "while reading ground pressure "
      write(nfecra,*) "unexpected read error (6) on line: ", &
           current_line(1:last)
      write(nfecra,*) "found in file : ",meteo_file
      stop
    endif
  endif
  do i = 1, nt, 1
    if (modulo(count_lines,section_length).eq.4+i) then
      call bounds(current_line,len(current_line),first,last)
      if (ippmod(iatmos).eq.2) then
        if (imbrication_verbose) write(nfecra,*) &
             "reading zt,tempC,qw,Nc in:",current_line(1:last)
        read(current_line,*,iostat=l_iostat)    &
             zt(i,section_count,file_count),   &
             tempC(i,section_count,file_count),&
             qw(i,section_count,file_count),   &
             Nc(i,section_count,file_count)
        if(l_iostat.ne.0)then
          write(nfecra,*) "while reading zt,tempC,qw,Nc "
          write(nfecra,*) "unexpected read error (7) on line: ", &
               current_line(1:last)
          write(nfecra,*) "found in file : ",meteo_file
          stop
        endif
      elseif (ippmod(iatmos).le.1) then
        if (imbrication_verbose) write(nfecra,*) &
             "reading zt,tempC,qw in:",current_line(1:last)
        read(current_line,*,iostat=l_iostat)    &
             zt(i,section_count,file_count),   &
             tempC(i,section_count,file_count),&
             qw(i,section_count,file_count)
        if (l_iostat.ne.0) then
          write(nfecra,*) "while reading zt,tempC,qw "
          write(nfecra,*) "unexpected read error (8) on line: ", &
               current_line(1:last)
          write(nfecra,*) "found in file : ",meteo_file
          stop
        endif
      endif
    endif
  enddo
  do i = 1, nd, 1
    if (modulo(count_lines,section_length) &
         .eq.modulo(5+nt+i,section_length)) then
      call bounds(current_line,len(current_line),first,last)
      if (imbrication_verbose) write(nfecra,*) &
           "reading u,v,tke,eps in:",current_line(1:last)
      read(current_line,*,iostat=l_iostat)       &
           zd(i,section_count,file_count),     &
           u(i,section_count,file_count),      &
           v(i,section_count,file_count),      &
           tke(i,section_count,file_count),    &
           eps(i,section_count,file_count)
      if (l_iostat.ne.0) then
        write(nfecra,*) "while reading u,v,tke,eps "
        write(nfecra,*) "unexpected read error (9) on line: ", &
             current_line(1:last)
        write(nfecra,*) "found in file : ",meteo_file
        stop
      endif
    endif
  enddo
  if (modulo(count_lines,section_length).eq.0) &
       section_count = section_count + 1
enddo read_loop
if(imbrication_verbose) write(nfecra,*) "read_loop is finished"
close(unilog)
end subroutine read_meteo_file


! *****************************************************************************
!> \brief    Checks the time variables to ensure the chronology
! *****************************************************************************
subroutine check_chronologies
use atincl, sim_year=>syear
use atincl, sim_ord=>squant
use atincl, sim_hour=>shour
use atincl, sim_min=>smin
use atincl, sim_sec=>ssec
implicit none
logical err_chronologies
integer i, j
double precision sim_time

! ------------------------------------------------
! some checking on chronologies must be done
! they must be synchronized
do i = 2, number_of_files, 1
  do j = 1, sections_per_file, 1
    err_chronologies = .false.
    if(years(j,i).ne.years(j,1)) err_chronologies = .true.
    if(ordinals(j,i).ne.ordinals(j,1)) err_chronologies = .true.
    if(hours(j,i).ne.hours(j,1)) err_chronologies = .true.
    if(minutes(j,i).ne.minutes(j,1)) err_chronologies = .true.
    if(seconds(j,i).ne.seconds(j,1)) err_chronologies = .true.
    if (err_chronologies) then
      write(nfecra,*) &
           "the chronologies of the different profiles are not synchronized"
      write(nfecra,*) &
           "faulty file:",imbrication_files(i)
      write(nfecra,*) &
           "faulty date:",years(j,i),ordinals(j,i), &
           hours(j,i),minutes(j,i),seconds(j,i)
      write(nfecra,*) &
           "should be equal to date:",years(j,1),   &
           ordinals(j,1),hours(j,1),minutes(j,1),seconds(j,1)
      write(nfecra,*) &
           "defined in file:",imbrication_files(1)
      write(nfecra,*) &
           "section:",j
      stop
   endif
 enddo
enddo

! they must be in the natural temporal order
allocate(times(sections_per_file,number_of_files))
if (sim_year.lt.0) then
  sim_year = years(1,1)
  sim_ord = ordinals(1,1)
  sim_hour = hours(1,1)
  sim_min = minutes(1,1)
  sim_sec = seconds(1,1)
  sim_time = yo2j(years(1,1),ordinals(1,1))*86400.d+00+ &
       hours(1,1)*3600.+ &
       minutes(1,1)*60.+ &
       seconds(1,1)
else
  sim_time=yo2j(sim_year,sim_ord)*86400.d+00+ &
       sim_hour*3600.+ &
       sim_min*60.+ &
       sim_sec
endif

do i = 1, number_of_files, 1
  do j = 1, sections_per_file, 1
    times(j,i) = yo2j(years(j,i),ordinals(j,i))*86400.d+00+ &
         hours(j,i)*3600.+ &
         minutes(j,i)*60.+ &
         seconds(j,i)
  enddo
enddo

do i = 1, number_of_files, 1
  do j = 1, sections_per_file, 1
    times(j,i) = times(j,i) - sim_time
    if (imbrication_verbose) write(nfecra,*) &
         "simulation times:",times(j,i)
  enddo
enddo

do i = 1, number_of_files, 1
  do j = 2, sections_per_file, 1
    if (times(j,i).le.times(1,i)) then
      write(nfecra,*) &
           "the chronologies of the different profiles are not in order"
      write(nfecra,*) &
           "faulty file:",imbrication_files(i)
      write(nfecra,*) &
           "faulty date:",years(j,i),ordinals(j,i), &
           hours(j,i),minutes(j,i),seconds(j,i)
      write(nfecra,*) &
           "defined in section ",j
      write(nfecra,*) &
           "should be posterior to date:",years(1,i), &
           ordinals(1,i),hours(1,i),minutes(1,i),seconds(1,i)
      write(nfecra,*) &
           "defined in section ",1
      stop
    endif
  enddo
enddo

do i = 1, number_of_files, 1
  do j = 1, sections_per_file,1
    if(imbrication_verbose) write(nfecra,*) &
         "simulation times:",times(j,i)
  enddo
enddo
end subroutine check_chronologies

! -----------------------------------------------------------------------------
!> \brief Check that the profiles position is the same over time
! -----------------------------------------------------------------------------
subroutine check_positions
implicit none
integer i,j,k
do i = 1, number_of_files, 1
  do j = 2, sections_per_file, 1
    if (xpos(j,i).ne.xpos(1,i)) then
      write(nfecra,*) &
           "the x-positions of the profiles in file ",imbrication_files(i)
      write(nfecra,*) &
           "are not consistent (vary in time)"
      write(nfecra,*) &
           "faulty section is :",j
      write(nfecra,*) &
           " xpos(1)=",xpos(1,i)
      write(nfecra,*) &
           " xpos(",j,")=",xpos(j,i)
      stop
    endif
    if (ypos(j,i).ne.ypos(1,i))then
      write(nfecra,*) &
           "the y-positions of the profiles in file ",imbrication_files(i)
      write(nfecra,*) &
           "are not consistent: they vary in time"
      write(nfecra,*) &
           "the faulty section is :",j
      write(nfecra,*) &
           " ypos(1)=",ypos(1,i)
      write(nfecra,*) &
           " ypos(",j,")=",ypos(j,i)
      stop
    endif
  enddo
enddo

do i = 1, number_of_files, 1
  do k = 1, number_of_files, 1
    if (k.ne.i) then
      if (xpos(1,i).eq.xpos(1,k) .and. ypos(1,i).eq.ypos(1,k)) then
        write(nfecra,*) &
             "the positions given of some profiles are not consistent"
        write(nfecra,*) &
             "The positions of the profiles in file ",imbrication_files(i)
        write(nfecra,*) &
             "and the positions of the profiles in file ",imbrication_files(k)
        write(nfecra,*) &
             "are equal."
        stop
      endif
    endif
  enddo
enddo
end subroutine check_positions


! -----------------------------------------------------------------------------
!> \brief Check that the profiles vertical grids heights
!>      are strictly increasing
! -----------------------------------------------------------------------------
subroutine check_altitudes
implicit none
integer i,j,k
integer kk
integer first,last

do i = 1, number_of_files, 1
  do j = 1, sections_per_file, 1
    do k = 2, thermal_profile_dim
      if (zt(k-1,j,i).ge.zt(k,j,i)) then
        write(nfecra,*) "the thermal profile in section ",j
        call bounds(imbrication_files(i),len(imbrication_files(i)),first,last)
        write(nfecra,*)"of the file '",imbrication_files(i)(first:last),"'"
        write(nfecra,*)"is not strictly increasing"
        write(nfecra,*)"erroneous level ",k," with zt =",zt(k,j,i)
        do kk = 1, thermal_profile_dim, 1
          write(nfecra,*)"k=",kk,"zt=",zt(kk,j,i)
        enddo
        stop
      endif
    enddo
  enddo
enddo

do i = 1, number_of_files, 1
  do j = 1, sections_per_file, 1
    do k = 2, dynamical_profile_dim
      if (zd(k-1,j,i).ge.zd(k,j,i)) then
        write(nfecra,*) "the dynamical profile in section ", j
        call bounds(imbrication_files(i),len(imbrication_files(i)),first,last)
        write(nfecra,*) "of the file '",imbrication_files(i)(first:last),"'"
        write(nfecra,*) "is not strictly increasing"
        write(nfecra,*) "erroneous level ",k," with zd =",zd(k,j,i)
        do kk = 1, dynamical_profile_dim, 1
          write(nfecra,*) "k=",kk,"zd=",zd(kk,j,i)
        enddo
        stop
      endif
    enddo
  enddo
enddo
end subroutine check_altitudes


! -----------------------------------------------------------------------------
!> \brief Compute the hydrostastic pressure by Laplace integration
! -----------------------------------------------------------------------------
subroutine hydrostatic_pressure

use cstphy, only: tkelvi
use cstphy, only: rair
use cstphy, only: gz
use atincl

implicit none
integer i,j,k
double precision tmoy
double precision q0,q1
double precision rho_moy
integer ih2o
data ih2o/0/
save ih2o
double precision rap

if (ippmod(iatmos).eq.2) ih2o=1
if(.not.allocated(pressure))then
  allocate(pressure(thermal_profile_dim,sections_per_file,number_of_files))
endif
integration_direction:if (ihpm.eq.0) then
  ! bottom to top integration
  do i = 1, number_of_files, 1
    if (imbrication_verbose) write(nfecra,*) &
         "hydrostatic_pressure::file:",imbrication_files(i)
    do j = 1, sections_per_file, 1
      if (imbrication_verbose) write(nfecra,*) &
           "hydrostatic_pressure::section:",j
      if (imbrication_verbose) write(nfecra,*) &
           "hydrostatic_pressure::thermal_profile_dim=",thermal_profile_dim
      pressure(1,j,i) = ground_pressure(j,i)
      do k = 2, thermal_profile_dim, 1
        if (imbrication_verbose) write(nfecra,*) &
             "hydrostatic_pressure::k=",k
        tmoy = 0.5d0*(tempC(k-1,j,i) + tempC(k,j,i)) + tkelvi
        if (ippmod(iatmos).eq.2) then
          q0 = min( qw(k-1,j,i), cs_air_yw_sat( tempC(k-1,j,i), &
               pressure(k-1,j,i)))
          q1 = min( qw(k  ,j,i), cs_air_yw_sat( tempC(k  ,j,i), &
               pressure(k-1,j,i)))
        else
          q0 = qw(k-1,j,i)
          q1 = qw(k  ,j,i)
        endif
        rho_moy = rair*(1.d+00+(rvsra-1.)*(q0 + q1)/2.d0*ih2o)
        rap = -abs(gz)*(zt(k,j,i)-zt(k-1,j,i))/rho_moy/tmoy
        pressure(k,j,i)=pressure(k-1,j,i)*exp(rap)
      enddo
    enddo
  enddo
else
  !top to bottom integration
  do i = 1, number_of_files, 1
    do j = 1, sections_per_file, 1
      pressure(thermal_profile_dim,j,i) = &
           101325d+00*(288.15/(288.15-6.5d-3*zt(thermal_profile_dim,j,i)))**( &
           -abs(gz)/rair/6.5d-3)
      do k = thermal_profile_dim, 2, -1
        tmoy = 0.5d0*(tempC(k-1,j,i) + tempC(k,j,i)) + tkelvi
        if (ippmod(iatmos).eq.2) then
          q0 = min( qw(k-1,j,i), cs_air_yw_sat( tempC(k-1,j,i), &
               pressure(k  ,j,i)))
          q1 = min( qw(k  ,j,i), cs_air_yw_sat( tempC(k  ,j,i), &
               pressure(k  ,j,i)))
        else
          q0 = qw(k-1,j,i)
          q1 = qw(k  ,j,i)
        endif
        rho_moy = rair*(1.d+00+(rvsra-1.)*(q0 + q1)/2.d0*ih2o)
        rap = abs(gz)*(zt(k,j,i)-zt(k-1,j,i))/rho_moy/tmoy
        pressure(k-1,j,i) = pressure(k  ,j,i)*exp(rap)
      enddo
    enddo
  enddo
endif integration_direction

do i = 1,number_of_files, 1
  if (imbrication_verbose) write(nfecra,*) &
       "hydrostatic_pressure::file:",imbrication_files(i)
  do j = 1,sections_per_file, 1
    if (imbrication_verbose) write(nfecra,*) &
         "hydrostatic_pressure::section:",j
    if (imbrication_verbose) write(nfecra,*)  &
         "hydrostatic_pressure::date:",years(j,i),ordinals(j,i), &
         hours(j,i),minutes(j,i),seconds(j,i)
    do k = 1, thermal_profile_dim, 1
      if(imbrication_verbose) write(nfecra,*) &
           "hydrostatic_pressure::z,t,p:",zt(k,j,i),tempC(k,j,i),&
           pressure(k,j,i)
    enddo
  enddo
enddo
end subroutine hydrostatic_pressure

! -----------------------------------------------------------------------------
!> \brief Computes the potential_temperature_and_density profiles
! -----------------------------------------------------------------------------
subroutine potential_temperature_and_density
use cstphy, only: rair
use cstphy, only: tkelvi
use cstphy, only: cp0
use atincl
implicit none
double precision rhum! I prefer bourbon
integer ih2o
data ih2o/0/
save ih2o
double precision rscp
integer i,j,k

if (.not.allocated(theta)) then
  allocate(theta(thermal_profile_dim,sections_per_file,number_of_files))
endif
if (.not.allocated(density)) then
  allocate(density(thermal_profile_dim,sections_per_file,number_of_files))
endif

if (ippmod(iatmos).eq.2) ih2o = 1
do i = 1, number_of_files, 1
  do j = 1,sections_per_file, 1
    do k = 1,thermal_profile_dim, 1
      rhum = rair*(1.d0+(rvsra-1.d0)*qw(k,j,i)*ih2o)
      if (ippmod(iatmos).eq.0) then
        ! constant density
        density(k,j,i) = pressure(1,j,i)/(tempC(k,j,i) + tkelvi)/rhum
      else
        ! variable density
        density(k,j,i) = pressure(k,j,i)/(tempC(k,j,i) + tkelvi)/rhum
      endif
      rscp = (rair/cp0)*(1.d0 + (rvsra-cpvcpa)*qw(k,j,i)*ih2o)
      theta(k,j,i) = (tempC(k,j,i)+tkelvi)*((ps/pressure(k,j,i))**rscp)
    enddo
  enddo
enddo
do i = 1, number_of_files, 1
   if (imbrication_verbose) &
   write(nfecra,*) &
   "potential_temperature_and_density::file:",imbrication_files(i)
   do j = 1,sections_per_file, 1
     if (imbrication_verbose) &
          write(nfecra,*) "potential_temperature_and_density::section:",j
     if (imbrication_verbose) &
          write(nfecra,*)     &
          "potential_temperature_and_density::date:",years(j,i), &
          ordinals(j,i),hours(j,i),minutes(j,i),seconds(j,i)
     do k = 1, thermal_profile_dim, 1
       if (imbrication_verbose) &
            write(nfecra,*) &
            "z,t,p,potential_temperature,density:::",zt(k,j,i), &
            tempC(k,j,i),pressure(k,j,i),theta(k,j,i),density(k,j,i)
     enddo
   enddo
 enddo
 end subroutine potential_temperature_and_density

! -----------------------------------------------------------------------------
!> \brief Search for the position of a value in an array,
!> assuming that the array is sorted in a strictly increasing order
!>
!> return if possible lower,upper such that :
!>-                    the_array(lower) <= the_value <= the_array(upper)
!>        otherwise :
!>          lower==upper if the_value<the_array(first)
!>                       or if the_value>the_array(last)
!>          lower> upper if none of the previous cases applies (anomaly)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]       the_array   an array
!> \param[in]       the_value   a particular value
!> \param[in]       lower       the index of the first membre of the array
!>                                  lower than the value
!> \param[in]       upper       the index of the first membre of the array
!>                                  greater than the value
!-------------------------------------------------------------------------------
subroutine get_index(the_array,the_value,lower,upper)
implicit none

double precision,dimension(:),intent(in) :: the_array
double precision,intent(in) :: the_value
integer,intent(out) ::upper,lower
integer dmin,dmax
integer i

dmin = lbound(the_array,1)
dmax = ubound(the_array,1)
do i = dmin, dmax-1, 1
  if (the_array(i).le.the_value .and. the_value.le.the_array(i+1)) then
    lower = i
    upper = i+1
    return
  endif
enddo

if (the_value.lt.the_array(dmin)) then
  lower = dmin
  upper = dmin
  return
endif

if (the_value.gt.the_array(dmax)) then
  lower = dmax
  upper = dmax
  return
endif
lower = dmax
upper = dmin
end subroutine get_index

! -----------------------------------------------------------------------------
!> \brief  Interpolates a "profile" at a given time.
!> Given a series of profiles varying in time
!> you get the profile interpolated from them at the given time.
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   the_time                current time
!> \param[in]   the_times               times array
!> \param[in]   the_profiles            input profiles
!> \param[out]   interpolated_profile    output profile
!-------------------------------------------------------------------------------
subroutine time_interpolation(the_time,the_times, &
     the_profiles,interpolated_profile)
implicit none
double precision,intent(in) :: the_time
double precision, dimension(:),intent(in) :: the_times
double precision, dimension(:,:),intent(in) :: the_profiles
double precision, dimension(:),intent(out) :: interpolated_profile
integer lower,upper
integer dmin,dmax
double precision weight
integer i

call get_index(the_times,the_time,lower,upper)
if (lower.lt.upper) then
  weight = (the_time-the_times(lower))/(the_times(upper)-the_times(lower))
  if (imbrication_verbose) write(nfecra,*) &
       "time_interpolation:: weight=",weight
  interpolated_profile = the_profiles(:,lower)*(1.-weight)  &
       + the_profiles(:,upper)*weight
elseif (lower.eq.upper) then
  interpolated_profile = the_profiles(:,lower)
else
  write(nfecra,*) &
       "time_interpolation:: the times array is not increasing"
  dmin = lbound(the_times,1)
  dmax = ubound(the_times,1)
  do i = dmin, dmax, 1
    write(nfecra,*) "time_interpolation:: the_times(",i,")=",the_times(i)
  enddo
  write(nfecra,*) "time_interpolation stops the calculations"
  stop
endif
end subroutine time_interpolation


!-------------------------------------------------------------------------------
!> \brief interpolates in a profile at a given altitude
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   the_altitude        current height
!> \param[in]   the_altitudes       height array
!> \param[in]   the_profile         the profile
!> \param[out]   interpolated_value  interpolated profile
!-------------------------------------------------------------------------------
subroutine altitude_interpolation(the_altitude,the_altitudes,the_profile, &
     interpolated_value)
implicit none
double precision,intent(in) :: the_altitude
double precision, dimension(:),intent(in) :: the_altitudes
double precision, dimension(:),intent(in) :: the_profile
double precision, intent(out) :: interpolated_value
integer lower,upper
integer dmin,dmax
double precision weight
integer i

call get_index(the_altitudes,the_altitude,lower,upper)
if (lower.lt.upper) then
  weight = (the_altitude-the_altitudes(lower)) &
       / (the_altitudes(upper)-the_altitudes(lower))
  if (imbrication_verbose) write(nfecra,*)     &
       "altitude_interpolation:: weight=", weight
  interpolated_value = the_profile(lower)*(1.-weight) &
       + the_profile(upper)*weight
elseif (lower.eq.upper) then
  interpolated_value = the_profile(lower)
else
  write(nfecra,*) &
       "altitude_interpolation:: the altitudes array is not increasing"
  dmin = lbound(the_altitudes,1)
  dmax = ubound(the_altitudes,1)
  do i = dmin, dmax, 1
    write(nfecra,*) "altitude_interpolation:: the_altitudes(",i,")=", &
         the_altitudes(i)
  enddo
  write(nfecra,*) "altitude_interpolation stops all the calculations"
  stop
endif
end subroutine altitude_interpolation


! -----------------------------------------------------------------------------
!> \brief  Compute radius of influence
! -----------------------------------------------------------------------------
subroutine red_tape
implicit none
integer i,j
allocate(coordinates_th(3,thermal_profile_dim,number_of_files))
allocate(coordinates_dyn(3,dynamical_profile_dim,number_of_files))
allocate(influence_param_th(3,thermal_profile_dim,number_of_files))
do i = 1, number_of_files, 1
  do j = 1, thermal_profile_dim, 1
    influence_param_th(1,j,i) = 1.d+00/horizontal_influence_radius
    influence_param_th(2,j,i) = 1.d+00/horizontal_influence_radius
    influence_param_th(3,j,i) = 1.d+00/vertical_influence_radius
  enddo
enddo
allocate(influence_param_dyn(3,dynamical_profile_dim,number_of_files))
do i = 1, number_of_files, 1
  do j = 1, dynamical_profile_dim, 1
    influence_param_dyn(1,j,i) = 1.d+00/horizontal_influence_radius
    influence_param_dyn(2,j,i) = 1.d+00/horizontal_influence_radius
    influence_param_dyn(3,j,i) = 1.d+00/vertical_influence_radius
  enddo
enddo
end subroutine red_tape

! -----------------------------------------------------------------------------
!> \brief  Identification of the first and last non white character
!>      of a string
! -----------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   string      the input string
!> \param[in]   length      its length
!> \param[out]  b           number of the first non white character
!> \param[out]  e           number of the last non white character
!-------------------------------------------------------------------------------
subroutine bounds(string,length,b,e)
!
! Being compelled to write such low level stuff because
! the designers of the language didn't bother to offer
! a decent suite of string manipulating tools really sucks.
!
implicit none
integer length
character*(length) string
integer b, e
integer i
b = 1
e = length
do i = 1, length, 1
  if ( string(i:i).eq." ") then
    b = b + 1
  else
    goto 1 ! simulate a "break"
  endif
enddo
1 continue
do i = length,1,-1
  if ( string(i:i).eq." ") then
    e = e-1
  else
    goto 2 ! simulate a "break"
  endif
enddo
2 continue
return
end subroutine bounds


! -----------------------------------------------------------------------------
!> \brief Prepare data for imbrication by reading meteo files
!>
!> Warning : the list of files is supposed to be "imbrication_files_list.txt"
!>
! -----------------------------------------------------------------------------
subroutine activate_imbrication
implicit none

integer i,j,k
integer first,last

write(nfecra,*)"*******************************"
write(nfecra,*)"Atmospheric Imbrication:       "
write(nfecra,*)"*******************************"

imbrication_files_list = "imbrication_files_list.txt"
call read_files_list(imbrication_files_list,imbrication_files)
write(nfecra,*)"number_of_files            : ",number_of_files
do i = 1, number_of_files
  call bounds(imbrication_files(i),len(imbrication_files(i)),first,last)
  if (imbrication_verbose) write(nfecra,*) &
       "file number ",i,"=","'",imbrication_files(i)(first:last),"'"
  call read_meteo_file(imbrication_files(i))
enddo

if (imbrication_verbose) then
  loop_on_files: do i = 1, number_of_files, 1
    call bounds(imbrication_files(i),len(imbrication_files(i)),first,last)
    write(nfecra,*) "file number ",i,"=","'", &
         imbrication_files(i)(first:last),"'"
    write(nfecra,*) "number of sections per file: ",sections_per_file
    loop_on_sections: do j = 1, sections_per_file, 1
      write(nfecra,*)"date:",years(j,i),ordinals(j,i), &
           hours(j,i),minutes(j,i),seconds(j,i)
      write(nfecra,*)"xpos,ypos:",xpos(j,i),ypos(j,i)
      write(nfecra,*)"ground_pressure:",ground_pressure(j,i)
      write(nfecra,*)"thermal profiles dim: ",thermal_profile_dim
      loop_on_thermal_profile: do k = 1, thermal_profile_dim, 1
        if (ippmod(iatmos).eq.2) then
          write (nfecra,*) "z,temp,qw,nc=", zt(k,j,i), &
               tempC(k,j,i),qw(k,j,i),Nc(k,j,i)
        elseif (ippmod(iatmos).eq.1) then
          write(nfecra,*) "z,temp,qw=", zt(k,j,i),     &
               tempC(k,j,i),qw(k,j,i)
        else
          write(nfecra,*) "z,temp=",zt(k,j,i),tempC(k,j,i)
        endif
      enddo loop_on_thermal_profile
      write(nfecra,*)"     dynamical profiles dim: ",dynamical_profile_dim
      loop_on_dynamical_profile: do k=1,dynamical_profile_dim,1
        write(nfecra,*)"z,u,v,k,eps=",zd(k,j,i),u(k,j,i),v(k,j,i), &
             tke(k,j,i),eps(k,j,i)
      enddo loop_on_dynamical_profile
    enddo loop_on_sections
  enddo loop_on_files
endif !imbrication_verbose

! ------------------------------------------------
! some checking on chronologies must be done
! ------------------------------------------------
call check_chronologies

! ------------------------------------------------
! some checking on positions of the profiles
! ------------------------------------------------
call check_positions
call check_altitudes

! ------------------------------------------------
! reading terminated: some calculations are done
! ------------------------------------------------
! calculating pressure by hydrostatic (Laplace) integration
! ------------------------------------------------
call hydrostatic_pressure

! ------------------------------------------------
! calculating potential temperature and density
! ------------------------------------------------
call potential_temperature_and_density

! at this point all meteorological profiles are:
! (1) with the same dimensions
! (2) synchronized
! (3) completed up to 11000 if iatra1==1
! one can interpolate in time to get values
!
times_sequence=>times(1:sections_per_file,1)
!
end subroutine activate_imbrication


! --------------------------------------------------------------
!> \brief Prepare for the cressman interpolation of the variables
! -----------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   the_time        current time
!-------------------------------------------------------------------------------
subroutine summon_cressman(the_time)
implicit none
double precision the_time
logical first_call
data first_call /.true./
save first_call
integer lb1, ub1
integer lb2, ub2
integer nbmes
integer i,j
integer, dimension(:,:), allocatable :: ones
if (first_call) then
  if (cressman_u) then
    call mestcr("u",len("u"),1,0,id_u)
  endif
  if (cressman_v) then
    call mestcr("v",len("v"),1,0,id_v)
  endif
  if (cressman_tke) then
    call mestcr("tke",len("tke"),1,0,id_tke)
  endif
  if (cressman_eps) then
    call mestcr("eps",len("eps"),1,0,id_eps)
  endif
  if (cressman_theta .and. ippmod(iatmos).ge.1) then
    call mestcr("theta",len("theta"),1,0,id_theta)
  endif
  if (cressman_qw .and. ippmod(iatmos).ge.2) then
    call mestcr("qw",len("qw"),1,0,id_qw)
  endif
  if (cressman_nc.and.ippmod(iatmos).ge.2) then
    call mestcr("nc",len("nc"),1,0,id_nc)
  endif
  call red_tape
  first_call = .false.
endif

! --------------------------------------------------------------
!
! --------------------------------------------------------------

call interpolate_all_profiles(the_time)

if (cressman_u) then
  lb1 = lbound(ti_u,1)
  lb2 = lbound(ti_u,2)
  ub1 = ubound(ti_u,1)
  ub2 = ubound(ti_u,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if(imbrication_verbose) &
       write(nfecra,*)"nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files &
       .or.lb1.ne.1 .or.ub1.ne.dynamical_profile_dim) then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated u are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:", &
         dynamical_profile_dim,",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      coordinates_dyn(1,j,i) = xpos(1,i)
      coordinates_dyn(2,j,i) = ypos(1,i)
      coordinates_dyn(3,j,i) = ti_zd(j,i)
      if(imbrication_verbose) &
           write(nfecra,*)"j,i=",j,i,"x,y,z=",coordinates_dyn(1,j,i), &
           coordinates_dyn(2,j,i), &
           coordinates_dyn(3,j,i), &
           "u=",ti_u(j,i)
    enddo
  enddo
  allocate(ones(dynamical_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_u,nbmes,ti_u,coordinates_dyn,ones,ones,influence_param_dyn)
  deallocate(ones)
endif

if (cressman_v) then
  lb1 = lbound(ti_v,1)
  lb2 = lbound(ti_v,2)
  ub1 = ubound(ti_v,1)
  ub2 = ubound(ti_v,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if (imbrication_verbose) &
       write(nfecra,*) "nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files &
       .or.lb1.ne.1 .or.ub1.ne.dynamical_profile_dim)then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated v are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:",dynamical_profile_dim, &
         ",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      coordinates_dyn(1,j,i) = xpos(1,i)
      coordinates_dyn(2,j,i) = ypos(1,i)
      coordinates_dyn(3,j,i) = ti_zd(j,i)
      if (imbrication_verbose) &
           write(nfecra,*)"j,i=",j,i,"x,y,z=",coordinates_dyn(1,j,i),&
           coordinates_dyn(2,j,i),&
           coordinates_dyn(3,j,i),&
           "v=",ti_v(j,i)
    enddo
  enddo
  allocate(ones(dynamical_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_v,nbmes,ti_v,coordinates_dyn,ones,ones,influence_param_dyn)
  deallocate(ones)
endif

if (cressman_tke) then
  lb1 = lbound(ti_tke,1)
  lb2 = lbound(ti_tke,2)
  ub1 = ubound(ti_tke,1)
  ub2 = ubound(ti_tke,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if (imbrication_verbose) &
       write(nfecra,*)"nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files &
       .or.lb1.ne.1 .or.ub1.ne.dynamical_profile_dim)then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated tke are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:",dynamical_profile_dim, &
         ",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      coordinates_dyn(1,j,i) = xpos(1,i)
      coordinates_dyn(2,j,i) = ypos(1,i)
      coordinates_dyn(3,j,i) = ti_zd(j,i)
      if (imbrication_verbose) &
           write(nfecra,*) &
           "j,i=",j,i,"x,y,z=",coordinates_dyn(1,j,i),&
           coordinates_dyn(2,j,i),&
           coordinates_dyn(3,j,i),&
           "tke=",ti_tke(j,i)
    enddo
  enddo
  allocate(ones(dynamical_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_tke,nbmes,ti_tke,coordinates_dyn,ones,ones,influence_param_dyn)
  deallocate(ones)
endif

if (cressman_eps) then
  lb1 = lbound(ti_eps,1)
  lb2 = lbound(ti_eps,2)
  ub1 = ubound(ti_eps,1)
  ub2 = ubound(ti_eps,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if (imbrication_verbose) &
       write(nfecra,*) "nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files &
       .or.lb1.ne.1 .or.ub1.ne.dynamical_profile_dim)then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated eps are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:", &
         dynamical_profile_dim,",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1,dynamical_profile_dim, 1
      coordinates_dyn(1,j,i) = xpos(1,i)
      coordinates_dyn(2,j,i) = ypos(1,i)
      coordinates_dyn(3,j,i) = ti_zd(j,i)
      if (imbrication_verbose) &
           write(nfecra,*)"j,i=",j,i,"x,y,z=",coordinates_dyn(1,j,i),&
           coordinates_dyn(2,j,i),&
           coordinates_dyn(3,j,i),&
           "eps=",ti_eps(j,i)
    enddo
  enddo
  allocate(ones(dynamical_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1, dynamical_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_eps,nbmes,ti_eps,coordinates_dyn,ones,ones,influence_param_dyn)
  deallocate(ones)
endif

if (cressman_theta.and.ippmod(iatmos).ge.1) then
  lb1 = lbound(ti_theta,1)
  lb2 = lbound(ti_theta,2)
  ub1 = ubound(ti_theta,1)
  ub2 = ubound(ti_theta,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if (imbrication_verbose) &
       write(nfecra,*)"nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files &
       .or.lb1.ne.1 .or.ub1.ne.thermal_profile_dim)then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated theta are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:", &
         thermal_profile_dim,",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1,thermal_profile_dim, 1
      coordinates_th(1,j,i) = xpos(1,i)
      coordinates_th(2,j,i) = ypos(1,i)
      coordinates_th(3,j,i) = ti_zt(j,i)
      if (imbrication_verbose) &
           write(nfecra,*)"j,i=",j,i,"x,y,z=",coordinates_th(1,j,i), &
           coordinates_th(2,j,i),&
           coordinates_th(3,j,i),&
           "theta=",ti_theta(j,i)
    enddo
  enddo
  allocate(ones(thermal_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1,thermal_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_theta,nbmes,ti_theta, &
       coordinates_th,ones,ones,influence_param_th)
  deallocate(ones)
endif

if (cressman_qw.and.ippmod(iatmos).ge.1)then
  lb1 = lbound(ti_qw,1)
  lb2 = lbound(ti_qw,2)
  ub1 = ubound(ti_qw,1)
  ub2 = ubound(ti_qw,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if (imbrication_verbose) &
       write(nfecra,*)"nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files  &
       .or.lb1.ne.1 .or.ub1.ne.thermal_profile_dim) then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated qw are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:", &
         thermal_profile_dim,",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1, thermal_profile_dim, 1
      coordinates_th(1,j,i) = xpos(1,i)
      coordinates_th(2,j,i) = ypos(1,i)
      coordinates_th(3,j,i) = ti_zt(j,i)
      if (imbrication_verbose) &
           write(nfecra,*)"j,i=",j,i,"x,y,z=",coordinates_th(1,j,i),&
           coordinates_th(2,j,i),&
           coordinates_th(3,j,i),&
           "qw=",ti_qw(j,i)
    enddo
  enddo
  allocate(ones(thermal_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1, thermal_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_qw,nbmes,ti_qw,coordinates_th,ones,ones,influence_param_th)
  deallocate(ones)
endif

if (cressman_nc.and.ippmod(iatmos).ge.2) then
  lb1 = lbound(ti_nc,1)
  lb2 = lbound(ti_nc,2)
  ub1 = ubound(ti_nc,1)
  ub2 = ubound(ti_nc,2)
  nbmes = (ub1-lb1+1)*(ub2-lb2+1)
  if (imbrication_verbose) &
       write(nfecra,*)"nbmes=",nbmes
  if (lb2.ne.1 .or. ub2.ne.number_of_files  &
       .or.lb1.ne.1 .or.ub1.ne.thermal_profile_dim)then
    write(nfecra,*) &
         "in module imbrication::summon_cressman"
    write(nfecra,*) &
         "dimensions of time interpolated nc are not consistent"
    write(nfecra,*) &
         "expected dimensions are: (1:", &
         thermal_profile_dim,",1:",number_of_files,")"
    write(nfecra,*) &
         "  actual dimensions are: (",lb1,":",ub1,",",lb2,":",ub2,")"
    write(nfecra,*) &
         "all calculations will be stopped"
    stop
  endif
  do i = 1, number_of_files, 1
    do j = 1, thermal_profile_dim, 1
      coordinates_th(1,j,i) = xpos(1,i)
      coordinates_th(2,j,i) = ypos(1,i)
      coordinates_th(3,j,i) = ti_zt(j,i)
      if (imbrication_verbose) &
           write(nfecra,*)"j,i=",j,i,"x,y,z=",coordinates_th(1,j,i),&
           coordinates_th(2,j,i),&
           coordinates_th(3,j,i),&
           "nc=",ti_nc(j,i)
    enddo
  enddo
  allocate(ones(thermal_profile_dim,number_of_files))
  do i = 1, number_of_files, 1
    do j = 1, thermal_profile_dim, 1
      ones(j,i) = 1
    enddo
  enddo
  call mesmap(id_nc,nbmes,ti_nc,coordinates_th,ones,ones,influence_param_th)
  deallocate(ones)
endif

end subroutine summon_cressman

end module atimbr

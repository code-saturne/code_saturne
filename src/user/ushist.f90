!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine ushist &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtpa   , rtp    , propce )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Non-standard monitoring point definition.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)

! Local variables

integer          ii, kk, node, ndrang, nvarpp, numcel, lng
double precision xx, yy, zz, xyztmp(3)

! Monitoring points number (lower than a maximum of 100)
integer          ncapmx
parameter       (ncapmx=100)
integer          icapt(ncapmx)
save             icapt
integer          ircapt(ncapmx)
save             ircapt

! Number of monitoring points
integer          ncapts
save             ncapts

! Current pass number
integer          ipass
data             ipass /0/
save             ipass

! Temporary array
double precision vacapt(ncapmx)

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initialization
!===============================================================================

! Memory management


! Current pass number in this subroutine

ipass = ipass + 1

!===============================================================================
! 2.  Search for the monitoring points
!===============================================================================

! The numbers stored in the 'ircapt' array give the processor rank on which
!   is the probes. The user should not have to care as soon as she/he is
!   using the 'findpt' subroutine to find the monitoring points.

! At the first pass:
!    Search for the cell number the centre of which is the closest of the
!    coordinates (xx, yy, zz).
!    In case of parallelism, the cell number 'icapt(ii)' is local to the
!    processor of rank 'ircapt(ii)' (from 0 to the number of processor - 1).
!    'ncapts' gives the total number of monitoring points.

if (ipass.eq.1) then

  ii = 0

  xx = 0.20d0
  yy = 0.15d0
  zz = 0.01d0
  call findpt &
  !==========
 ( ncelet , ncel   , xyzcen ,                 &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  xx = 0.70d0
  yy = 0.15d0
  zz = 0.01d0
  call findpt &
  !==========
 ( ncelet , ncel   , xyzcen ,                  &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  xx = 0.20d0
  yy = 0.75d0
  zz = 0.01d0
  call findpt &
  !==========
 ( ncelet , ncel   , xyzcen ,                  &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  xx = 0.70d0
  yy = 0.75d0
  zz = 0.01d0
  call findpt                                                     &
  !==========
 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node  , ndrang)
  ii = ii + 1
  icapt(ii) = node
  ircapt(ii) = ndrang

  ncapts = ii

  if(ii.gt.ncapmx) then
    write(nfecra,*) ' ushist: ncapmx should at least be', ii
    call csexit(1)
    !==========
  endif

endif


!===============================================================================
! 3.  Open files: example for a variable per file
!===============================================================================

! Number of variables = number of files

nvarpp = nvar


! At the first pass: open files and write a header

if (ipass.eq.1) then

  ! Test the maximum number of user files

  if (nvarpp.gt.nushmx) then
    write(nfecra,*) &
      ' ushist: no more than ', nushmx,' monitoring files are allowed'
    call csexit(1)
    !==========
  endif

  do ii = 1, nvarpp

    ! Open the files with the availabe Fortran 'file units'

    if (irangp.le.0) then
      open(file=ficush(ii), unit=impush(ii))
    endif

    ! Print the (global) cell number and the center coordinates

    do kk = 1, ncapts

      ! Cell number (in a parallel run: local to the current processor)
      numcel = icapt(kk)

      if (irangp.lt.0 .or. irangp.eq.ircapt(kk)) then
        ! Cell coordinates (in a parallel run: only one processor gives values)
        xyztmp(1) = xyzcen(1,numcel)
        xyztmp(2) = xyzcen(2,numcel)
        xyztmp(3) = xyzcen(3,numcel)
      else
        ! Fake values on the other processors
        xyztmp(1) = 0.d0
        xyztmp(2) = 0.d0
        xyztmp(3) = 0.d0
      endif

      ! In case of parallelism, the processor on which the cell was found
      ! sends its global number and coordinates to the others.
      if (irangp.ge.0) then
        lng = 3
        call parbcr(ircapt(kk), lng, xyztmp)
        !==========
      endif

      ! Write information
      !   (only rank 0 works in a parallel run: only one file is needed)
      if (irangp.le.0) then
        write(impush(ii),1000) &
          '#', ' Coord ', xyztmp(1), xyztmp(2), xyztmp(3)
      endif

    enddo

  enddo

endif

1000 format(a,a9,3e14.5)

!===============================================================================
! 4.  Write values: example for a variable per file
!===============================================================================

! Write the time step number,
!       the physical time value
!       the variable at each monitoring points
! In a serial run:   the value is merely 'rtp(icapt(kk),ii)'
! In a parallel run: the value may come from a different processor, to be
!                    determined in 'vacapt(kk)' with the 'parhis' subroutine)

do ii = 1 , nvarpp
  do kk = 1, ncapts
    if (irangp.lt.0) then
      vacapt(kk) = rtp(icapt(kk),ii)
    else
      call parhis(icapt(kk), ircapt(kk), rtp(1,ii), vacapt(kk))
      !==========
    endif
  enddo
  if (irangp.le.0) then
    write(impush(ii),1010) ntcabs, ttcabs, (vacapt(kk),kk=1,ncapts)
  endif
enddo

! WARNING: The format must be modified in case of more than 9 monitoring points

1010 format(i10,10e17.9)

!===============================================================================
! 4.  Close files
!===============================================================================

if (ntcabs.eq.ntmabs .and. irangp.le.0) then
  do ii = 1, nvarpp
    close(impush(ii))
  enddo
endif


!----
! End
!----

return
end subroutine ushist

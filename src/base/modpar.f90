!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine modpar &
!================

 ( ntcabs , ntmabs )

!===============================================================================
! Purpose:
! -------

!    Modify ntmabs during the calculation for interactive stop.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ntcabs           ! i  ! <-- ! absolute current time step number              !
! ntmabs           ! i  ! <-> ! absolute final time step number                !
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
use parall

!===============================================================================

implicit none

! Arguments

integer ntcabs , ntmabs

! Local variables

integer irangs, lng, itmp(1)
logical exstp

!===============================================================================

! Only one rank needs to test this (and broadcast later).

if (irangp.le.0) then

  !---> Emergency stop

  inquire (file=ficstp, exist=exstp)

  ! If a ficstp file is present

  if (exstp) then

    ! Read the (absolute) number of iterations

    open(file=ficstp, unit=impstp)
    read(impstp, *, err=5200, end=5200)
 5200     read (impstp,*,err=5100,end=5100)ntmabs
 5100     continue
    close (impstp,status='delete')

    ! Compare elapsed and maximum available time;
    ! modify ficstp if necessary.

    if(ntcabs.gt.ntmabs)then
      ntmabs = ntcabs
    endif

    ! Output

    write (nfecra,1000) ntcabs,ntmabs

    open (file=ficstp//'_updated', unit=impstp)
    write (impstp,1000) ntcabs,ntmabs
    close (impstp)
  endif

endif

! In parallel, broadcast
if (irangp.ge.0) then
  irangs  = 0
  lng     = 1
  itmp(1) = ntmabs
  call parbci(irangs,lng,itmp)
  ntmabs = itmp(1)
endif

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'=============================================================',/,&
'            NTCABS COURANT  = ', i10,                          /,&
'            NTMABS RESET TO = ', i10,                          /,&
'=============================================================',/,&
                                                                /)
#else

 1000 format(/,                                                   &
'=============================================================',/,&
'            NTCABS CURRENT  = ', i10,                          /,&
'            NTMABS RESET TO = ', i10,                          /,&
'=============================================================',/,&
                                                                /)
#endif

end subroutine

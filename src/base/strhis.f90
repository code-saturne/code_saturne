!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine strhis &
!================

 ( modhis )

!===============================================================================
! Purpose:
! -------

! Write time data for structures

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! modhis           ! i  ! <-- ! 0 or 1: initialize/output; 2: finalize         !
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
use optcal
use alstru
use parall

!===============================================================================

implicit none

! Arguments

integer          modhis

! Local variables

integer          nbname
parameter        (nbname=12)
character(len=300) :: nompre, nenvar
character(len=80) :: namevr(nbname)
integer          ii, jj, ii1, ii2, lpre, lnam, tplnum
double precision, dimension(:), allocatable :: vartmp

! Time plot number shift (in case multiple routines define plots)

integer  nptpl
data     nptpl /0/
save     nptpl

! Number of passes in this routine

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 0. Local initializations
!===============================================================================

ipass = ipass + 1

!--> Only rank 0 outputs, and nothing to do if we have no structure

if (irangp.gt.0 .or. nbstru.le.0) return

!--> If no history was output
if (ipass.eq.1.and.modhis.eq.2) return

if (ipass.eq.1) then
  call tplnbr(nptpl)
  !==========
endif

!===============================================================================
! 2. Initialize output
!===============================================================================

! Create directory if required
if (ipass.eq.1 .and. irangp.le.0) then
  call csmkdr(emphis, len(emphis))
  !==========
endif

if (ipass.eq.1 .and. irangp.le.0) then

  namevr(1 ) = "displacement x"
  namevr(2 ) = "displacement y"
  namevr(3 ) = "displacement z"
  namevr(4 ) = "velocity x"
  namevr(5 ) = "velocity y"
  namevr(6 ) = "velocity z"
  namevr(7 ) = "acceleration x"
  namevr(8 ) = "acceleration y"
  namevr(9 ) = "acceleration z"
  namevr(10) = "force x"
  namevr(11) = "force y"
  namevr(12) = "force z"

  ! --> write one file per structure info type

  do ii = 1, nbname

    ! plot prefix
    nompre = ' '
    call verlon(emphis, ii1, ii2, lpre)
    !==========
    nompre(1:lpre) = emphis(ii1:ii2)
    call verlon(prehis, ii1, ii2, lnam)
    !==========
    nompre(lpre+1:lpre+lnam) = prehis(ii1:ii2)
    call verlon(nompre, ii1, ii2, lpre)
    !==========
    nompre(lpre+1:lpre+4) = 'str_'
    call verlon(nompre, ii1, ii2, lpre)
    !==========

    ! plot name
    nenvar = namevr(ii)
    call verlon(nenvar,ii1,ii2,lnam)
    !==========

    tplnum = nptpl + ii

    call tpsini(tplnum, nenvar, nompre, tplfmt, idtvar, &
    !==========
                nbstru, xmstru, xcstru, xkstru, lnam, lpre)

  enddo

endif

!===============================================================================
! 3. Output results
!===============================================================================

if ((modhis.eq.0 .or. modhis.eq.1) .and. irangp.le.0) then

  allocate(vartmp(nbstru))

  tplnum = nptpl + 1
  do jj = 1, nbstru
    vartmp(jj) = xstr(1, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 2
  do jj = 1, nbstru
    vartmp(jj) = xstr(2, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 3
  do jj = 1, nbstru
    vartmp(jj) = xstr(3, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 4
  do jj = 1, nbstru
    vartmp(jj) = xpstr(1, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 5
  do jj = 1, nbstru
    vartmp(jj) = xpstr(2, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 6
  do jj = 1, nbstru
    vartmp(jj) = xpstr(3, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 7
  do jj = 1, nbstru
    vartmp(jj) = xppstr(1, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 8
  do jj = 1, nbstru
    vartmp(jj) = xppstr(2, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 9
  do jj = 1, nbstru
    vartmp(jj) = xppstr(3, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 10
  do jj = 1, nbstru
    vartmp(jj) = forstr(1, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 11
  do jj = 1, nbstru
    vartmp(jj) = forstr(2, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  tplnum = nptpl + 12
  do jj = 1, nbstru
    vartmp(jj) = forstr(3, jj)
  enddo
  call tplwri(tplnum, tplfmt, nbstru, ntcabs, ttcabs, vartmp)

  deallocate(vartmp)

endif

!===============================================================================
! 4. Close output
!===============================================================================

if (modhis.eq.2 .and. irangp.le.0) then

  do ii = 1, nbname
    tplnum = nptpl + ii
    call tplend(tplnum, tplfmt)
    !==========
  enddo

endif

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine

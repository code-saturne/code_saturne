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

subroutine ecrhis &
!================

 ( ndim   , ncelet , ncel , modhis , xyzcen , ra )

!===============================================================================
! Purpose:
! -------

! Write plot data

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! modhis           ! i  ! <-- ! 0 or 1: initialize/output; 2: finalize         !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! ra(*)            ! ra ! <-- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstnum
use optcal
use parall

!===============================================================================

implicit none

! Arguments

integer          ndim, ncelet, ncel
integer          modhis
double precision xyzcen(ndim, ncelet)
double precision, dimension(*), target :: ra

! Local variables

character        nompre*300, nomhis*300
integer          tplnum, ii, ii1, ii2, lpre, lnom, lng
integer          icap, ncap, ipp, ira
integer          idivdt, iel
integer          nbcap(nvppmx)
double precision varcap(ncaptm)

integer, dimension(:), allocatable :: lsttmp
double precision, dimension(:), allocatable, target :: momtmp
double precision, dimension(:,:), allocatable :: xyztmp

double precision, dimension(:), pointer :: varptr => null()

! Time plot number shift (in case multiple routines define plots)

integer  nptpl
data     nptpl /0/
save     nptpl

! Number of passes in this routine

integer  ipass
data     ipass /0/
save     ipass

!===============================================================================
! 0. Local initializations
!===============================================================================

ipass = ipass + 1

!--> If no probe data has been output or there are no probes, return
if ((ipass.eq.1 .and. modhis.eq.2) .or. ncapt.eq.0) return

if (ipass.eq.1) then
  call tplnbr(nptpl)
  !==========
endif

!===============================================================================
! 1. Search for neighboring nodes -> nodcap
!===============================================================================

if (ipass.eq.1) then

  do ii = 1, ncapt
    call findpt                                                        &
    !==========
    (ncelet, ncel, xyzcen,                                             &
     xyzcap(1,ii), xyzcap(2,ii), xyzcap(3,ii), nodcap(ii), ndrcap(ii))
  enddo

endif

!===============================================================================
! 2. Initialize output
!===============================================================================

! Create directory if required
if (ipass.eq.1 .and. irangp.le.0) then
  call csmkdr(emphis, len(emphis))
  !==========
endif

if (ipass.eq.1) then

  ! Number of probes per variable
  do ipp = 2, nvppmx
    nbcap(ipp) = ihisvr(ipp,1)
    if (nbcap(ipp).lt.0) nbcap(ipp) = ncapt
  enddo

  allocate(lsttmp(ncapt))
  allocate(xyztmp(3, ncapt))

  ! Initialize one output per variable

  do ipp = 2, nvppmx

    if (ihisvr(ipp,1) .gt. 0) then
      do ii=1, ihisvr(ipp,1)
        lsttmp(ii) = ihisvr(ipp,ii+1)
        if (irangp.lt.0 .or. irangp.eq.ndrcap(ihisvr(ipp, ii+1))) then
          xyztmp(1, ii) = xyzcen(1, nodcap(ihisvr(ipp, ii+1)))
          xyztmp(2, ii) = xyzcen(2, nodcap(ihisvr(ipp, ii+1)))
          xyztmp(3, ii) = xyzcen(3, nodcap(ihisvr(ipp, ii+1)))
        endif
        if (irangp.ge.0) then
          lng = 3
          call parbcr(ndrcap(ihisvr(ipp,ii+1)), lng , xyztmp(1, ii))
          !==========
        endif
      enddo
    else
      do ii = 1, nbcap(ipp)
        lsttmp(ii) = ii
        if (irangp.lt.0 .or. irangp.eq.ndrcap(ii)) then
          xyztmp(1, ii) = xyzcen(1, nodcap(ii))
          xyztmp(2, ii) = xyzcen(2, nodcap(ii))
          xyztmp(3, ii) = xyzcen(3, nodcap(ii))
        endif
        if (irangp.ge.0) then
          lng = 3
          call parbcr(ndrcap(ii), lng , xyztmp(1, ii))
          !==========
        endif
      enddo
    endif

    if (nbcap(ipp) .gt. 0) then

      if (irangp.le.0) then

        ! plot prefix
        nompre = ' '
        call verlon(emphis, ii1, ii2, lpre)
        !==========
        nompre(1:lpre) = emphis(ii1:ii2)
        call verlon(prehis, ii1, ii2, lnom)
        !==========
        nompre(lpre+1:lpre+lnom) = prehis(ii1:ii2)
        call verlon(nompre, ii1, ii2, lpre)
        !==========

        ! plot name
        nomhis = ' '
        call verlon(nomvar(ipp), ii1, ii2, lnom)
        !==========
        nomhis(1:lnom) = nomvar(ipp)(ii1:ii2)

        tplnum = nptpl + ipp
        call tppini(tplnum, nomhis, nompre, tplfmt, idtvar, nthsav, tplflw, &
        !==========
                    nbcap(ipp), lsttmp(1), xyzcap(1,1), lnom, lpre)

      endif ! (irangp.le.0)

    endif

  enddo

  deallocate(lsttmp)
  deallocate(xyztmp)

endif

!===============================================================================
! 3. Output results
!===============================================================================

if (modhis.eq.0 .or. modhis.eq.1) then

  do ipp = 2, nvppmx
    if (ihisvr(ipp,1).ne.0) then
      ira = abs(ipp2ra(ipp))

      ! For moments, we must divide by the cumulative time
      idivdt = ippmom(ipp)
      if (idivdt.eq.0) then
        varptr => ra(ira:ira+ncel)
      else
        allocate(momtmp(ncel))
        varptr => momtmp
        if (idivdt.gt.0) then
          do iel = 1, ncel
            momtmp(iel) = ra(ira+iel-1)/max(ra(idivdt+iel-1),epzero)
          enddo
        elseif (idivdt.lt.0) then
          do iel = 1, ncel
            momtmp(iel) = ra(ira+iel-1)/max(dtcmom(-idivdt),epzero)
          enddo
        endif
      endif

      if (ihisvr(ipp,1).lt.0) then
        do icap = 1, ncapt
          if (irangp.lt.0) then
            varcap(icap) = varptr(nodcap(icap))
          else
            call parhis(nodcap(icap), ndrcap(icap), varptr, varcap(icap))
            !==========
          endif
        enddo
        ncap = ncapt
      else
        do icap = 1, ihisvr(ipp,1)
          if (irangp.lt.0) then
            varcap(icap) = varptr(nodcap(ihisvr(ipp,icap+1)))
          else
            call parhis(nodcap(ihisvr(ipp,icap+1)), &
            !==========
                        ndrcap(ihisvr(ipp,icap+1)), &
                        varptr, varcap(icap))
          endif
        enddo
        ncap = ihisvr(ipp,1)
      endif

      if (idivdt.ne.0) then
        deallocate(momtmp)
      endif

      if (irangp.le.0 .and. ncap.gt.0) then
        tplnum = nptpl + ipp
        call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
        !==========
      endif
    endif
  enddo

endif

!===============================================================================
! 4. Close output
!===============================================================================

if (modhis.eq.2) then

  do ipp = 2, nvppmx
    if (ihisvr(ipp,1).lt.0) then
      ncap = ncapt
    else
      ncap = ihisvr(ipp,1)
    endif
    if (irangp.le.0 .and. ncap.gt.0) then
      tplnum = nptpl + ipp
      call tplend(tplnum, tplfmt)
      !==========
    endif
  enddo

endif

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine

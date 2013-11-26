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

!> \file ecrhis.f90
!> \brief Write plot data
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in]     nproce        total number of physical properties
!> \param[in]     modhis        0 or 1: initialize/output; 2: finalize
!> \param[in,out] rtp           calculated variables at cell centers
!>                              (at current time step)
!______________________________________________________________________________

subroutine ecrhis &
!================

 ( nvar , nproce , modhis , rtp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstnum
use optcal
use parall
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar, nproce
integer          modhis
double precision, dimension(*), target :: rtp(ncelet,*)

! Local variables

character        nompre*300, nomhis*300
logical          lprev
integer          tplnum, ii, ii1, ii2, lpre, lnom, lng
integer          icap, ncap, ipp, ivar, iprop
integer          iel, isou, keymom, mom_id, f_id
integer          nbcap(nvppmx)
double precision varcap(ncaptm)

integer, dimension(:), allocatable :: lsttmp
double precision, dimension(:), allocatable, target :: momtmp
double precision, dimension(:,:), allocatable :: xyztmp

double precision, dimension(:,:), pointer :: val_v
double precision, dimension(:), pointer :: val_s
double precision, dimension(:), pointer :: num_s, div_s

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

  ! Loop on variables

  do ivar = 1, nvar

    ipp = ipprtp(ivar)

    if (ihisvr(ipp,1).ne.0) then

      val_s => rtp(:, ivar)

      if (ihisvr(ipp,1).lt.0) then
        do icap = 1, ncapt
          if (irangp.lt.0) then
            varcap(icap) = rtp(nodcap(icap), ivar)
          else
            call parhis(nodcap(icap), ndrcap(icap), val_s, varcap(icap))
            !==========
          endif
        enddo
        ncap = ncapt
      else
        do icap = 1, ihisvr(ipp,1)
          if (irangp.lt.0) then
            varcap(icap) = val_s(nodcap(ihisvr(ipp,icap+1)))
          else
            call parhis(nodcap(ihisvr(ipp,icap+1)), &
            !==========
                        ndrcap(ihisvr(ipp,icap+1)), &
                        val_s, varcap(icap))
          endif
        enddo
        ncap = ihisvr(ipp,1)
      endif

      if (irangp.le.0 .and. ncap.gt.0) then
        tplnum = nptpl + ipp
        call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
        !==========
      endif
    endif
  enddo

  ! Loop on physical properties

  call field_get_key_id('moment_dt', keymom)

  do iprop= 1, nproce

    ipp = ipppro(iprop)

    if (ihisvr(ipp,1).ne.0) then

      f_id = iprpfl(iprop)
      lprev = .false.

      if (f_id .eq. -1 .and. iroma.gt.0) then
        if (iprop.eq.ipproc(iroma)) then
          f_id = iprpfl(ipproc(irom))
          lprev = .true.
        endif
      endif
      if (f_id.eq.-1) cycle

      call field_get_key_int(f_id, keymom, mom_id)

      ! For moments, we must divide by the cumulative time

      if (mom_id.eq.-1) then
        if (.not. lprev) then
          call field_get_val_s(f_id, val_s)
        else
          call field_get_val_prev_s(f_id, val_s)
        endif
      else
        allocate(momtmp(ncel))
        val_s => momtmp
        if (.not. lprev) then
          call field_get_val_s(f_id, num_s)
        else
          call field_get_val_prev_s(f_id, num_s)
        endif
        if (mom_id.ge.0) then
          call field_get_val_s(mom_id, div_s)
          do iel = 1, ncel
            momtmp(iel) = num_s(iel)/max(div_s(iel),epzero)
          enddo
        else
          do iel = 1, ncel
            momtmp(iel) = num_s(iel)/max(dtcmom(-mom_id -1),epzero)
          enddo
        endif
      endif

      if (ihisvr(ipp,1).lt.0) then
        do icap = 1, ncapt
          if (irangp.lt.0) then
            varcap(icap) = val_s(nodcap(icap))
          else
            call parhis(nodcap(icap), ndrcap(icap), val_s, varcap(icap))
            !==========
          endif
        enddo
        ncap = ncapt
      else
        do icap = 1, ihisvr(ipp,1)
          if (irangp.lt.0) then
            varcap(icap) = val_s(nodcap(ihisvr(ipp,icap+1)))
          else
            call parhis(nodcap(ihisvr(ipp,icap+1)), &
            !==========
                        ndrcap(ihisvr(ipp,icap+1)), &
                        val_s, varcap(icap))
          endif
        enddo
        ncap = ihisvr(ipp,1)
      endif

      if (mom_id.ne.-1) then
        deallocate(momtmp)
      endif

      if (irangp.le.0 .and. ncap.gt.0) then
        tplnum = nptpl + ipp
        call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
        !==========
      endif
    endif
  enddo

  ! Other scalar fields

  ipp = ippdt
  if (ihisvr(ipp,1).ne.0) then

    call field_get_id('dt', f_id)
    call field_get_val_s(f_id, val_s)

    if (ihisvr(ipp,1).lt.0) then
      do icap = 1, ncapt
        if (irangp.lt.0) then
          varcap(icap) = val_s(nodcap(icap))
        else
          call parhis(nodcap(icap), ndrcap(icap), val_s, varcap(icap))
          !==========
        endif
      enddo
      ncap = ncapt
    else
      do icap = 1, ihisvr(ipp,1)
        if (irangp.lt.0) then
          varcap(icap) = val_s(nodcap(ihisvr(ipp,icap+1)))
        else
          call parhis(nodcap(ihisvr(ipp,icap+1)), &
          !==========
                      ndrcap(ihisvr(ipp,icap+1)), &
                      val_s, varcap(icap))
        endif
      enddo
      ncap = ihisvr(ipp,1)
    endif

    if (irangp.le.0 .and. ncap.gt.0) then
      tplnum = nptpl + ipp
      call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
      !==========
    endif
  endif

  ! Other vector fields (currently, only tpucou)

  if (     ihisvr(ipptx,1).ne.0 .or. ihisvr(ippty,1).ne.0  &
      .or. ihisvr(ipptz,1).ne.0) then
    f_id = idtten
  else
    f_id = -1
  endif

  if (f_id .ge. 0) then

    call field_get_val_v(f_id, val_v)

    do isou = 1, 3

      if (isou.eq.1) then
        ipp = ipptx
      else if (isou.eq.2) then
        ipp = ippty
      else if (isou.eq.3) then
        ipp = ipptz
      endif

      if (ihisvr(ipp,1).lt.0) then
        do icap = 1, ncapt
          if (irangp.lt.0 .or. ndrcap(icap).eq.irangp) then
            varcap(icap) = val_v(isou, nodcap(icap))
          endif
          if (irangp.ge.0) then
            lng = 1
            call parbcr(ndrcap(icap), lng, varcap(icap))
            !==========
          endif
        enddo
        ncap = ncapt
      else if (ihisvr(ipp,1).gt.0) then
        do icap = 1, ihisvr(ipp,1)
          if (irangp.lt.0 .or. ndrcap(icap).eq.irangp) then
            varcap(icap) = val_v(isou, nodcap(ihisvr(ipp,icap+1)))
          endif
          if (irangp.ge.0) then
            lng = 1
            call parbcr(ndrcap(icap), lng, varcap(icap))
            !==========
          endif
        enddo
        ncap = ihisvr(ipp,1)
      else
        ncap = 0
      endif

      if (irangp.le.0 .and. ncap.gt.0) then
        tplnum = nptpl + ipp
        call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
        !==========
      endif

    enddo

  endif

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

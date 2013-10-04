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

subroutine laghis &
!================

 ( ndim   , ncelet , ncel   , modhis , nvlsta ,  &
   xyzcen , volume , statis , stativ )

!===============================================================================
! Purpose:
! -------

! Write plot data for the Lagrangian module

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
use parall
use cstnum
use optcal
use lagpar
use lagran

!===============================================================================

implicit none

! Arguments

integer          ndim, ncelet, ncel
integer          nvlsta
integer          modhis
double precision xyzcen(ndim,ncelet) , volume(ncelet)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)

! Local variables

character        nomhis*300, nompre*300
integer          ii, ii1, ii2, lpre, lnom, inam1, inam2, lng
integer          icap, ncap, ipp
integer          iel, ivarl, tplnum
integer          ipas, ilpd1, il, ilfv1, icla
integer          iokhis

double precision varcap(ncaptm)
double precision dmoy

integer, dimension(:), allocatable :: lsttmp
double precision, dimension(:), allocatable, target :: moytmp
double precision, dimension(:,:), allocatable :: xyztmp

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

!--> If no probe data has been output or there are no probes, return
if ((ipass.eq.1 .and. modhis.eq.2) .or. ncapt.eq.0) return

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

if (ipass.eq.1) then

  allocate(lsttmp(ncapt))
  allocate(xyztmp(3, ncapt))

  do ii=1, ncapt
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

  ! Initialize one output per variable

  do ipas  = 1, 1+nbclst

    do ipp = 1, 2*nvlsta

      if (ipp .le. nvlsta) then
        ivarl = (ipas-1)*nvlsta+ipp
        ilpd1 = (ipas-1)*nvlsta+ilpd
        ilfv1 = (ipas-1)*nvlsta+ilfv
        icla  = ipas -1
      else
        ivarl = (ipas-1)*nvlsta+(ipp-nvlsta)
        ilpd1 = (ipas-1)*nvlsta+ilpd
        ilfv1 = (ipas-1)*nvlsta+ilfv
        icla  = ipas -1
      endif

      iokhis = 0
      if (ipp.le.nvlsta) then
        if (ihslag(ipp).ge.1) iokhis = 1
      else
        if ((ipp-nvlsta).ne.ilpd .and. ihslag(ipp-nvlsta).eq.2) iokhis = 1
      endif

      if (iokhis.eq.1) then

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
          nompre(lpre+1:lpre+4) = 'Lag_'
          call verlon(nompre, ii1, ii2, lpre)
          !==========

          ! plot name
          nomhis = ' '
          if (ipas.eq.1) then
            if (ipp.le.nvlsta) then
              nomhis = nomlag(ipp)
            else
              nomhis = nomlav(ipp-nvlsta)
            endif
          else
            if (ipp.le.nvlsta) then
              write(nomhis, '(a8, a4, i3)') nomlag(ipp), '_grp', icla
            else
              write(nomhis, '(a8, a4, i3)') nomlav(ipp-nvlsta), '_grp', icla
            endif
          endif
          call verlon(nomhis, inam1, inam2, lnom)
          !==========
          call undscr(inam1, inam2, nomhis)
          !==========

          tplnum = nptpl +(ipas-1)*2*nvlsta + ipp + 1

          call tppini(tplnum, nomhis, nompre, tplfmt, idtvar, nthsav, tplflw, &
          !==========
                      ncapt, lsttmp(1), xyzcap(1,1), lnom, lpre)

        endif

      endif

    enddo

  enddo

  deallocate(lsttmp)
  deallocate(xyztmp)

endif

!===============================================================================
! 3. Output results
!===============================================================================

if (modhis.eq.0 .or. modhis.eq.1) then

  allocate(moytmp(ncel))

  do ipas  = 1, 1+nbclst

    ! Mean

    do il = 1, nvlsta

      ivarl = (ipas-1)*nvlsta+il
      ilpd1 = (ipas-1)*nvlsta+ilpd
      ilfv1 = (ipas-1)*nvlsta+ilfv
      icla  = ipas -1

      ! Pour l'instant on fait des chrono sur toutes les variables Stat. Lag.
      ! et sur tout les capteurs

      if (ihslag(ivarl).ge. 1) then

        if (ivarl .ne. ilpd1 .and. ivarl .ne. ilfv1) then
          do iel = 1, ncel
            if (statis(iel, ilpd1) .gt. seuil) then
              moytmp(iel) = statis(iel, ivarl) / statis(iel, ilpd1)
            else
              moytmp(iel) = 0.d0
            endif
          enddo

        else if (ivarl.eq.ilpd1) then
          do iel=1, ncel
            moytmp(iel) = statis(iel, ivarl)
          enddo
        else
          do iel=1, ncel
            if (npst.gt.0) then
              moytmp(iel) = statis(iel, ivarl) / (dble(npst)*volume(iel))
            else
              moytmp(iel) = 0.d0
            endif
          enddo
        endif

        do icap = 1, ncapt
          if (irangp.lt.0) then
            varcap(icap) = moytmp(nodcap(icap))
          else
            call parhis(nodcap(icap), ndrcap(icap), moytmp, varcap(icap))
            !==========
          endif
        enddo
        ncap = ncapt

        if (irangp.le.0 .and. ncap.gt.0) then
          tplnum = nptpl + (ipas-1)*2*nvlsta + il + 1
          call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
          !==========
        endif

      endif
    enddo

    ! Variance

    do il = 1, nvlsta-1

      ivarl = (ipas-1)*nvlsta+il
      ilpd1 = (ipas-1)*nvlsta+ilpd
      ilfv1 = (ipas-1)*nvlsta+ilfv
      icla  = ipas -1

      ! Pour l'instant on fait des chrono sur toutes les variables Stat. Lag.
      ! et sur tout les capteurs

      if (ihslag(ivarl).eq. 2) then
        do iel = 1, ncel

          if (ivarl.ne.ilfv) then
            if (statis(iel, ilpd1).gt.seuil) then
              moytmp(iel) =   stativ(iel, ivarl)/statis(iel, ilpd1)     &
                            -(statis(iel, ivarl)/statis(iel, ilpd1)     &
                              *statis(iel, ivarl)/statis(iel, ilpd1))
            else
              moytmp(iel) = zero
            endif
          else
            if (statis(iel, ilpd1).gt.seuil .and. npst.gt.0) then
              dmoy = statis(iel, ivarl) /(dble(npst)*volume(iel))
              moytmp(iel) =  stativ(iel, ivarl) / (dble(npst) * volume(iel))  &
                            -dmoy*dmoy

            else if (statis(iel, ilpd1).gt.seuil .and. iplas.ge.idstnt) then
              dmoy =  statis(iel, ivarl) / volume(iel)
              moytmp(iel) = stativ(iel, ilfv) / volume(iel) -dmoy*dmoy
            else
              moytmp(iel) = zero
            endif
          endif
          moytmp(iel) = sqrt(max(zero, moytmp(iel)))
        enddo

        do icap = 1, ncapt
          if (irangp.lt.0) then
            varcap(icap) = moytmp(nodcap(icap))
          else
            call parhis(nodcap(icap), ndrcap(icap), moytmp, varcap(icap))
            !==========
          endif
        enddo
        ncap = ncapt

        if (irangp.le.0 .and. ncap.gt.0) then
          tplnum = nptpl + (ipas-1)*2*nvlsta + il + nvlsta + 1
          call tplwri(tplnum, tplfmt, ncap, ntcabs, ttcabs, varcap)
          !==========
        endif

      endif
    enddo

  enddo

  deallocate(moytmp)

endif

!===============================================================================
! 4. Close output
!===============================================================================

if (modhis.eq.2) then

  if (irangp.le.0 .and. ncapt.gt.0) then

    ! --> nombre de pas de temps enregistres

    do ipas  = 1, 1+nbclst

      do ipp = 1, 2*nvlsta

        iokhis = 0
        if (ipp.le.nvlsta) then
          if (ihslag(ipp).ge.1) iokhis = 1
        else
          if ((ipp-nvlsta).ne.ilpd .and. ihslag(ipp-nvlsta).eq.2) iokhis = 1
        endif

        if (iokhis.eq.1) then
          tplnum = nptpl + (ipas-1)*2*nvlsta + ipp + 1
          call tplend(tplnum, tplfmt)
          !==========
        endif

      enddo
    enddo

  endif

endif

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine

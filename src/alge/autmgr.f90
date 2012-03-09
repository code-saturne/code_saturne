!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine autmgr &
!================

    ( igr    , isym   , nagmax ,                                  &
      ncelf  , ncelfe , nfacf  , iwarnp , ifacef ,                &
      daf    , xaf    , surfaf , volumf , xyzfin ,                &
      iordf  , irscel ,                                           &
      indic  , inombr , irsfac , indicf , w1     , w2 )

!===============================================================================
! Purpose:
! --------

!  Algebraic multigrid:
!  build a coarse grid level from the previous level using
!  an automatic criterion.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! igr              ! i  ! <-- ! coarse grid level                              !
! isym             ! i  ! <-- ! 1: symmetric matrix; 2: non-symmetric matrix   !
! nagmax           ! i  ! <-- ! max. fine cells per coarse cell                !
! ncelf            ! i  ! <-- ! number of cells in fine grid                   !
! ncelfe           ! i  ! <-- ! extended number of cells in fine grid          !
! nfacf            ! i  ! <-- ! number of interior faces in fine grid          !
! iwarnp           ! i  ! <-- ! verbosity                                      !
! ifacef(2, nfacf) ! ia ! <-- ! fine grid interior face -> cells connectivity  !
! daf(ncelf)       ! ra ! <-- ! fine grid matrix diagonal terms                !
! xaf(nfacf*isym)  ! ra ! <-- ! fine grid matrix extra-diagonal terms          !
! surfaf(3, nfacf) ! ra ! <-- ! fine grid face surfaces                        !
! volumf(ncelf)    ! ra ! <-- ! fine grid cell volumes                         !
! xyzfin(3, ncelf) ! ra ! <-- ! fine grid cell centers                         !
! irscel(ncelfe)   ! ia ! --> ! fine cell -> coarse cell                       !
! indic(ncelfe)    ! ia ! --- ! work array                                     !
! inombr(ncelfe)   ! ia ! --- ! work array                                     !
! irsfac(nfacf)    ! ia ! --- ! fine face -> coarse face (work array)          !
! indicf(nfacf)    ! ia ! --- ! face merging indicator                         !
! icelfa(2*nfacf)  ! ia ! --- ! fine grid cells -> faces connectivity          !
! icelce(2*nfacf)  ! ia ! --- ! fine mesh cells -> cells connectivity          !
! rw(ncelf)        ! ra ! --- ! work array                                     !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use entsor
use parall

!===============================================================================

implicit none

! Arguments

integer          igr, isym, nagmax
integer          ncelf, ncelfe, nfacf
integer          iwarnp

integer          ifacef(2, nfacf)
integer          iordf(nfacf), irscel(ncelfe)
integer          indic(ncelfe), inombr(ncelfe)
integer          indicf(nfacf), irsfac(nfacf)

double precision daf(ncelfe), xaf(*) ! xaf of size nfacf or nfacf*2
double precision surfaf(3,nfacf), volumf(ncelfe)
double precision xyzfin(3,ncelfe)
double precision w1(ncelfe), w2(ncelfe)

! Local variables

integer          ncelg, ncelgg, icel, ifac , ifac1, ifac2, icelg
integer          nfacn,nfacnr,npass,npasmx
integer          inditt, noaglo, ngros, incvoi, iagmax
integer          i, j, imin, imax

double precision epslon, xaf1, xaf2

integer, allocatable, dimension(:) :: ihist
double precision, allocatable, dimension(:) :: critr

!===============================================================================

! Default parameters

epslon = +1.d-6
ngros  = 8
npasmx = 10
incvoi = 1

!===============================================================================
! 1. Initialization
!===============================================================================

iagmax = 1

do icel = 1, ncelfe
  indic(icel) = -1
  irscel(icel) = 0
  inombr(icel) = 1
enddo

do ifac = 1, nfacf
  indicf(ifac) = ifac
  irsfac(ifac) = 0
enddo

! Compute cardinality (number of neighbors for each cell -1)

do ifac = 1, nfacf
  i = ifacef(1,ifac)
  j = ifacef(2,ifac)
  indic(i) = indic(i) + 1
  indic(j) = indic(j) + 1
enddo

ncelg  = 0
nfacnr = nfacf
npass  = 0
noaglo = ncelf

allocate(critr(nfacf))

critr = 2.d0

! Passes

do while (npass .lt. npasmx)

  npass = npass+1
  nfacn = nfacnr
  iagmax= iagmax +1
  iagmax= min(iagmax, nagmax)

  do ifac=1,nfacn
    irsfac(ifac) = indicf(ifac)
    indicf(ifac) = 0
  enddo
  if (nfacn .lt. nfacf) then
    do ifac = nfacn+1, nfacf
      indicf(ifac) = 0
      irsfac(ifac) = 0
    enddo
  endif

  if (iwarnp .gt. 3) then
    write(nfecra,2001) npass, nfacnr, noaglo
  endif

  ! Increment number of neighbors

  do icel = 1, ncelf
    indic(icel) = indic(icel) + incvoi
  enddo

  ! Initialize non-eliminated faces

  nfacnr = 0

  do ifac1 = 1, nfacn

    ifac = irsfac(ifac1)
    i = ifacef(1,ifac)
    j = ifacef(2,ifac)

    ! Exclude faces on parallel or periodic boundary, so as not to
    ! coarsen the grid across those boundaries (which would change
    ! the communication pattern and require a more complex algorithm).

    if (i.le.ncelf .and. j.le.ncelf) then
      xaf1 = xaf((ifac-1)*isym + 1)
      xaf2 = xaf(ifac*isym)
      xaf1 = max(-xaf1, 1.d-15)
      xaf2 = max(-xaf2, 1.d-15)
      critr(ifac1) = (daf(i)/indic(i))*(daf(j)/indic(j))/(xaf1*xaf2)
    endif

  enddo

  ! order faces by criteria (0 to n-1 numbering)

  call clmlgo(nfacn, critr, iordf)
  !==========

  ! Loop on non-eliminated faces

  do ifac1 = 1, nfacn

    ifac2 = iordf(ifac1) + 1
    ifac = irsfac(ifac2)
    i = ifacef(1,ifac)
    j = ifacef(2,ifac)

    ! Exclude faces on parallel or periodic boundary, so as not to
    ! coarsen the grid across those boundaries (which would change
    ! the communication pattern and require a more complex algorithm).

    if (i.le.ncelf .and. j.le.ncelf) then

      inditt = 0

      if (critr(ifac2).lt.(1.d0-epslon) .and. irscel(i)*irscel(j).le.0) then

        if (irscel(i).gt.0 .and. irscel(j).le.0) then
          if(inombr(irscel(i)) .le. iagmax) then
            irscel(j) = irscel(i)
            inombr(irscel(i)) = inombr(irscel(i)) +1
            inditt = inditt +1
          endif
        else if (irscel(i).le.0 .and. irscel(j).gt.0) then
          if (inombr(irscel(j)).le.iagmax) then
            irscel(i) = irscel(j)
            inombr(irscel(j)) = inombr(irscel(j)) + 1
            inditt = inditt +1
          endif
        else if (irscel(i).le.0 .and. irscel(j).le.0) then
          ncelg = ncelg+1
          irscel(i) = ncelg
          irscel(j) = ncelg
          inombr(ncelg) = inombr(ncelg) +1
          inditt = inditt +1
        endif

      endif

      if (inditt.eq.0 .and. irscel(i)*irscel(j).le.0) then
        nfacnr = nfacnr + 1
        indicf(nfacnr) = ifac
        iordf(nfacnr) = nfacnr - 1
      endif

    endif

  enddo

  ! Check the number of coarse cells created

  noaglo = 0
  do icel=1,ncelf
    if (irscel(icel).le.0) noaglo = noaglo+1
  enddo

  ! Additional passes if agglomeration is insufficient

  if (     (noaglo.eq.0) &
      .or. ((ncelg+noaglo)*ngros .lt. ncelf) &
      .or. (nfacnr.eq.0)) then
    npasmx = npass
  endif

enddo ! Loop on passes

! Finish assembly

do icel = 1, ncelf
  if (irscel(icel).le.0) then
    ncelg = ncelg+1
    irscel(icel) = ncelg
  endif
enddo

! Various checks

imax = 0
imin = 2*ncelf
do icelg =1,ncelg
  imax = max(imax, inombr(icelg))
  imin = min(imin, inombr(icelg))
enddo

ncelgg = ncelg
if (irangp .ge. 0) then
  call parcmn(imin)
  call parcmx(imax)
  call parcpt(ncelgg)
endif

if (iwarnp.gt.3) then

  write(nfecra,2002) imin, imax, nagmax
  write(nfecra,2003)
  noaglo=imax-imin+1
  if (noaglo.gt.0) then
    allocate(ihist(noaglo))
    do i = 1, noaglo
      ihist(i) = 0
    enddo
    do icelg = 1, ncelg
      do i = 1, noaglo
        if (inombr(icelg).eq.(imin+i-1))then
          ihist(i)=ihist(i)+1
        endif
      enddo
    enddo
    if (irangp .ge. 0) then
      call parism(noaglo, ihist)
    endif
    do i = 1, noaglo
      epslon = 100.d0*ihist(i)/ncelgg
      write(nfecra,2004) imin+i-1, epslon
    enddo

    deallocate(ihist)

  endif

endif

do icel = 1, ncelf
  indic(icel) = 0
enddo
do icel = 1, ncelf
  icelg = irscel(icel)
  indic(icelg) = indic(icelg) +1
enddo

i=0
j=2*ncelf
noaglo = 0
do icelg = 1, ncelg
  i = max(i, indic(icelg))
  j = min(j, indic(icelg))
  noaglo = noaglo + indic(icelg)
enddo

if (irangp .ge. 0) then
  call parcmn(j)
  call parcmx(j)
endif

if (iwarnp.gt.3) then
  write(nfecra,2005) j, i
endif

if (noaglo .ne. ncelf) then
  write(nfecra,*) ' Bug in autmgr, contact support.'
  call csexit(1)
endif

deallocate(critr)

!--------
! Formats
!--------

 2001 format(&
  '    autmgr: pass ', i3, ' nfacnr = ', i10, ' noaglo = ', i10)
 2002 format(&
  '    autmgr: inombr min = ', i10, ' max = ', i10, ' target = ', i10)
 2003 format(&
  '      histogram ')
 2004 format(&
  '        regroupment ', i10,' = ', e12.5,' %')
 2005 format(&
  '    autmgr: agglomeration min = ', i10, ' max= ', i10)
!==============================================================================

return
end subroutine


!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!===============================================================================
!-------------------------------------------------------------------------------
!> \file porcel.f90
!> \brief  This routine permits to calculate the porosity in wall-normal
!>  cells from the mean deposit height (which is only known
!>  at the boundary faces)
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    mdiam         mean particle diameter
!> \param[out]    porosi        deposition porosity
!_______________________________________________________________________________

subroutine porcel &
  ( mdiam, porosi)

  !===============================================================================

  !===============================================================================
  ! Module files
  !===============================================================================

  use paramx
  use numvar
  use entsor
  use optcal
  use cstphy
  use cstnum
  use lagran
  use lagdim
  use lagpar
  use pointe, only: parbor, itypfb
  use ppppar
  use coincl
  use parall
  use period
  use mesh
  use field
  use cs_c_bindings

  !===============================================================================

  implicit none

  ! Arguments

  double precision mdiam(ncelet)
  double precision porosi(ncelet)

  ! Local variables

  integer          f_id
  integer          ifac  , iel   , iel1, iel2, ii
  integer          inc   , iccocg
  integer          isou

  double precision xnorme

  integer, allocatable, dimension(:) :: itytmp
  double precision, allocatable, dimension(:) :: coefap, coefbp
  double precision, allocatable, dimension(:) :: masflu, depvol, distpw
  double precision, allocatable, dimension(:,:) :: q
  double precision, dimension(:), pointer :: crom
  double precision, dimension(:), pointer :: viscl

  integer nn , indic
  double precision dt, prod1, prod2, epsi

  !===============================================================================

  !===============================================================================
  ! 1. Initialization
  !===============================================================================

  ! Compute distance to deposition wall, and check deposition exists on
  ! at least one face.

  indic = 0

  allocate(itytmp(nfabor))

  do ifac = 1, nfabor
    if (      (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) &
        .and. parbor(ifac,ihdepm).gt.0.d0) then
      indic = 1
      itytmp(ifac) = iparoi
    else
      itytmp(ifac) = iindef
    endif
  enddo

  if (irangp.ge.0.or.iperio.eq.1) then
    call parcpt(indic)
  endif

  ! mean particle diameter and porosity value due to the deposit in each cell

  do iel = 1, ncelet
    porosi(iel) = 1.d0
    mdiam(iel) = 0.d0
  enddo

  ! Nothing more to do if no deposition at thi stage

  if (indic.eq.0) then
    deallocate(itytmp)
    return
  endif

  ! Allocate temporary arrays for the distance resolution
  allocate(distpw(ncelet))

  allocate(coefap(nfabor), coefbp(nfabor))
  allocate(q(3,ncelet))
  allocate(masflu(ncelet), depvol(ncelet))

  call field_get_val_s(icrom, crom)
  call field_get_val_s(iprpfl(iviscl), viscl)

  if (indic.gt.0) then
    call distpr(itytmp, distpw)
  endif

  !===============================================================================
  ! 2. Compute  n = Grad(DISTPW)/|Grad(DISTPW)|
  !===============================================================================

  ! Distance to wall is 0 at the wall by definition, zero flux elsewhere

  do ifac = 1, nfabor

    if (itytmp(ifac).eq.iparoi) then

      ! Dirichlet Boundary Condition for gradients
      !-------------------------------------------

      coefap(ifac) = 0.0d0
      coefbp(ifac) = 0.0d0

    else

      ! Neumann Boundary Condition for gradients
      !-----------------------------------------

      coefap(ifac) = 0.0d0
      coefbp(ifac) = 1.0d0

    endif
  enddo

  ! Compute the gradient of the distance to the wall

  inc    = 1
  iccocg = 1
  f_id   = -1

  call gradient_s                                                   &
   ( f_id   , imrgra , inc    , iccocg , nswrgy , imligy , iwarny , &
     epsrgy , climgy , extray ,                                     &
     distpw , coefap , coefbp ,                                     &
     q      )

  ! Normalisation (attention, le gradient peut etre nul, parfois)

  do iel = 1, ncel
    xnorme = max(sqrt(q(1,iel)**2+q(2,iel)**2+q(3,iel)**2),epzero)
    do isou = 1, 3
      q(isou,iel) = q(isou,iel)/xnorme
    enddo
  enddo

  ! Paralellism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(q)
  endif

  !===============================================================================
  ! 3. Compute  porosity
  !===============================================================================

  epsi = 1.d-6

  ! time step
  dt = 1.d-4

  ! volume of the deposited particles (calculated from mean deposit height)
  do iel = 1, ncelet
    depvol(iel) = 0.d0
  enddo

  indic = 0

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    depvol(iel) = depvol(iel) + parbor(ifac,ihdepm) * surfbn(ifac)
    mdiam(iel) = mdiam(iel) + parbor(ifac,ihdiam)
  enddo

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    porosi(iel) = (volume(iel) - (1.d0 - mporos) * depvol(iel)) / volume(iel)
    if (porosi(iel) .lt. mporos) indic = 1
  enddo

  ! Paralellism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call parcpt(indic)
  endif

  nn = 0
  do while (indic .gt. 0)

    do iel = 1, ncelet
      masflu(iel) = 0.d0
    enddo

    do ifac = 1 , nfac

      prod1 = 0.d0
      prod2 = 0.d0

      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)

      do ii = 1, 3
        prod1 = prod1 + q(ii,iel1) * surfac(ii,ifac)
      enddo

      do ii = 1, 3
        prod2 = prod2 - q(ii,iel2) * surfac(ii,ifac)
      enddo

      if ((porosi(iel1) .ge. mporos) .or. (prod1 .le. epsi)) then
        masflu(iel1) = masflu(iel1)
      else
        masflu(iel1) = masflu(iel1) - (porosi(iel1)-mporos)/dt * volume(iel1)
        masflu(iel2) = masflu(iel2) + (porosi(iel1)-mporos)/dt * volume(iel1)

        mdiam(iel2) = mdiam(iel1)
      endif

      if ((porosi(iel2) .ge. mporos) .or. (prod2 .le. epsi)) then
        masflu(iel2) = masflu(iel2)
      else
        masflu(iel2) = masflu(iel2) - (porosi(iel2)-mporos)/dt * volume(iel2)
        masflu(iel1) = masflu(iel1) + (porosi(iel2)-mporos)/dt * volume(iel2)

        mdiam(iel1) = mdiam(iel2)
      endif

    enddo

    indic = 0
    do iel = 1 , ncel
      porosi(iel) = porosi(iel) + (dt/volume(iel)) * masflu(iel)
    enddo

    ! Paralellism and periodicity
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(porosi)
      call synsca(mdiam)
    endif

    do iel = 1 , ncel
      if (porosi(iel) .lt. mporos) indic = 1
    enddo

    ! Paralellism and periodicity
    if (irangp.ge.0.or.iperio.eq.1) then
      call parcpt(indic)
    endif

    nn = nn + 1

    if (nn .ge. 100) then
      write(nfecra,*) '=========================================='
      write(nfecra,*) ' Error nn > 100'
      write(nfecra,*) ' Stop inside the porcel subroutine '
      write(nfecra,*) '=========================================='
      call csexit(1)
    endif

  enddo

  ! Free memory

  deallocate(masflu, depvol)
  deallocate(q)
  deallocate(coefap, coefbp)
  deallocate(distpw)
  deallocate(itytmp)

  !----
  ! End
  !----

  return

end subroutine porcel

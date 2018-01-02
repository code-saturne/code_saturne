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

!===============================================================================
! Function:
! ---------
!> \file rijech.f90
!> \brief Terms of wall echo for R_{ij}
!>        \f$var  = R_{11} \: R_{22} \: R_{33} \: R_{12} \: R_{13} \: R_{23}\f$
!>        \f$isou =  1 \:  2 \:  3 \:  4 \:  5 \:  6\f$

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     isou          passage numbero
!> \param[in]     produc        production
!> \param[in,out] smbr          work array for second member
!______________________________________________________________________________!

subroutine rijech &
 ( isou   ,                                                       &
   produc , smbr   )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          isou

double precision produc(6,ncelet)
double precision smbr(ncelet)

! Local variables

integer          iel   , ii    , jj    , kk    , mm
integer          iskm  , iski  , iskj
integer          inc   , iccocg, f_id

double precision cmu075, distxn, d2s3  , trrij , xk
double precision vnk   , vnm   , vni   , vnj
double precision deltki, deltkj, deltkm, deltij
double precision aa    , bb    , xnorme

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: produk, epsk
double precision, allocatable, dimension(:) :: w2, w3, w4, w6
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:), pointer :: w_dist
double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:), pointer :: cvara_rkm, cvara_rki, cvara_rkj
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays
allocate(produk(ncelet), epsk(ncelet))
allocate(w2(ncelet), w3(ncelet), w4(ncelet), w6(ncelet))

! Initialize variables to avoid compiler warnings

ii = 0
jj = 0
iski = 0
iskj = 0
iskm = 0

vni = 0.d0
vnj = 0.d0
vnk = 0.d0
vnm = 0.d0

! Memory

call field_get_val_s(icrom, crom)
if (isto2t.gt.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
call field_get_val_prev_s(ivarfl(ir33), cvara_r33)

deltij = 1.0d0
if(isou.gt.3) then
  deltij = 0.0d0
endif

cmu075 = cmu**0.75d0
d2s3   = 2.d0/3.d0

!===============================================================================
! 2.Calculation in the orthogonal straights cells in corresponding walls
!===============================================================================

!     The orthogonal straght is defined as -gradient of the distance
!       to the wall

allocate(grad(3,ncelet))

inc    = 1
iccocg = 1

call field_get_id("wall_distance", f_id)
call field_get_val_s(f_id, w_dist)

! Current gradient: iprev = 0
call field_gradient_scalar(f_id, 0, imrgra, inc, iccocg, grad)

! Normalization (warning, the gradient may be sometimes equal to 0)

do iel = 1 ,ncel
  xnorme = max(sqrt(grad(1,iel)**2+grad(2,iel)**2+grad(3,iel)**2),epzero)
  w2(iel) = -grad(1,iel)/xnorme
  w3(iel) = -grad(2,iel)/xnorme
  w4(iel) = -grad(3,iel)/xnorme
enddo

! Free memory
deallocate(grad)

!===============================================================================
! 3. Calculation of work variables
!===============================================================================

! ---> Production and k

do iel = 1 , ncel
  produk(iel) = 0.5d0 * (produc(1,iel)  + produc(2,iel)  + produc(3,iel))
  xk          = 0.5d0 * (cvara_r11(iel) + cvara_r22(iel) + cvara_r33(iel))
  epsk(iel)   = cvara_ep(iel)/xk
enddo

! ---> Tension rating

if     ((isou.eq.1).or.(isou.eq.4).or.(isou.eq.5)) then
  ii = 1
elseif ((isou.eq.2).or.(isou.eq.6)) then
  ii = 2
elseif  (isou.eq.3) then
  ii = 3
endif

if     ((isou.eq.3).or.(isou.eq.5).or.(isou.eq.6)) then
  jj = 3
elseif ((isou.eq.2).or.(isou.eq.4)) then
  jj = 2
elseif  (isou.eq.1) then
  jj = 1
endif

! ---> Loop for source terms construction

do iel = 1, ncel
  w6(iel) = 0.d0
enddo

do kk = 1, 3

  ! ---> Sum on m

  do mm = 1, 3

    !   --> Delta km

    if(kk.eq.mm) then
      deltkm = 1.0d0
    else
      deltkm = 0.0d0
    endif

    !  --> R km

    if     ((kk*mm).eq.1) then
      call field_get_val_prev_s(ivarfl(ir11), cvara_rkm)
      iskm = 1
    elseif ((kk*mm).eq.4) then
      call field_get_val_prev_s(ivarfl(ir22), cvara_rkm)
      iskm = 2
    elseif ((kk*mm).eq.9) then
      call field_get_val_prev_s(ivarfl(ir33), cvara_rkm)
      iskm = 3
    elseif ((kk*mm).eq.2) then
      call field_get_val_prev_s(ivarfl(ir12), cvara_rkm)
      iskm = 4
    elseif ((kk*mm).eq.6) then
      call field_get_val_prev_s(ivarfl(ir23), cvara_rkm)
      iskm = 5
    elseif ((kk*mm).eq.3) then
      call field_get_val_prev_s(ivarfl(ir13), cvara_rkm)
      iskm = 6
    endif

    !  --> Terms with R km and Phi km

    do iel = 1, ncel

      if    (kk.eq.1) then
        vnk    = w2(iel)
      elseif(kk.eq.2) then
        vnk    = w3(iel)
      elseif(kk.eq.3) then
        vnk    = w4(iel)
      endif

      if    (mm.eq.1) then
        vnm    = w2(iel)
      elseif(mm.eq.2) then
        vnm    = w3(iel)
      elseif(mm.eq.3) then
        vnm    = w4(iel)
      endif

      w6(iel) = w6(iel) + vnk*vnm*deltij*(                        &
             crijp1*cvara_rkm(iel)*epsk(iel)                      &
            -crijp2                                               &
             *crij2*(produc(iskm,iel)-d2s3*produk(iel)*deltkm) )
    enddo

  enddo

  ! ---> End of sum on m


  !  --> R ki

  if     ((kk*ii).eq.1) then
    call field_get_val_prev_s(ivarfl(ir11), cvara_rki)
    iski = 1
  elseif ((kk*ii).eq.4) then
    call field_get_val_prev_s(ivarfl(ir22), cvara_rki)
    iski = 2
  elseif ((kk*ii).eq.9) then
    call field_get_val_prev_s(ivarfl(ir33), cvara_rki)
    iski = 3
  elseif ((kk*ii).eq.2) then
    call field_get_val_prev_s(ivarfl(ir12), cvara_rki)
    iski = 4
  elseif ((kk*ii).eq.6) then
    call field_get_val_prev_s(ivarfl(ir23), cvara_rki)
    iski = 5
  elseif ((kk*ii).eq.3) then
    call field_get_val_prev_s(ivarfl(ir13), cvara_rki)
    iski = 6
  endif

  !  --> R kj

  if     ((kk*jj).eq.1) then
    call field_get_val_prev_s(ivarfl(ir11), cvara_rkj)
    iskj = 1
  elseif ((kk*jj).eq.4) then
    call field_get_val_prev_s(ivarfl(ir22), cvara_rkj)
    iskj = 2
  elseif ((kk*jj).eq.9) then
    call field_get_val_prev_s(ivarfl(ir33), cvara_rkj)
    iskj = 3
  elseif ((kk*jj).eq.2) then
    call field_get_val_prev_s(ivarfl(ir12), cvara_rkj)
    iskj = 4
  elseif ((kk*jj).eq.6) then
    call field_get_val_prev_s(ivarfl(ir23), cvara_rkj)
    iskj = 5
  elseif ((kk*jj).eq.3) then
    call field_get_val_prev_s(ivarfl(ir13), cvara_rkj)
    iskj = 6
  endif

  !   --> Delta ki

  if (kk.eq.ii) then
    deltki = 1.d0
  else
    deltki = 0.d0
  endif

  !   --> Delta kj

  if (kk.eq.jj) then
    deltkj = 1.d0
  else
    deltkj = 0.d0
  endif

  do iel = 1, ncel

      if    (kk.eq.1) then
        vnk    = w2(iel)
      elseif(kk.eq.2) then
        vnk    = w3(iel)
      elseif(kk.eq.3) then
        vnk    = w4(iel)
      endif

      if    (ii.eq.1) then
        vni    = w2(iel)
      elseif(ii.eq.2) then
        vni    = w3(iel)
      elseif(ii.eq.3) then
        vni    = w4(iel)
      endif

      if    (jj.eq.1) then
        vnj    = w2(iel)
      elseif(jj.eq.2) then
        vnj    = w3(iel)
      elseif(jj.eq.3) then
        vnj    = w4(iel)
      endif

    w6(iel) = w6(iel) + 1.5d0*vnk*(                               &
    -crijp1*(cvara_rki(iel)*vnj+cvara_rkj(iel)*vni)*epsk(iel)     &
    +crijp2                                                       &
     *crij2*((produc(iski,iel)-d2s3*produk(iel)*deltki)*vnj       &
            +(produc(iskj,iel)-d2s3*produk(iel)*deltkj)*vni) )

  enddo

enddo

! ---> Distance to the wall and amortization function: W3
!      For each calculation mode: same code, test
!      Apart from the loop

do iel = 1, ncel
  distxn =  max(w_dist(iel),epzero)
  trrij  = 0.5d0 * (cvara_r11(iel) + cvara_r22(iel) + cvara_r33(iel))
  aa = 1.d0
  bb = cmu075*trrij**1.5d0/(xkappa*cvara_ep(iel)*distxn)
  w3(iel) = min(aa, bb)
enddo

! ---> Increment of source term

do iel = 1, ncel
  smbr(iel) = smbr(iel) + cromo(iel)*volume(iel)*w6(iel)*w3(iel)
enddo

! Allocate temporary arrays
deallocate(produk, epsk)
deallocate(w2, w3, w4, w6)

return
end subroutine rijech

!===============================================================================
! Function:
! ---------
!> \file rijech.f90
!> \brief Terms of wall echo for \f$ R_{ij} \f$
!>        \f$var  = R_{11} \: R_{22} \: R_{33} \: R_{12} \: R_{13} \: R_{23}\f$
!>        \f$comp =  1 \:  2 \:  3 \:  4 \:  5 \:  6\f$

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! ARGUMENTS
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     produc        production
!> \param[in,out] smbr          work array for second member
!______________________________________________________________________________!

subroutine rijech2 &
 (produc , smbr   )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          isou

double precision produc(6,ncelet)
double precision smbr(6,ncelet)

! Local variables

integer          iel   , ii    , jj    , kk    , mm
integer          iskm  , iski  , iskj
integer          inc   , iccocg, f_id

double precision cmu075, distxn, d2s3  , trrij , xk
double precision vnk   , vnm   , vni   , vnj
double precision deltki, deltkj, deltkm, deltij(6)
double precision bb    , xnorme

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: produk, epsk
double precision, allocatable, dimension(:) :: w6
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:), pointer :: w_dist
double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:,:), pointer :: cvara_var
double precision, allocatable, dimension(:,:,:) :: cvara_rij
!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays
allocate(produk(ncelet), epsk(ncelet))
allocate(w6(ncelet))
allocate(grad(3,ncelet))
allocate(cvara_rij(3,3,ncel))
! Initialize variables to avoid compiler warnings

ii = 0
jj = 0
iski = 0
iskj = 0
iskm = 0

vni = 0.d0
vnj = 0.d0
vnk = 0.d0
vnm = 0.d0

! Memory

call field_get_val_s(icrom, crom)
if (isto2t.gt.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_prev_v(ivarfl(irij), cvara_var)

do iel = 1, ncel
  cvara_rij(1,1,iel) = cvara_var(1,iel)
  cvara_rij(2,2,iel) = cvara_var(2,iel)
  cvara_rij(3,3,iel) = cvara_var(3,iel)
  cvara_rij(1,2,iel) = cvara_var(4,iel)
  cvara_rij(2,3,iel) = cvara_var(5,iel)
  cvara_rij(1,3,iel) = cvara_var(6,iel)
  cvara_rij(2,1,iel) = cvara_var(4,iel)
  cvara_rij(3,2,iel) = cvara_var(5,iel)
  cvara_rij(3,1,iel) = cvara_var(6,iel)
enddo

do isou = 1, 6
  deltij(isou) = 1.0d0
  if (isou.gt.3) then
    deltij(isou) = 0.0d0
  endif
enddo

cmu075 = cmu**0.75d0
d2s3   = 2.d0/3.d0

!===============================================================================
! 2.Calculation in the orthogonal straights cells in corresponding walls
!===============================================================================

!     The orthogonal straght is defined as -gradient of the distance
!       to the wall

!       Calculation of gradient

inc    = 1
iccocg = 1

call field_get_id("wall_distance", f_id)

call field_get_val_s(f_id, w_dist)

! Current gradient: iprev = 0
call field_gradient_scalar(f_id, 0, imrgra, inc, iccocg, grad)

! Normalization (warning, the gradient may be sometimes equal to 0)
do iel = 1 ,ncel
  xnorme = max(sqrt(grad(1,iel)**2+grad(2,iel)**2+grad(3,iel)**2), epzero)
  grad(1, iel) = -grad(1,iel)/xnorme
  grad(2, iel) = -grad(2,iel)/xnorme
  grad(3, iel) = -grad(3,iel)/xnorme
enddo

!===============================================================================
! 3. Calculation of work variables
!===============================================================================

! ---> Production and k

do iel = 1 , ncel
  produk(iel) = 0.5d0 * (produc(1,iel)  + produc(2,iel)  + produc(3,iel))
  xk          = 0.5d0 * (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel))
  epsk(iel)   = cvara_ep(iel)/xk
enddo

! ---> Tension rating
do isou = 1, 6
  if     ((isou.eq.1).or.(isou.eq.4).or.(isou.eq.5)) then
    ii = 1
  elseif ((isou.eq.2).or.(isou.eq.6)) then
    ii = 2
  elseif  (isou.eq.3) then
    ii = 3
  endif

  if     ((isou.eq.3).or.(isou.eq.5).or.(isou.eq.6)) then
    jj = 3
  elseif ((isou.eq.2).or.(isou.eq.4)) then
    jj = 2
  elseif  (isou.eq.1) then
    jj = 1
  endif

  ! ---> Loop for source terms construction

  do iel = 1, ncel
    w6(iel) = 0.d0
  enddo

  do kk = 1, 3

    ! ---> Sum on m

    do mm = 1, 3

      !   --> Delta km

      if(kk.eq.mm) then
        deltkm = 1.0d0
      else
        deltkm = 0.0d0
      endif

      !  --> R km

      if     ((kk*mm).eq.1) then
        iskm = 1
      elseif ((kk*mm).eq.4) then
        iskm = 2
      elseif ((kk*mm).eq.9) then
        iskm = 3
      elseif ((kk*mm).eq.2) then
        iskm = 4
      elseif ((kk*mm).eq.6) then
        iskm = 5
      elseif ((kk*mm).eq.3) then
        iskm = 6
      endif

      !  --> Terms with R km and Phi km

      do iel = 1, ncel
        vnk    = grad(kk, iel)
        vnm    = grad(mm, iel)

        w6(iel) = w6(iel) + vnk*vnm*deltij(isou)*(                        &
               crijp1*cvara_rij(kk,mm,iel)*epsk(iel)                      &
              -crijp2                                               &
               *crij2*(produc(iskm,iel)-d2s3*produk(iel)*deltkm) )
      enddo

    enddo

    ! ---> End of sum on m


    !  --> R ki

    if     ((kk*ii).eq.1) then
      iski = 1
    elseif ((kk*ii).eq.4) then
      iski = 2
    elseif ((kk*ii).eq.9) then
      iski = 3
    elseif ((kk*ii).eq.2) then
      iski = 4
    elseif ((kk*ii).eq.6) then
      iski = 5
    elseif ((kk*ii).eq.3) then
      iski = 6
    endif

    !  --> R kj

    if     ((kk*jj).eq.1) then
      iskj = 1
    elseif ((kk*jj).eq.4) then
      iskj = 2
    elseif ((kk*jj).eq.9) then
      iskj = 3
    elseif ((kk*jj).eq.2) then
      iskj = 4
    elseif ((kk*jj).eq.6) then
      iskj = 5
    elseif ((kk*jj).eq.3) then
      iskj = 6
    endif

    !   --> Delta ki

    if (kk.eq.ii) then
      deltki = 1.d0
    else
      deltki = 0.d0
    endif

    !   --> Delta kj

    if (kk.eq.jj) then
      deltkj = 1.d0
    else
      deltkj = 0.d0
    endif

    do iel = 1, ncel

        vnk    = grad(kk, iel)
        vni    = grad(ii, iel)
        vnj    = grad(jj, iel)

      w6(iel) = w6(iel) + 1.5d0*vnk*(                               &
      -crijp1*(cvara_rij(kk,ii,iel)*vnj+cvara_rij(kk,jj,iel)*vni)*epsk(iel)     &
      +crijp2                                                       &
       *crij2*((produc(iski,iel)-d2s3*produk(iel)*deltki)*vnj       &
              +(produc(iskj,iel)-d2s3*produk(iel)*deltkj)*vni) )

    enddo

  enddo
enddo

! ---> Distance to the wall and amortization function: W3
!      For each calculation mode: same code, test
!      Apart from the loop

do iel = 1, ncel
  distxn =  max(w_dist(iel), epzero)!FIXME
  trrij  = 0.5d0 * (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel))
  bb = cmu075*trrij**1.5d0/(xkappa*cvara_ep(iel)*distxn)
  bb = min(bb, 1.d0)

  ! ---> Increment of source term
  do isou = 1, 6
    smbr(isou,iel) = smbr(isou,iel) + cromo(iel)*cell_f_vol(iel)*w6(iel)*bb
  enddo
enddo

! Free memory
deallocate(produk, epsk)
deallocate(w6)
deallocate(grad)

return
end subroutine rijech2

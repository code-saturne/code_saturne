!-------------------------------------------------------------------------------

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


!===============================================================================
! Function :
! ----------
!> \file resalp.f90
!> \brief Solving the equation on alpha in the framwork of the Rij-EBRSM model.
!>        (written from the equation of \f$ \overline{f})\f$

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           nom           role
!______________________________________________________________________________!
!> \param[in]     nvar          total numver of variables
!> \param[in]     dt            step time
!> \param[in]     rtp, rtpa     calculation variable in cell centers
!>                              (current or previous instant)
!______________________________________________________________________________!

subroutine resalp &
 ( nvar   ,                                                       &
   rtp    , rtpa   )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstnum
use cstphy
use pointe
use period
use parall
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar

double precision rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)

! Local variables

integer          ivar  , iel
integer          ipcvlo
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp, ipp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          imucpp, idftnp, iswdyp
integer          icvflb
integer          ivoid(1)
double precision blencp, epsilp, epsrsp, epsrgp, climgp, extrap, relaxp
double precision thetv , thetap
double precision d1s4, d3s2, d1s2
double precision xk, xnu, l2
double precision xllke, xllkmg, xlldrb

double precision rvoid(1)

character(len=80) :: label
double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: cvara_ep, cvara_al
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: viscl, visct

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

allocate(smbr(ncelet), rovsdt(ncelet), w1(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(dpvar(ncelet))

call field_get_val_s(icrom, crom)
call field_get_val_s(iprpfl(iviscl), viscl)
call field_get_val_s(iprpfl(ivisct), visct)

call field_get_val_prev_s(ivarfl(iep), cvara_ep)
call field_get_val_prev_s(ivarfl(ial), cvara_al)

call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
call field_get_val_prev_s(ivarfl(ir33), cvara_r33)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

d1s2 = 1.d0/2.d0
d1s4 = 1.d0/4.d0
d3s2 = 3.d0/2.d0

!  test on alpha which must not be above 1
if (iwarni(ial).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Resolving the equation of alpha
!===============================================================================

ivar = ial
ipp    = ipprtp(ivar)

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

if(iwarni(ivar).ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1100) label
endif

thetv  = thetav(ivar)

ipcvlo = ipproc(iviscl)
if(isto2t.gt.0) then
  if(iviext.gt.0) then
    ipcvlo = ipproc(ivisla)
  endif
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 2.2 Source term of alpha
!     \f$ smbr = \dfrac{1}{L^2 (\alpha)} - \dfrac{1}{L^2}\f$
!  In fact there is a mark "-" because the solved equation is
!    \f$-\div{\grad {alpha}} = smbr \f$
!===============================================================================

! ---> Matrix

if (isto2t.gt.0) then
  thetap = thetv
else
  thetap = 1.d0
endif

!FIXME the source term extrapolation is not well done!!!!
do iel=1,ncel

  xk = d1s2*(cvara_r11(iel)+cvara_r22(iel)+cvara_r33(iel))
  xnu  = viscl(iel)/crom(iel)

  ! Integral length scale
  xllke = xk**d3s2/cvara_ep(iel)

  ! Kolmogorov length scale
  xllkmg = xceta*(xnu**3/cvara_ep(iel))**d1s4

  ! Durbin length scale
  xlldrb = xcl*max(xllke,xllkmg)

  l2      = xlldrb**2

  ! Explicit term
  smbr(iel) = volume(iel)*(1.d0 -cvara_al(iel)) / l2

  ! Implicit term
  rovsdt(iel) = (rovsdt(iel) + volume(iel)*thetap) / l2

enddo

! Calculation of viscf and viscb for codits

do iel = 1, ncel
  w1(iel) = 1.d0
enddo

call viscfa                                                       &
!==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

!===============================================================================
! 2.3 Effective resolution of the equation of alpha
!===============================================================================

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imucpp = 0
idftnp = idften(ivar)
iswdyp = iswdyn(ivar)
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
! all boundary convective flux with upwind
icvflb = 0

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 3. Clipping
!===============================================================================

call clpalp(ncelet, ncel, nvar, rtp)
!==========

! Free memory
deallocate(smbr, rovsdt, w1)
deallocate(viscf, viscb)
deallocate(dpvar)

!--------
! Formats
!--------

 1000    format(/,                                                &
'   ** Solving alpha                                          ',/,&
'      -----------------------------------------------        ',/)
 1100    format(/,'           Solving the variable ',A8,/)

!----
! End
!----

return

end

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file rafsck.f90
!>
!> \brief This subroutine is part of the FSCK radiation model
!>
!> Determination of the radiation coeffcients of the FSCK model as well as the
!> corresponding weights.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     pco2          CO2 volume fraction
!> \param[in]     ph2o          H2O volume fraction
!> \param[in]     fv            Soot volume fraction
!> \param[in]     teloc         Gas temperature
!> \param[out]    kloc          Radiation coeffcient of the i different gases
!> \param[out]    aloc          Weights of the i different gases in cells
!> \param[out]    alocbo        Weights of the i different gases at boundaries
!_______________________________________________________________________________

subroutine rafsck &
(pco2 , ph2o , fv, teloc, kloc, aloc, alocbo)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field
!===============================================================================

implicit none

! Arguments

double precision      pco2(ncelet), ph2o(ncelet), fv(ncelet), teloc(ncelet)
double precision      kloc(ncelet,nwsgg), aloc(ncelet,nwsgg)
double precision      alocbo(nfabor,nwsgg)

! Local variables

character(len=256) :: pathdatadir

integer               ipass
integer               cco2, ch2o, it, itrad, ig, iwvnb
integer               its, j, m, n, i
integer               iel
integer               interp_method, k, ifac

integer, parameter :: ng = 100, nt = 38, nconc = 5, nband = 450
integer, parameter :: maxit = 1000
integer, parameter :: imlinear = 0, im3dspline = 7, im1dspline = 2
integer, parameter :: im2dspline = 6

double precision      wvnb, dwvnb
double precision      x1, x2, xl, xm
double precision      tref, pco2ref, ph2oref, sum1, sum2, kp, kp1, kp2, kpt4dv

double precision, parameter :: eps = 3.0d-14
double precision, parameter :: x_kg(nconc)=(/1.d-4,0.25d0,0.50d0,0.75d0,1.00d0/)

logical, allocatable, dimension(:)          :: unfinished
double precision, allocatable, dimension(:) :: gi, tt, kpco2, kph2o, wv, dwv
double precision, allocatable, dimension(:,:,:,:,:) :: kmfs
double precision, allocatable, dimension(:) :: gq, p1, p2, p3, pp, z, z1, arth
double precision, allocatable, dimension(:) ::  gFSKref, kFSKref, kFSK, gFSK, gg1
double precision, allocatable, dimension(:) ::  kg1, as, ag, aw, kloctmp

double precision, dimension(:), pointer     ::  tpfsck

data ipass /0/

save ipass, gi, kmfs, tt, kpco2, kph2o, wv, dwv
save gq

!===============================================================================
! 0 - Initialization
!===============================================================================

allocate(gfskref(ng),kfskref(ng),gfsk(ng),kfsk(ng),gg1(ng),kg1(ng),as(ng))
allocate(ag(nwsgg),aw(nwsgg),kloctmp(nwsgg))

call field_get_val_s(itempb, tpfsck)

do i = 1, nwsgg
  kloctmp(i) = 0.d0
enddo

!===============================================================================
!  1 - Read the data base files
!===============================================================================

ipass = ipass + 1
if (ipass.eq.1) then

  call csdatadir(len(pathdatadir), pathdatadir)

  allocate(gi(ng),tt(nt),kpco2(nt),kph2o(nt),wv(nband),dwv(nband))
  allocate(kmfs(nconc,nconc,nt,nt,ng))

  ! Read k-distributions
  open(unit=100, file=trim(pathdatadir)// '/data/thch/dp_radiat_MFS')
  do cco2 = 1, nconc
    do ch2o = 1, nconc
      do it = 1, nt
        do itrad = 1, nt
          read(100,140)
          read(100,140)
          do ig = 1, ng
            read(100,140) gi(ig), kmfs(ch2o,cco2,it,itrad,ig)
          enddo
        enddo
      enddo
    enddo
  enddo
  close(100)

  ! Read the Planck coefficients
  open(unit=101, file=trim(pathdatadir)//'/data/thch/dp_radiat_Planck_CO2')
  do it = 1, nt
    read(101,140) tt(it), kpco2(it)
  enddo
  close(101)
  open(unit=102, file=trim(pathdatadir)//'/data/thch/dp_radiat_Planck_H2O')
  do it = 1, nt
    read(102,140) tt(it), kph2o(it)
  enddo
  close(102)

  ! Read the wavelength intervall
  open(unit=103, file=trim(pathdatadir)//'/data/thch/dp_radiat_wave')
  do iwvnb = 1, nband
    read(103,140) wv(iwvnb), dwv(iwvnb)
  enddo
  close(103)

!===============================================================================
!  2 - Quadrature Gaussian
!===============================================================================
! Allocation
  allocate(gq(nwsgg))
  allocate(p1((size(gq)+1)/2),p2((size(gq)+1)/2),p3((size(gq)+1)/2))
  allocate(pp((size(gq)+1)/2),z((size(gq)+1)/2), &
           z1((size(gq)+1)/2),arth((size(gq)+1)/2))
  allocate(unfinished((size(gq)+1)/2))
! Boundaries
  x1 = 0.d0
  x2 = 1.d0
! Given the lower and upper limits of integration x1 and x2, this routine
! returns arrays x[0..n-1]
! and w[0..n-1] of length n, containing the abscissas and weights of the
! Gauss-Legendre n-point
! quadrature formula. The parameter EPS is the relative precision.
! Note that internal computations are done in double precision.
  n = nwsgg
! The roots are symmetric in the interval, so we only have to find half of them.
  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)

  do i = 1, m
    arth(i) = i
  enddo
! initial aproximation of the roots
  z = cos(acos(-1.0d0)*(arth-0.25d0)/(n+0.5d0))

  unfinished =.true.

! Newton's method carried out simultaneously on the roots
  do its = 1, maxit
    where (unfinished)
      p1 = 1.0d0
      p2 = 0.0d0
    end where

    do j = 1, n
      where (unfinished)
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
      end where
    end do

  ! p1 now contains the desired Legendre polynomials.
  ! We next compute pp, its derivative, by a standard relation involving
  ! also p2, the polynomial of one lower order.

    where (unfinished)
      pp = n*(z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z = z1-p1/pp
      unfinished = (abs(z-z1) .gt. eps)
    end where
    if(.not. any(unfinished)) exit
  end do

  if (its == maxit+1) write(nfecra,*)"Maximum number of iterations during GAULEG"

! Scale the root to the desired interval, and put in its symmetric counterpart.
  gq(1:m) = xm-xl*z
  gq(n:n-m+1:-1) = xm+xl*z

! Compute the weight and its symmetric counterpart.
  wq(1:m) = 2.0d0*xl/((1.0d0-z**2.d0)*pp**2.d0)
  wq(n:n-m+1:-1) = wq(1:m)

  deallocate(p1, p2, p3, pp, z, z1, arth, unfinished)

endif

!===============================================================================
!  3 - Calculate the reference state
!===============================================================================
pco2ref = 0.d0
ph2oref = 0.d0
tref    = 0.d0
sum1    = 0.d0
sum2    = 0.d0

do iel = 1, ncel

  ! Calculation of pco2ref and ph2oref
  pco2ref = pco2ref + pco2(iel)*volume(iel)
  ph2oref = ph2oref + ph2o(iel)*volume(iel)

  ! Calculation of tref
  ! Interplolation du coefficient de Planck pour  mélange
  if (teloc(iel).le.tt(1)) then
    kp = pco2(iel)*kpco2(1) + ph2o(iel)*kph2o(1)
    elseif (teloc(iel).ge.tt(nt))then
    kp = pco2(iel)*kpco2(nt) + ph2o(iel)*kpco2(nt)
  else
    do it = 1, nt-1
      if ((teloc(iel).ge.tt(it)).and.(teloc(iel).lt.tt(it+1))) then
        kp1 = pco2(iel)*kpco2(it)   + ph2o(iel)*kph2o(it)
        kp2 = pco2(iel)*kpco2(it+1) + ph2o(iel)*kph2o(it+1)
        kp  = (kp2-kp1)/(tt(it+1)-tt(it))*(teloc(iel)-tt(it)) + kp1
        exit
      endif
    enddo
  endif

  kpt4dv = kp*teloc(iel)*teloc(iel)*teloc(iel)*teloc(iel)*volume(iel)
  sum1   = sum1 + kpt4dv*teloc(iel)
  sum2   = sum2 + kpt4dv

enddo

if (irangp.ge.0) then
  call parsom (pco2ref)
  call parsom (ph2oref)
  call parsom (sum1)
  call parsom (sum2)
endif

pco2ref = pco2ref/voltot
ph2oref = ph2oref/voltot
tref    = sum1/sum2

!===============================================================================
!  4 - Main program
!===============================================================================

! Determination of the k-distribution at the reference state
interp_method = imlinear
call interpolation4d(tref,tref,pco2ref,ph2oref,gi,kmfs,interp_method,gfskref,  &
                     kfskref)
kfskref = kfskref*100.d0  ! [m^-1]

do iel = 1, ncel
  ! Determination of the local absorbtion coeffcient
  kfsk = 0.d0
  gfsk = 0.d0
  call interpolation4d(tref,teloc(iel),pco2(iel),ph2o(iel),gi,kmfs,            &
                       interp_method,gfsk,kfsk)
  kfsk = kfsk*100.d0  ! [m^-1]
  call simple_interpg(ng,gfsk,kfsk,nwsgg,gq,kloctmp)
  do i = 1, nwsgg
    kloc(iel,i) = kloctmp(i)
  enddo
  ! Determination of the local weights
  kfsk = 0.d0
  gfsk = 0.d0
  call interpolation4d(teloc(iel),tref,pco2ref,ph2oref,gi,kmfs,interp_method,  &
                       gg1,kg1)
  kg1 = kg1*100.d0 ! [m^-1]
  call simple_interpg(ng,kg1,gg1,ng,kfskref,gfsk)
  as(1) = (gfsk(2)-gfsk(1))/(gfskref(2)-gfskref(1)+1.d-15)
  as(ng) = (gfsk(ng)-gfsk(ng-1))/(gfskref(ng)-gfskref(ng-1)+1.d-15)
  do k = 2, ng-1
    as(k) = (gfsk(k+1)-gfsk(k-1))/(gfskref(k+1)-gfskref(k-1)+1.d-15)
  enddo
  call simple_interpg(ng,gfskref,as,nwsgg,gq,ag)
  do i = 1, nwsgg
    aloc(iel,i) = ag(i)
  enddo
enddo

do ifac = 1, nfabor
  iel = ifabor(ifac)
  call interpolation4d(tpfsck(ifac),tref,pco2ref,ph2oref,gi,kmfs,interp_method,&
                       gg1,kg1)
  kg1 = kg1*100.d0
  call simple_interpg(ng,kg1,gg1,ng,kfskref,gfsk)
  as(1)  = (gfsk(2)-gfsk(1))/(gfskref(2)-gfskref(1)+1.d-15)
  as(ng) = (gfsk(ng)-gfsk(ng-1))/(gfskref(ng)-gfskref(ng-1)+1.d-15)
  do k = 2, ng-1
    as(k) = (gfsk(k+1)-gfsk(k-1))/(gfskref(k+1)-gfskref(k-1)+1.d-15)
  enddo
  call simple_interpg(ng,gfskref,as,nwsgg,gq,aw)
  do i = 1, nwsgg
    alocbo(ifac,i) = aw(i)
  enddo
enddo

!Deallocation
if (ntcabs.eq.ntmabs) then
  deallocate(gi,tt,kpco2,kph2o,wv,dwv,kmfs)
  deallocate(gq)
endif

deallocate(kloctmp)

!Formats
140    format (2(E11.4,1X))

return

contains

!===============================================================================
!  Interpolation routines
!===============================================================================

!===============================================================================
!  interpolation4d
!===============================================================================
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     trad          Reference temperature
!> \param[in]     t             Local temperature
!> \param[in]     xco2          Reference CO2 volume fraction
!> \param[in]     xh2o          Reference H2O volume fraction
!> \param[in]     gi
!> \param[in]     kmfs
!> \param[in]     interp_method Interpolation method
!> \param[out]    gdb
!> \param[out]    kdb
!_______________________________________________________________________________

subroutine interpolation4d(trad,t,xco2,xh2o,gi,kmfs,interp_method,gdb,kdb)

implicit none

! Arguments
integer           interp_method

double precision  trad, t, xco2, xh2o

double precision  kmfs(nconc,nconc,nt,nt,ng)
double precision  gi(ng), kdb(ng), gdb(ng)

! Local variables
integer           itx(4,4), itrada(4), ita(4), ico2a(4), ih2oa(4)
integer           i, ico2, ih2o, it, itrad, ig, nix, nit

double precision  wx, wt

double precision, allocatable, dimension(:,:,:,:,:) :: karray
double precision, allocatable, dimension(:,:,:,:)   :: kint1
double precision, allocatable, dimension(:,:,:)     :: kint2
double precision, allocatable, dimension(:,:)       :: kint3
double precision, allocatable, dimension(:)         :: b, c, d, kg_t2, kg_x2

! Table allocation
allocate(karray(4,4,4,4,ng), kint1(4,4,4,ng), kint2(4,4,ng), kint3(4,ng))
allocate(b(4), c(4), d(4), kg_t2(4), kg_x2(4))

! Détermination des positions en x et T
! dans le NB database pour l'interpolation.
gdb = gi

! Détermination des positions en x et T
! dans le NB database pour l'interpolation.
!! nbre de points d'interpolation suivant t & x: 2 linear ou 4 spline
call gridposnbsg1(trad,t,xco2,xh2o,interp_method,itx)

if (interp_method == im2dspline) then !!!! spline sur x
  nix = 4
else
  nix = 2
endif

if (interp_method/= 0) then  !!! spline sur t
  nit = 4
else
  nit = 2
endif

! Attribution des indices des points d'interpolation suivant T & x

itrada = 0
ita = 0
ico2a = 0
ih2oa = 0

do i = 1, nit
  itrada(i) = itx(1,i)
enddo

do i = 1, nit
  ita(i) = itx(2,i)
enddo

do i = 1, nix
  ico2a(i) = itx(3,i)
enddo

do i = 1, nix
  ih2oa(i) = itx(4,i)
enddo

do ih2o = 1, nix
  do ico2 = 1, nix
    do it = 1, nit
      do itrad = 1, nit
        do ig = 1, ng
          karray(ih2o,ico2,it,itrad,ig) = &
          kmfs(ih2oa(ih2o),ico2a(ico2),ita(it),itrada(itrad),ig)
        enddo
      enddo
    enddo
  enddo
enddo

! INTERPOLATION SUR XH2O
if (interp_method == im2dspline) then  ! spline sur x
  do ico2 = 1, nix
    do it = 1, nit
      do itrad = 1, nit
        do ig = 1, ng
          kg_x2(1) = x_kg(ih2oa(1)); kg_x2(2) = x_kg(ih2oa(2))
          kg_x2(3) = x_kg(ih2oa(3)); kg_x2(4) = x_kg(ih2oa(4))
          call splmi(4,kg_x2,karray(1:4,ico2,it,itrad,ig),b,c,d)
          kint1(ico2,it,itrad,ig) =  &
            seval(4,xh2o,kg_x2,karray(1:4,ico2,it,itrad,ig),b,c,d)
        enddo
      enddo
    enddo
  enddo
else
  wx = (xh2o-x_kg(ih2oa(1)))/(x_kg(ih2oa(2))-x_kg(ih2oa(1)))
  do ico2 = 1, nix
    do it = 1, nit
      do itrad = 1, nit
        do ig = 1, ng
          kint1(ico2,it,itrad,ig) = &
            wx*karray(2,ico2,it,itrad,ig)+(1.d0-wx)*karray(1,ico2,it,itrad,ig)
        enddo
      enddo
    enddo
  enddo
endif

! INTERPOLATION SUR XCO2
if (interp_method == im2dspline) then  ! spline sur x
  do it = 1, nit
    do itrad = 1, nit
      do ig = 1, ng
        kg_x2(1) = x_kg(ico2a(1)); kg_x2(2) = x_kg(ico2a(2))
        kg_x2(3) = x_kg(ico2a(3)); kg_x2(4) = x_kg(ico2a(4))
        call splmi(4,kg_x2,kint1(1:4,it,itrad,ig),b,c,d)
        kint2(it,itrad,ig) = seval(4,xco2,kg_x2,kint1(1:4,it,itrad,ig),b,c,d)
      enddo
    enddo
  enddo
else
  wx = (xco2-x_kg(ico2a(1)))/(x_kg(ico2a(2))-x_kg(ico2a(1)))
  do it = 1, nit
    do itrad = 1, nit
      do ig = 1, ng
        kint2(it,itrad,ig) = wx*kint1(2,it,itrad,ig)+(1.d0-wx)*kint1(1,it,itrad,ig)
      enddo
    enddo
  enddo
endif

! INTERPOLATION SUR T
if (interp_method/= 0) then !spline sur t
  do itrad = 1, nit
    do ig = 1, ng
      kg_t2(1) = tt(ita(1)); kg_t2(2) = tt(ita(2))
      kg_t2(3) = tt(ita(3)); kg_t2(4) = tt(ita(4))
      call splmi(4,kg_t2,kint2(1:4,itrad,ig),b,c,d)
      kint3(itrad,ig) = seval(4,t,kg_t2,kint2(1:4,itrad,ig),b,c,d)
    enddo
  enddo
else
  wt = (t-tt(ita(1)))/(tt(ita(2))-tt(ita(1)))
  do itrad = 1, nit
    do ig = 1, ng
      kint3(itrad,ig) = wt*kint2(2,itrad,ig)+(1.d0-wt)*kint2(1,itrad,ig)
    enddo
  enddo
endif

! INTERPOLATION SUR Trad
if (interp_method/= 0) then !spline sur trad
  do ig = 1, ng
    kg_t2(1) = tt(itrada(1)); kg_t2(2) = tt(itrada(2))
    kg_t2(3) = tt(itrada(3)); kg_t2(4) = tt(itrada(4))
    call splmi(4,kg_t2,kint3(1:4,ig),b,c,d)
    kdb(ig) = seval(4,trad,kg_t2,kint3(1:4,ig),b,c,d)
  enddo
else
  wt = (trad-tt(itrada(1)))/(tt(itrada(2))-tt(itrada(1)))
  do ig = 1, ng
    kdb(ig) = wt*kint3(2,ig)+(1.d0-wt)*kint3(1,ig)
  enddo
endif

! Free memory
deallocate(karray,kint1,kint2,kint3,b,c,d,kg_t2,kg_x2)

return
end subroutine interpolation4d

!=============================================================================
! gridposnbsg1
!=============================================================================

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     trad          Reference temperature
!> \param[in]     t             Local Temperature
!> \param[in]     xco2          CO2 volume fraction
!> \param[in]     xh2o          H2O volume fraction
!> \param[in]     interp_method Interpolation method
!> \param[in]     itx           itx(1,:): TPlanck; iTx(2,:): Tloc;
!>                              iTx(3,:): xCO2; iTx(4,:) = xH2O
!_______________________________________________________________________________

subroutine gridposnbsg1(trad,t,xco2,xh2o,interp_method,itx)

!===============================================================================

implicit none

! Arguments
integer               interp_method
integer               itx(4,4) !
double precision      trad, t, xco2, xh2o

! Local variable
integer               itrad(4), ita(4), ico2a(4), ih2oa(4), ip, it, ix, i, j

!===============================================================================

itrad = 0
ita   = 0
ico2a = 0
ih2oa = 0

! Interpolation mole fraction
! Deterination de = xCO2 dans la base de données la plus proche

i  = 1
j  = nconc
ip = (i+j)/2
do while(j-i.gt.1)
  if (xco2 .lt. x_kg(ip)) then
    j = ip
  else
    i = ip
  endif
  ip = (i+j)/2
enddo

if (interp_method == im2dspline) then   !  spline on x
  if ( (i.gt.1).and.(i.lt.nconc-1) ) then
    ico2a(1) = i-1
    ico2a(2) = i
    ico2a(3) = i+1
    ico2a(4) = i+2
    elseif (i == 1) then
    ico2a = (/1,2,3,4/)
    elseif (i == nconc-1) then
    ico2a(1) = nconc-3
    ico2a(2) = nconc-2
    ico2a(3) = nconc-1
    ico2a(4) = nconc
  else
    write(nfecra,*) 'x grid failure,i =', i
  endif
else
  if (i.lt.1) then
    i = 1
    elseif(i.gt.nconc-1) then
    i = nconc-1
  endif
  ico2a(1) = i
  ico2a(2) = i+1
endif

! Interpolation mole fraction
! Deterination de = xH20 dans la base de données la plus proche

i  = 1
j  = nconc
ip = (i+j)/2
do while(j-i.gt.1)
  if (xh2o .lt. x_kg(ip)) then
    j = ip
  else
    i = ip
  endif
  ip = (i+j)/2
enddo

if (interp_method == im2dspline) then   !  spline on x
  if ( (i.gt.1).and.(i.lt.nconc-1) ) then
    ih2oa(1) = i-1
    ih2oa(2) = i
    ih2oa(3) = i+1
    ih2oa(4) = i+2
    elseif (i == 1) then
    ih2oa = (/1,2,3,4/)
    elseif (i == nconc-1) then
    ih2oa(1) = nconc-3
    ih2oa(2) = nconc-2
    ih2oa(3) = nconc-1
    ih2oa(4) = nconc
  else
    write(nfecra,*) 'x grid failure,i =', i
  endif
else
  if (i.lt.1) then
    i = 1;
    elseif(i.gt.nconc-1) then
    i = nconc-1
  endif
  ih2oa(1) = i
  ih2oa(2) = i+1
endif

! Interpolation temperature
! Determination de la température locale dans la base de données la plus proche

i  = 1
j  = nt
ip = (i+j)/2
do while(j-i.gt.1)
  if (t .lt.tt(ip)) then
    j = ip
  else
    i = ip
  endif
  ip = (i+j)/2
enddo

if (interp_method/= 0) then       !spline in t
  if ( (i.gt.1).and.(i.lt.nt-1) ) then
    ita(1) = i-1
    ita(2) = i
    ita(3) = i+1
    ita(4) = i+2
    elseif (i == 1) then
    ita = (/1,2,3,4/)
    elseif (i == nt) then
    ita(1) = nt-3
    ita(2) = nt-2
    ita(3) = nt-1
    ita(4) = nt
  else
    write(nfecra,*) 't grid failure,i =', i
  endif
else ! linear in t
  if (i.lt.1) then
    i = 1
    elseif(i.gt.nt-1) then
    i = nt-1
  endif
  ita(1) = i
  ita(2) = i+1
endif

! Deterination de la température de pLanck dans la base de données la plus proche
! inférieure à la température locale

i  = 1
j  = nt
ip = (i+j)/2
do while(j-i.gt.1)
  if (trad .lt. tt(ip)) then
    j = ip
  else
    i = ip
  endif
  ip = (i+j)/2
enddo

if (interp_method/= 0) then       !spline in t
  if ( (i.gt.1).and.(i.lt.nt-1) ) then
    itrad(1) = i-1
    itrad(2) = i
    itrad(3) = i+1
    itrad(4) = i+2
    elseif (i == 1) then
    itrad = (/1,2,3,4/)
    elseif (i == nt-1) then
    itrad(1) = nt-3
    itrad(2) = nt-2
    itrad(3) = nt-1
    itrad(4) = nt
  else
    write(nfecra,*) 't grid failure,i =', i
  endif
else ! linear in t
  if (i.lt.1) then
    i = 1
    elseif(i.gt.nt-1) then
    i = nt-1
  endif
  itrad(1) = i
  itrad(2) = i+1
endif

! attribution itx
itx(1,:) = itrad(:)
itx(2,:) = ita(:)
itx(3,:) = ico2a(:)
itx(4,:) = ih2oa(:)

return
end subroutine gridposnbsg1

!===============================================================================
! Function Seval
!===============================================================================
!
!  This subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
!-------------------------------------------------------------------------------
! Arguments
!_______________________________________________________________________________
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n            Number of data points
!> \param[in]     u            Abscissa at which the spline is to be evaluated
!> \param[in]     x,y          Arrays of data abscissas and ordinates
!> \param[in]     b,c,d        Arrays of spline coefficients computed by spline
!_______________________________________________________________________________

double precision function seval(n, u, x, y, b, c, d)

! Arguments

integer           n
double precision  u, x(n), y(n), b(n), c(n), d(n)

! Local variable

Integer i, j, k
double precision dx
data i/1/

if (i .ge. n) i = 1
if (u .lt. x(i)) go to 10
if (u .le. x(i+1)) go to 30

!  Binary search

10 i = 1
j = n+1
20 k = (i+j)/2
if ( u .lt. x(k) ) j = k
if ( u .ge. x(k) ) i = k
if ( j .gt. i+1 ) go to 20

!  Evaluate spline

30 dx = u - x(i)
seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function

!===============================================================================
! Function Splmi
!===============================================================================
! 1-d monotonic spline interpolation: coefficients
!  The coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
!  for a monotonically varying cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!    for  x(i) .le. x .le. x(i+1)
!    with y(i+1).ge.y(i) (all i) or y(i+1).lt.y(i) (all i)
!
!*******************************************************************************
!
! THEORY FROM 'MONOTONE PIECEWISE CUBIC INTERPOLATION',
! BY F.N. FRITSCH AND R.E. CARLSON IN SIAM J.NUMER.ANAL.,V.17,P.238
!
!*******************************************************************************
!
!    Y(I) = S(X(I))
!    B(I) = SP(X(I))
!    C(I) = SPP(X(I))/2
!    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!
!  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
!  TO EVALUATE THE SPLINE.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n           Number of data points or knots (n.ge.2)
!> \param[in]     x           Abscissa of the knots in strictly increasing order
!> \param[in]     y           Ordinates of the knots
!> \param[out]    b,c,d       Arrays of spline coefficients as defined above
!_______________________________________________________________________________

subroutine splmi (n, x, y, b, c, d)

implicit none

! Arguments
integer            n
double precision   x(n), y(n), b(n), c(n), d(n)

! Local variables
integer            nm1, i, nm2
double precision   phi, ti, delta(n), h(n), al(n), be(n)

nm1 = n-1
if (n .lt. 2) return
if (n .lt. 3) go to 100
!
! Calculate the h(i) and delta(i)
!
do 10 i = 1, nm1
  al(i) = 0.
  be(i) = 0.
  h(i) = x(i+1)-x(i)
10 delta(i) = (y(i+1)-y(i))/h(i)

! calculate first values for al and be by 3-point difference

if (delta(1).eq.0) goto 15

al(1) = ((h(1)+h(2))**2*y(2)-h(1)**2*y(3)-h(2)*(2.*h(1)+h(2))&
  *y(1))/(h(2)*(h(1)+h(2))*(y(2)-y(1)))

15 do 20 i = 2, nm1
  if (delta(i).eq.0) goto 20

  al(i) = (h(i-1)**2*y(i+1)+(h(i)**2-h(i-1)**2)*y(i)-h(i)**2*&
  y(i-1))/(h(i-1)*(h(i)+h(i-1))*(y(i+1)-y(i)))

20 continue

nm2 = n-2
do 30 i = 1, nm2
  if (delta(i).eq.0.) goto 30
  be(i) = (h(i)**2*y(i+2)+(h(i+1)**2-h(i)**2)*y(i+1)-h(i+1)**2*&
  y(i))/(h(i+1)*(h(i)+h(i+1))*(y(i+1)-y(i)))

30 continue

if (delta(n-1).eq.0.) goto 35

be(n-1) = (h(n-2)*(2.*h(n-1)+h(n-2))*y(n)-(h(n-1)+h(n-2))**2&
*y(n-1)+h(n-1)**2*y(n-2))/(h(n-2)*(h(n-1)+h(n-2))&
*(y(n)-y(n-1)))

! Correct values for al and be

35 do 40 i = 1, nm1
  if (al(i)+be(i).le.2.) goto 40
  if (2.*al(i)+be(i).le.3.) goto 40
  if(al(i)+2.*be(i).le.3.) goto 40
  phi = al(i)-(2.*al(i)+be(i)-3.)**2/(al(i)+be(i)-2.)/3.
  if(phi.ge.0.) goto 40
  ti = 3./sqrt(al(i)**2+be(i)**2)
  al(i) = ti*al(i)
  be(i) = ti*be(i)

40 continue

! Calculate spline coefficients

do 50 i = 1, nm1
  d(i) = (al(i)+be(i)-2.)*delta(i)/(h(i)**2)
  c(i) = (3.-2.*al(i)-be(i))*delta(i)/(h(i))
  b(i) = al(i)*delta(i)
  if(b(i)*delta(i).ge.0.) goto 50
  b(i) = delta(i)
  c(i) = 0.
  d(i) = 0.

50 continue
b(n) = be(n-1)*delta(n-1)
c(n) = 0.
d(n) = 0.

return
100 b(1) = (y(2)-y(1))/(x(2)-x(1))
c(1) = 0.
d(1) = 0.
b(2) = b(1)
c(2) = 0.
d(2) = 0.

return
end subroutine splmi

!===============================================================================
! Function simple_interpg
!===============================================================================
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nxy
!> \param[in]     xx
!> \param[in]     yy
!> \param[in]     ni
!> \param[in]     xi
!> \param[out]    yi
!_______________________________________________________________________________
subroutine simple_interpg(nxy,xx,yy,ni,xi,yi)

implicit none

! Arguments

integer          nxy, ni

double precision xx(nxy), yy(nxy), xi(ni), yi(ni)

! Local variables
integer          i, iq, ibgn

ibgn = 1
do iq = 1, ni
  do i = ibgn, nxy
    if (xi(iq).lt.xx(1)) then
      yi(iq) = yy(1)*xi(iq)/max(1.0d-09,xx(1))
      exit
    else if (xi(iq).gt.xx(nxy)) then
      yi(iq) = yy(nxy)
    end if
    ! interpolate
    if (abs(xi(iq)-xx(i))/(xx(i)+1.d-15).lt.1.0d-3 ) then
      yi(iq) = yy(i)
      ibgn = i
      exit
    else if (xx(i).gt.xi(iq) ) then
      yi(iq) = yy(i-1) + ( yy(i)-yy(i-1) ) * ( xi(iq)-xx(i-1) ) / &
      max(1.0d-09,(xx(i)-xx(i-1)))
      ibgn = i
      exit
    end if
  end do ! i
end do ! iq
if (abs(xi(ni)-xx(nxy))/(xx(nxy)+1.d-15).lt.1.0d-03 ) then
  yi(ni) = yy(nxy)
end if

return
end subroutine simple_interpg


end subroutine rafsck

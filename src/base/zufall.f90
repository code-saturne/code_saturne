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

!     LIBRAIRIE DE SOUS-PROGRAMMES DE GENERATION DE NOMBRES ALEATOIRES
!=======================================================================

!===============================================================================
! Function:
! ---------

!> \file zufall.f90
!>
!> \brief This package downloaded at  http://www.netlib.org/random/zufall.f.
!> It generates random numbers.
!>
!> This package contains a portable random number generator set
!> for: uniform (u in [0,1)), normal (<g> = 0, <g^2> = 1), and
!> Poisson distributions. The basic module, the uniform generator,
!> uses a lagged Fibonacci series generator:
!> \f{eqnarray*}{
!>              t   & = & u(n-273) + u(n-607) \\
!>              u(n)& = & t - float(int(t))
!>
!> \f}
!> where each number generated, u(k), is floating point. Since
!> the numbers are floating point, the left end boundary of the
!> range contains zero. This package is nearly portable except
!> for the following. (1) It is written in lower case, (2) the
!> test package contains a timer (second) which is not portable,
!> and (3) there are cycle times (in seconds) in data statements
!> for NEC SX-3, Fujitsu VP2200, and Cray Y-MP. Select your
!> favorite and comment out the others. Replacement functions
!> for 'second' are included - comment out the others. Otherwise
!> the package is portable and returns the same set of floating
!> point numbers up to word precision on any machine. There are
!> compiler directives ($cdir for Cray, *vdir for SX-3, and VOCL
!> for Fujitsu VP2200) which should be otherwise ignored.
!>
!> To compile this beast, note that all floating point numbers
!> are declared 'double precision'. On Cray X-MP, Y-MP, and C-90
!> machines, use the cft77 (cf77) option -dp to run this in 64
!> bit mode (not 128 bit double).
!>
!> External documentation, "Lagged Fibonacci Random Number Generators
!> for the NEC SX-3," is to be published in the International
!> Journal of High Speed Computing (1994). Otherwise, ask the
!> author:
!>
!>         W. P. Petersen
!>         IPS, RZ F-5
!>         ETHZ
!>         CH 8092, Zurich
!>         Switzerland
!>
!> e-mail:  wpp@ips.ethz.ch.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             seed
!> \param[in]     a             random number
!_______________________________________________________________________________

!===============================================================================

subroutine zufall(n,a)
implicit none

! portable lagged Fibonacci series uniform random number
! generator with "lags" -273 und -607:

!       t    = u(i-273)+buff(i-607)  (floating pt.)
!       u(i) = t - float(int(t))

double precision a(*)
double precision buff(607)
double precision t
integer i,k,ptr,vl,k273,k607
integer buffsz,nn,n,left,q,qq
integer aptr,aptr0,bptr

common /klotz0/buff,ptr
data buffsz/607/

aptr = 0
nn   = n

1     continue

if(nn .le. 0) return

! factor nn = q*607 + r

q    = (nn-1)/607
left = buffsz - ptr

if(q .le. 1) then

! only one or fewer full segments

   if(nn .lt. left) then
      do 2 i=1,nn
         a(i+aptr) = buff(ptr+i)
2           continue
      ptr  = ptr + nn
      return
   else
      do 3 i=1,left
         a(i+aptr) = buff(ptr+i)
3           continue
      ptr  = 0
      aptr = aptr + left
      nn   = nn - left
!  buff -> buff case
      vl   = 273
      k273 = 334
      k607 = 0
      do 4 k=1,3
!dir$ ivdep
         do 5 i=1,vl
            t            = buff(k273+i) + buff(k607+i)
            buff(k607+i) = t - float(int(t))
5              continue
         k607 = k607 + vl
         k273 = k273 + vl
         vl   = 167
         if(k.eq.1) k273 = 0
4           continue

      goto 1
   endif
else

! more than 1 full segment

    do 6 i=1,left
       a(i+aptr) = buff(ptr+i)
6         continue
    nn   = nn - left
    ptr  = 0
    aptr = aptr+left

! buff -> a(aptr0)

    vl   = 273
    k273 = 334
    k607 = 0
    do 7 k=1,3
       if(k.eq.1)then
          do 8 i=1,vl
             t         = buff(k273+i) + buff(k607+i)
             a(aptr+i) = t - float(int(t))
8               continue
          k273 = aptr
          k607 = k607 + vl
          aptr = aptr + vl
          vl   = 167
       else
!dir$ ivdep
          do 9 i=1,vl
             t         = a(k273+i) + buff(k607+i)
             a(aptr+i) = t - float(int(t))
9               continue
          k607 = k607 + vl
          k273 = k273 + vl
          aptr = aptr + vl
       endif
7         continue
    nn = nn - 607

! a(aptr-607) -> a(aptr) for last of the q-1 segments

    aptr0 = aptr - 607
    vl    = 607

    do 10 qq=1,q-2
       k273 = 334 + aptr0
!dir$ ivdep
       do 11 i=1,vl
          t         = a(k273+i) + a(aptr0+i)
          a(aptr+i) = t - float(int(t))
11           continue
       nn    = nn - 607
       aptr  = aptr + vl
       aptr0 = aptr0 + vl
10        continue

! a(aptr0) -> buff, last segment before residual

    vl   = 273
    k273 = 334 + aptr0
    k607 = aptr0
    bptr = 0
    do 12 k=1,3
       if(k.eq.1) then
          do 13 i=1,vl
             t            = a(k273+i) + a(k607+i)
             buff(bptr+i) = t - float(int(t))
13              continue
          k273 = 0
          k607 = k607 + vl
          bptr = bptr + vl
          vl   = 167
       else
!dir$ ivdep
          do 14 i=1,vl
             t            = buff(k273+i) + a(k607+i)
             buff(bptr+i) = t - float(int(t))
14              continue
          k607 = k607 + vl
          k273 = k273 + vl
          bptr = bptr + vl
       endif
12        continue
    goto 1
endif
end subroutine

!===============================================================================

!> \brief Generates initial seed buffer by linear congruential
!>  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
!>  variable seed should be 0 < seed <31328
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     seed          seed
!_______________________________________________________________________________

subroutine zufalli(seed)
implicit none

integer seed
integer ptr
double precision s,t
double precision buff(607)
integer ij,kl,i,ii,j,jj,k,l,m
common /klotz0/buff,ptr
data ij/1802/,kl/9373/

if(seed.ne.0) ij = seed

i = mod(ij/177,177) + 2
j = mod(ij,177) + 2
k = mod(kl/169,178) + 1
l = mod(kl,169)
do 1 ii=1,607
   s = 0.0
   t = 0.5
   do 2 jj=1,24
      m = mod(mod(i*j,179)*k,179)
      i = j
      j = k
      k = m
      l = mod(53*l+1,169)
      if(mod(l*m,64).ge.32) s = s+t
      t = .5*t
2        continue
   buff(ii) = s
1     continue
return
end subroutine

!===============================================================================

!> \brief Box-Muller method for Gaussian random numbers
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             seed
!> \param[out]    x             random number
!_______________________________________________________________________________

subroutine normalen(n,x)
implicit none

double precision x(*)
double precision xbuff(1024)
integer i,ptr,xptr,first
integer buffsz,nn,n,left
common /klotz1/xbuff,first,xptr
data buffsz/1024/

nn   = n
if(nn .le. 0) return
if(first.eq.0)then
   call normal00
   first = 1
endif
ptr = 0

1     continue
left = buffsz - xptr
if(nn .lt. left) then
   do 2 i=1,nn
      x(i+ptr) = xbuff(xptr+i)
2        continue
   xptr = xptr + nn
   return
else
   do 3 i=1,left
      x(i+ptr) = xbuff(xptr+i)
3        continue
   xptr = 0
   ptr  = ptr+left
   nn   = nn - left
   call normal00
   goto 1
endif
end subroutine

!===============================================================================

subroutine normal00
implicit none
double precision pi,twopi
parameter(pi=3.141592653589793)
double precision xbuff(1024),r1,r2,t1,t2
integer first,xptr,i
common /klotz1/xbuff,first,xptr

twopi = 2.*pi
call zufall(1024,xbuff)
do 1 i=1,1024,2
   r1         = twopi*xbuff(i)
   t1         = cos(r1)
   t2         = sin(r1)
   r2         = sqrt(-2.*log(1.-xbuff(i+1)))
   xbuff(i)   = t1*r2
   xbuff(i+1) = t2*r2
1     continue
return
end subroutine

!===============================================================================

!> \brief Poisson generator for distribution function of p's:
!>  \f$ q(mu,p) = exp(-mu) mu**p/p \f$
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n             seed
!> \param[in]     mu            parameter
!> \param[out]    p             random number
!_______________________________________________________________________________

subroutine fische(n,mu,p)
implicit none
integer p(*)
integer indx(1024)
integer n,i,ii,jj,k,left,nl0,nsegs,p0
double precision u(1024),q(1024)
double precision q0,pmu,mu

! initialize arrays, pointers

if (n.le.0) return

pmu = exp(-mu)
p0  = 0

nsegs = (n-1)/1024
left  = n - nsegs*1024
nsegs = nsegs + 1
nl0   = left

do 2 k = 1,nsegs

   do 3 i=1,left
      indx(i)    = i
      p(p0+i)    = 0
      q(i)       = 1.0
3        continue

! Begin iterative loop on segment of p's

1        continue

! Get the needed uniforms

   call zufall(left,u)

   jj = 0

!dir$ ivdep
   do 4 i=1,left
      ii    = indx(i)
      q0    = q(ii)*u(i)
      q(ii) = q0
      if( q0.gt.pmu ) then
         jj       = jj + 1
         indx(jj) = ii
         p(p0+ii) = p(p0+ii) + 1
      endif
4        continue

! any left in this segment?

   left = jj
   if(left.gt.0)then
      goto 1
   endif

   p0    = p0 + nl0
   nl0   = 1024
   left  = 1024

2     continue

return
end subroutine

!===============================================================================

block data
implicit none

! globally accessable, compile-time initialized data

integer ptr,xptr,first
double precision buff(607),xbuff(1024)
common /klotz0/buff,ptr
common /klotz1/xbuff,first,xptr
data ptr/0/,xptr/0/,first/0/
end block data

!===============================================================================

!----------------------------
! End of the library
!----------------------------


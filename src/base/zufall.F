c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
C     LIBRAIRIE DE SOUS-PROGRAMMES DE GENERATION DE NOMBRES ALEATOIRES
C     ****************************************************************
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC
CFONC Ce package de sous-programmes a ete telechargee sur
CFONC http://www.netlib.org/random/zufall.f
CFONC
CFONC
CFONC README for zufall random number package
CFONC ------ --- ------ ------ ------ -------
CFONC This package contains a portable random number generator set
CFONC for: uniform (u in [0,1)), normal (<g> = 0, <g^2> = 1), and
CFONC Poisson distributions. The basic module, the uniform generator,
CFONC uses a lagged Fibonacci series generator:
CFONC
CFONC               t    = u(n-273) + u(n-607)
CFONC               u(n) = t - float(int(t))
CFONC
CFONC where each number generated, u(k), is floating point. Since
CFONC the numbers are floating point, the left end boundary of the
CFONC range contains zero. This package is nearly portable except
CFONC for the following. (1) It is written in lower case, (2) the
CFONC test package contains a timer (second) which is not portable,
CFONC and (3) there are cycle times (in seconds) in data statements
CFONC for NEC SX-3, Fujitsu VP2200, and Cray Y-MP. Select your
CFONC favorite and comment out the others. Replacement functions
CFONC for 'second' are included - comment out the others. Otherwise
CFONC the package is portable and returns the same set of floating
CFONC point numbers up to word precision on any machine. There are
CFONC compiler directives ($cdir for Cray, *vdir for SX-3, and VOCL
CFONC for Fujitsu VP2200) which should be otherwise ignored.
CFONC
CFONC To compile this beast, note that all floating point numbers
CFONC are declared 'double precision'. On Cray X-MP, Y-MP, and C-90
CFONC machines, use the cft77 (cf77) option -dp to run this in 64
CFONC bit mode (not 128 bit double).
CFONC
CFONC External documentation, "Lagged Fibonacci Random Number Generators
CFONC for the NEC SX-3," is to be published in the International
CFONC Journal of High Speed Computing (1994). Otherwise, ask the
CFONC author:
CFONC
CFONC          W. P. Petersen
CFONC          IPS, RZ F-5
CFONC          ETHZ
CFONC          CH 8092, Zurich
CFONC          Switzerland
CFONC
CFONC e-mail:  wpp@ips.ethz.ch.
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C
C***********************************************************************
C
      subroutine zufall(n,a)
      implicit none
c
c portable lagged Fibonacci series uniform random number
c generator with "lags" -273 und -607:
c
c       t    = u(i-273)+buff(i-607)  (floating pt.)
c       u(i) = t - float(int(t))
c
      double precision a(*)
      double precision buff(607)
      double precision t
      integer i,k,ptr,VL,k273,k607
      integer buffsz,nn,n,left,q,qq
      integer aptr,aptr0,bptr
c
      common /klotz0/buff,ptr
      data buffsz/607/
c
      aptr = 0
      nn   = n
c
1     continue
c
      if(nn .le. 0) return
c
c factor nn = q*607 + r
c
      q    = (nn-1)/607
      left = buffsz - ptr
c
      if(q .le. 1) then
c
c only one or fewer full segments
c
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
c  buff -> buff case
            VL   = 273
            k273 = 334
            k607 = 0
            do 4 k=1,3
cdir$ ivdep
               do 5 i=1,VL
                  t            = buff(k273+i) + buff(k607+i)
                  buff(k607+i) = t - float(int(t))
5              continue
               k607 = k607 + VL
               k273 = k273 + VL
               VL   = 167
               if(k.eq.1) k273 = 0
4           continue
c
            goto 1
         endif
      else
c
c more than 1 full segment
c
          do 6 i=1,left
             a(i+aptr) = buff(ptr+i)
6         continue
          nn   = nn - left
          ptr  = 0
          aptr = aptr+left
c
c buff -> a(aptr0)
c
          VL   = 273
          k273 = 334
          k607 = 0
          do 7 k=1,3
             if(k.eq.1)then
                do 8 i=1,VL
                   t         = buff(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
8               continue
                k273 = aptr
                k607 = k607 + VL
                aptr = aptr + VL
                VL   = 167
             else
cdir$ ivdep
                do 9 i=1,VL
                   t         = a(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
9               continue
                k607 = k607 + VL
                k273 = k273 + VL
                aptr = aptr + VL
             endif
7         continue
          nn = nn - 607
c
c a(aptr-607) -> a(aptr) for last of the q-1 segments
c
          aptr0 = aptr - 607
          VL    = 607
c
          do 10 qq=1,q-2
             k273 = 334 + aptr0
cdir$ ivdep
             do 11 i=1,VL
                t         = a(k273+i) + a(aptr0+i)
                a(aptr+i) = t - float(int(t))
11           continue
             nn    = nn - 607
             aptr  = aptr + VL
             aptr0 = aptr0 + VL
10        continue
c
c a(aptr0) -> buff, last segment before residual
c
          VL   = 273
          k273 = 334 + aptr0
          k607 = aptr0
          bptr = 0
          do 12 k=1,3
             if(k.eq.1) then
                do 13 i=1,VL
                   t            = a(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
13              continue
                k273 = 0
                k607 = k607 + VL
                bptr = bptr + VL
                VL   = 167
             else
cdir$ ivdep
                do 14 i=1,VL
                   t            = buff(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
14              continue
                k607 = k607 + VL
                k273 = k273 + VL
                bptr = bptr + VL
             endif
12        continue
          goto 1
      endif
      end
c
C=======================================================================
c
      subroutine zufalli(seed)
      implicit none
c
c  generates initial seed buffer by linear congruential
c  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
c  variable seed should be 0 < seed <31328
c
      integer seed
      integer ptr
      double precision s,t
      double precision buff(607)
      integer ij,kl,i,ii,j,jj,k,l,m
      common /klotz0/buff,ptr
      data ij/1802/,kl/9373/
c
      if(seed.ne.0) ij = seed
c
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
      end
c
C=======================================================================
c
      subroutine normalen(n,x)
      implicit none
c
c Box-Muller method for Gaussian random numbers
c
      double precision x(*)
      double precision xbuff(1024)
      integer i,ptr,xptr,first
      integer buffsz,nn,n,left
      common /klotz1/xbuff,first,xptr
      data buffsz/1024/
c
      nn   = n
      if(nn .le. 0) return
      if(first.eq.0)then
         call normal00
         first = 1
      endif
      ptr = 0
c
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
      end
c
C=======================================================================
c
      subroutine normal00
      implicit none
      double precision pi,twopi
      parameter(pi=3.141592653589793)
      double precision xbuff(1024),r1,r2,t1,t2
      integer first,xptr,i
      common /klotz1/xbuff,first,xptr
c
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
      end
c
C=======================================================================
c
      subroutine fische(n,mu,p)
      implicit none
      integer p(*)
      integer indx(1024)
      integer n,i,ii,jj,k,left,nl0,nsegs,p0
      double precision u(1024),q(1024)
      double precision q0,pmu,mu
c
c Poisson generator for distribution function of p's:
c
c    q(mu,p) = exp(-mu) mu**p/p!
c
c initialize arrays, pointers
c
      if (n.le.0) return
c
      pmu = exp(-mu)
      p0  = 0
c
      nsegs = (n-1)/1024
      left  = n - nsegs*1024
      nsegs = nsegs + 1
      nl0   = left
c
      do 2 k = 1,nsegs
c
         do 3 i=1,left
            indx(i)    = i
            p(p0+i)    = 0
            q(i)       = 1.0
3        continue
c
c Begin iterative loop on segment of p's
c
1        continue
c
c Get the needed uniforms
c
         call zufall(left,u)
c
         jj = 0
c
cdir$ ivdep
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
c
c any left in this segment?
c
         left = jj
         if(left.gt.0)then
            goto 1
         endif
c
         p0    = p0 + nl0
         nl0   = 1024
         left  = 1024
c
2     continue
c
      return
      end
c
C=======================================================================
c
      block data
      implicit none
c
c globally accessable, compile-time initialized data
c
      integer ptr,xptr,first
      double precision buff(607),xbuff(1024)
      common /klotz0/buff,ptr
      common /klotz1/xbuff,first,xptr
      data ptr/0/,xptr/0/,first/0/
      end
c
C=======================================================================
C
C--------------------
C FIN DE LA LIBRAIRIE
C--------------------
C
c@z

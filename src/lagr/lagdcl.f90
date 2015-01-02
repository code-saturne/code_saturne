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

subroutine lagdcl &
!================

 ( dx    , vvue   , vpart  , marko  ,                                &
   tempf , depint ,                                                  &
   dtl   , tstruc , tdiffu , ttotal , vstruc ,                       &
   romp  , taup   , kdif   , tlag2  , yplus  ,                       &
   lvisq ,                                                           &
   unif1 , unif2  , dintrf , rpart  ,                                &
   kdifcl, indint , gnorm  , vnorm  , grpn   , piiln)

!===============================================================================
!
! Purpose:
! --------
!
!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------
!
!
!   Deposition submodel (Guingo & Minier, 2008) :
!
!   Management of the diffusion in the internal zone (y^+ < dintrf)
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! dx               ! r  ! --> ! wall-normal displacement                       !
! vvue             ! r  ! <-> ! wall-normal velocity of the flow seen          !
! vpart            ! r  ! <-> ! particle wall-normal velocity                  !
! marko            ! r  ! <-> ! state of the jump process                      !
! tempf            ! r  ! <-- ! temperature of the fluid                       !
! dtl              ! r  ! --> ! Lagrangian time step                           !
! tstruc           ! r  ! <-> ! coherent structure mean duration               !
! tdiffu           ! r  ! <-> ! diffusion phase mean duration                  !
! ttotal           ! r  ! <-> ! tdiffu + tstruc                                !
! vstruc           ! r  ! <-- ! coherent structure velocity                    !
! romp             ! r  ! --> ! particle density                               !
! taup             ! r  ! <-> ! particle relaxation time                       !
! kdif             ! r  ! <-> ! diffusion phase diffusion coefficient          !
! tlag2            ! r  ! <-> ! diffusion relaxation timescale                 !
! lvisq            ! r  ! <-- ! wall-unit lengthscale                          !
! yplus            ! r  ! <-- ! wall-normal velocity of the flow seen          !
! unif1            ! r  ! <-> ! random number (uniform law)                    !
! unif2            ! r  ! <-> ! random number (uniform law)                    !
! dintrf           ! r  ! <-- ! extern-intern interface location               !
! rpart            ! r  ! <-- ! particle radius                                !
! kdifcl           ! r  ! <-> ! internal zone diffusion coefficient            !
! indint           ! r  ! <-> ! interface indicator                            !
! gnorm            ! r  ! <-- ! wall-normal gravity component                  !
! vnorm            ! r  ! <-- ! wall-normal fluid (Eulerian) velocity          !
! grpn             ! r  ! <-- ! wall-normal pressure gradient                  !
! piiln            ! r  ! <-- ! SDE integration auxiliary term                 !
!__________________!____!_____!________________________________________________!
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array


!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use cstphy
use lagran
use ppthch

!===============================================================================

implicit none

! Arguments

integer          marko , indint
double precision tempf, depint
double precision dx   , vvue  , vpart
double precision dtl,tstruc, tdiffu, ttotal, vstruc
double precision romp , taup  , kdif  , tlag2 , yplus
double precision unif1, unif2 , dintrf, rpart, grpn, piiln
double precision lvisq,  kdifcl, gnorm, vnorm

! Local variables

double precision vvue0 , vpart0
double precision argt , kaux, tci
double precision mpart,kdifbr,kdifbrtp
double precision tlmtp, tlptp, tltp, tl2, tp2
double precision thet, the2, etl, etp
double precision l1l, l1p, l2l, l2p, l3
double precision kaux2
double precision aa1, bb1, cc1, dd1, ee1
double precision tdciu,dtstl,dtstp,k2the2
double precision ubr,xiubr, ucarbr, xcarbr
double precision pgam2, ggam2, ome2, pgagga, pgaome, ggaome
double precision p11, p21, p22, p31, p32, p33
double precision p21br, p11br, p22br
double precision terx, terp, terf, terpbr, terxbr
double precision dtp1 , dxaux , yplusa
double precision argtn1,kauxn1,tcin1,pox1, pox2
double precision aa2, bb2, c2c, a2c, b2c, a22, b22
double precision ketoi,ketoi2
double precision vagaus(3), vagausbr(2), force

!===============================================================================

call normalen(3,vagaus)
call normalen(2,vagausbr)

force = gnorm * taup
tdciu = 0.d0

vvue0  = vvue
vpart0 = vpart
!
if (yplus.lt.5.d0) then

  argt   = pi * yplus / 5.d0
  kaux = kdifcl * 0.5d0 * ( 1.d0 - cos(argt))
  tci = - tlag2**2 * 0.5d0                                        &
        * kdifcl**2 * pi * sin(argt) * (1.d0-cos(argt))           &
        / (2.d0*5.d0)/lvisq
else

  kaux = kdifcl
! Interpolation of the decreasing normal fluid velocity around zero:
  tci = vnorm * yplus / dintrf

endif


!------------------------------------------------
!  Brownian motion
!------------------------------------------------

mpart    = 4.d0/3.d0*pi*rpart**3*romp
kdifbr   = sqrt(2.d0*kboltz*tempf / (mpart * taup))
kdifbrtp = kdifbr * taup

!------------------------------------------------


dtstl = dtl / tlag2
dtstp = dtl / taup

tlmtp = tlag2 - taup
tlptp = tlag2 + taup
tltp  = tlag2*taup
tl2   = tlag2**2
tp2   = taup**2

thet = tlag2 / tlmtp
the2 = thet**2

etl = exp(-dtstl)
etp = exp(-dtstp)

l1l = 1.d0 - etl
l1p = 1.d0 - etp
l2l = 1.d0 - etl * etl
l2p = 1.d0 - etp * etp
l3  = 1.d0 - etl * etp

kaux2  = kaux**2
k2the2 = kaux2*the2

aa1 = taup * l1p
bb1 = thet*(tlag2*l1l - aa1)
cc1 = dtl - aa1 - bb1
dd1 = thet*(etl - etp)
ee1 = l1p

!---------------------------------------------------
! Auxiliary terms for Brownian motion
!---------------------------------------------------

xiubr  = 0.5d0 * (kdifbrtp * l1p) ** 2
ucarbr = kdifbrtp * kdifbr * 0.5d0 * l2p
xcarbr = kdifbrtp**2 * (dtl - l1p*(2.d0 + l1p)*0.5d0*taup)
ubr    = sqrt(max(ucarbr, zero))

!----------------------------------------------------
! Deterministic terms computation
!----------------------------------------------------

vvue  = vvue0 * etl  + tci * l1l
vpart = vpart0*etp+dd1*vvue0+tci*(ee1 - dd1) + force*ee1
dx    = aa1*vpart0 + bb1*vvue0 + cc1*tci + (dtl-aa1)*force

!----------------------------------------------------
! Correlation matrix
!----------------------------------------------------

pgam2 = 0.5d0*kaux2*tlag2*l2l

ggam2 = the2*pgam2                                                &
      + k2the2*(l3*(-2*tltp / tlptp) + l2p*(taup*0.5d0))

ome2  = k2the2*( dtl*tlmtp**2                                     &
        + l2l*(tl2*tlag2*0.5d0)                                   &
        + l2p*(tp2*taup *0.5d0)                                   &
        + l1l*(-2.*tl2*tlmtp)                                     &
        + l1p*(2.*tp2*tlmtp)                                      &
        + l3*(-2.*(tltp**2)/tlptp) )

pgagga = thet*(pgam2 - kaux2*tltp/tlptp*l3)

pgaome = thet*tlag2*(-pgam2 + kaux2*(l1l*tlmtp+l3*tp2/tlptp))

ggaome = k2the2*( tlmtp*(tlag2*l1l + l1p*(-taup))                 &
         + l2l*(-tl2*0.5d0)                                       &
         + l2p*(-tp2*0.5d0)                                       &
         + l3*tltp )


!------------------------------------------------------
!  Choleski decomposition
!------------------------------------------------------

!  P11 Computation
!
p11 = sqrt(max(zero,pgam2))

!  P21 and P22 computations

if (abs(p11).gt.epzero) then
  p21 = pgagga / p11
else
  p21 = 0.d0
endif

p22 = sqrt(max(zero,ggam2 - p21**2))

!  P31, P32 and P33 computations

if (abs(p11).gt.epzero) then
  p31 = pgaome / p11
else
  p31 = 0.d0
endif

if (abs(p22).gt.epzero) then
  p32 = (ggaome - p21 * p31) / p22
else
  p32 = 0.d0
endif

p33 = sqrt(max(zero,ome2 - p31**2 - p32**2))

!-----------------------------------------------------------
!  Brownian motion term
!-----------------------------------------------------------
p11br = ubr

if (abs(p11br).gt.epzero) then
  p21br = xiubr / p11br
else
  p21br = 0.d0
endif

p22br = sqrt(max(xcarbr - p21br **2, zero))


!-----------------------------------------------------------
!  The random terms are consequently:
!-----------------------------------------------------------

terf = p11*vagaus(1)
terp = p21*vagaus(1)+ p22*vagaus(2)
terx = p31*vagaus(1)+ p32*vagaus(2)+p33*vagaus(3)

!-----------------------------------------------------------
!  Brownian motion
!-----------------------------------------------------------

terpbr = p11br * vagausbr(1)
terxbr = p21br * vagausbr(1) + p22br * vagausbr(2)

!------------------------------------------------------------
!  Finalization (First order)
!-----------------------------------------------------------

vvue   = vvue  + terf
vpart  = vpart + terp + terpbr
dx     = dx    + terx + terxbr

yplusa = yplus - dx / lvisq

if ((yplusa*lvisq).lt.rpart) then

  dx = dx + 2*rpart
  goto 65

endif

if ((yplusa.gt.dintrf).and.(indint.ne.1)) then

  marko = 2
  vvue  = -sqrt((kdifcl*(ttotal/tdiffu))**2*tlag2/2.d0)           &
              *sqrt(2.d0*pi)*0.5d0

  dx    = dx*(dintrf - yplus)/(yplusa - yplus)
  vpart =  (yplus - yplusa)*lvisq / dtl

  dxaux = dx
  dtp1  = dtl*(dintrf - yplusa)/(yplus - yplusa)
  yplus = dintrf

  call lagdif                                                     &
  !==========
        ( dx     , vvue   , vpart  , marko  ,                     &
          tempf  , depint ,                                       &
          dtp1   , tstruc , tdiffu , ttotal , vstruc ,            &
          romp   , taup   , kdif   , tlag2  , lvisq  , yplus ,    &
          unif1  , unif2  , dintrf , rpart  ,                     &
          kdifcl , indint , gnorm  , vnorm  , grpn   , piiln)

  dx = dxaux + dx

else

  if (yplusa.gt.0.d0) then

    if (yplusa.lt.5.) then

      argtn1 = pi * yplusa / 5.d0
      kauxn1 = kdifcl * 0.5d0 * ( 1.d0 - cos(argtn1))
      tcin1  = - tlag2**2 * 0.5d0                                 &
             * kdifcl**2 * pi * sin(argtn1) * (1.d0 - cos(argtn1))&
             / (2.d0*5.d0)/lvisq

    else

      kauxn1 = kdifcl
      tcin1 = 0.d0

    endif

!--------------------------------------------------------------------------
!   Auxiliary computation
!--------------------------------------------------------------------------

    pox1 = l1l / dtstl
    pox2 = tlptp / dtl * l1p

    aa2 = - etl + pox1

    bb2 = 1.d0 - pox1

    c2c = tlag2 / tlmtp * (etl - etp)

    a2c = - etp + pox2 - (1.d0 + tlag2 / dtl)*c2c

    b2c = 1.d0 - pox2 + (tlag2 / dtl)*c2c

    a22 = l2l + l2l/(2*dtstl) - 1.d0

    b22 = 1.d0 - l2l /(2*dtstl)

!--------------------------------------------------------------------------
!   Deterministic terms computation
!--------------------------------------------------------------------------

    vvue  = vvue0 * etl + aa2 * tci + bb2 * tcin1
    vpart =  vpart0*etp + vvue0*c2c                           &
         + a2c*tci    + b2c*tcin1                             &
         + force * (1.d0 - (etp - 1.d0) / (-dtstp))


!---------------------------------------------------
!  Diffusion coefficient computation
!---------------------------------------------------

    ketoi = (a22 * kaux + b22 * kauxn1) / l2l
    ketoi2 = ketoi**2

!---------------------------------------------------
!  Correlation matrix computation
!---------------------------------------------------

    pgam2 = 0.5d0*ketoi2*tlag2*l2l

    ggam2 = the2*(pgam2 + ketoi2*(l3*(-2.*tltp / tlptp)           &
           +l2p*taup*0.5d0))

    pgagga = thet*(pgam2 - ketoi2*tltp/tlptp*l3)

!------------------------------------------------------
!  Choleski decomposition
!------------------------------------------------------

!  P11 computation
!
    p11 = sqrt(max(zero,pgam2))

!  P21 and P22 computation

    if (abs(p11).gt.epzero) then
      p21 = pgagga / p11
    else
      p21 = 0.d0
    endif

    p22 = sqrt(max(zero,ggam2 - p21**2))

!-----------------------------------------------------------
!  The random terms are:
!-----------------------------------------------------------

    terf = p11 * vagaus(1)
    terp = p21 * vagaus(1)  + p22 * vagaus(2)

!------------------------------------------------------------
!  Finalization (Second order)
!-----------------------------------------------------------
    vvue  = vvue  + terf
    vpart = vpart + terp + terpbr

  endif
!
endif

65 continue

end subroutine

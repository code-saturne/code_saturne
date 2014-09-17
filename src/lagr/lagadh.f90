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

subroutine lagadh &
!================

 ( ip     ,                                                       &
   rtp    , adhesion_energ)

!===============================================================================

! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------


!   Calculation of the adhesion force and adhesion energy
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ip               ! i  ! <-- ! particle number                                !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! adhesion_energ   ! r  ! --> ! particle adhesion energy                       !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================


!===============================================================================


!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: nvar
use cstphy
use cstnum
use lagpar
use lagran
use ppthch
use entsor
use mesh
use optcal
use numvar

!===============================================================================

implicit none

! Arguments

integer          ip

double precision rtp (ncelet,nflown:nvar)
double precision adhesion_energ

! Local variables

integer nbasg, nbasp, np, ntmp(1)

double precision step, rpart, rtmp(1)
double precision paramh, nmoyap, nmoyag, scovag, scovap
double precision dismin, distcc, distp
double precision udlvor(2), uvdwsp, uvdwss, uedlsp, uedlss
double precision fadhes


! Variables for the adhesion moment
double precision dismom, omsurf

integer iel, mode
double precision tempf
! ==========================================================================
! 0.    initialization
! ==========================================================================

!     step = step used to calculate the adhesion force following
!                         F = U(dcutof+step)-U(dcutof-step)/(2*step)

step = 1.0d-11

scovap = denasp * pi * rayasp**2
scovag = pi * rayasg**2 / espasg**2

! Determination of the temperature


  iel = ipepa(jisor,ip)

  if (iscalt.gt.0) then
    if (itherm.eq.1) then
      if (itpscl.eq.2) then
        tempf = rtp(iel,isca(iscalt)) + tkelvi
      else if (itpscl.eq. 1) then
        tempf = rtp(iel,isca(iscalt))
      endif
    else if (itherm.eq.2) then
      mode = 1
      call usthht(mode,rtp(iel,isca(iscalt)),tempf)
    endif
  else
    tempf = t0
  endif


! ==========================================================================
! 3.    calculation of the adhesion force
! ==========================================================================


!     determination of the number of contacts with asperities
!     =======================================================

! Number of large-scale asperities

rpart = 0.5d0 * eptp(jdp,ip)

nmoyag = (2.0d0 * rpart + rayasg) / rayasg * scovag

if (nmoyag.gt.600.d0) then
   call normalen(1,ntmp)
   ipepa(jnbasg,ip) = nmoyag + sqrt(nmoyag)*ntmp(1)
   ipepa(jnbasg,ip) = max(0,ipepa(jnbasg,ip))
else
   call fische1(nmoyag, ipepa(jnbasg,ip))
endif


if (ipepa(jnbasg,ip).gt.1) then

   nmoyag = 1 + 2.0d0 * dcutof*(2.0d0*rpart + 2.0d0 * rayasg+4.0d0*dcutof)       &
        / rayasg**2 * scovag

   if (nmoyag.gt.600.d0) then
     call normalen(1,ntmp)
     nbasg = nmoyag + sqrt(nmoyag)*ntmp(1)
     nbasg = max(0,nbasg)
   else
     call fische(1, nmoyag, ntmp)
     nbasg = ntmp(1)
   endif

   nbasg = max(1,nbasg)

else
   nbasg = ipepa(jnbasg,ip)
endif

! Nb of small-scale asperities : 1st case: no large-scale asperities

if (nbasg.eq.0) then

   nmoyap = (2.0d0 * rpart + rayasp) / rayasp * scovap

   if (nmoyap.gt.600.d0) then
      call normalen(1,ntmp)
      ipepa(jnbasp,ip) = nmoyap + sqrt(nmoyap)*ntmp(1)
      ipepa(jnbasp,ip) = max(0,ipepa(jnbasp,ip))
   else
      call fische1(nmoyap, ipepa(jnbasp,ip))
   endif

   if (ipepa(jnbasp,ip).gt.1) then

      nmoyap = 1 + 2.0d0*dcutof*(2.0d0*rpart+2.0d0*rayasp+4.0d0*dcutof)    &
           / rayasp**2 * scovap

      if (nmoyap.gt.600.d0) then
         call normalen(1,ntmp)
         ipepa(jnbasp,ip) = nmoyap + sqrt(nmoyap)*ntmp(1)
         ipepa(jnbasp,ip) = max(0,ipepa(jnbasp,ip))
      else
         call fische(1, nmoyap, ntmp)
         nbasp = ntmp(1)
      endif
      nbasp = max(1,nbasp)

   else
      nbasp = ipepa(jnbasp,ip)
   endif

   ! Determination of the minimal distance between the particle and the plate
   dismin = rayasp * min(1.0d0,ipepa(jnbasp,ip)*1.0d0)

   ! 2nd case: contact with large-scale asperities

else

   paramh = 0.5d0*(2.0d0*rpart+rayasp)*rayasp / (rpart + rayasg)

   nmoyap = paramh*(2*rayasg-paramh) / rayasp**2 * scovap

   if (nmoyap.gt.600.d0) then
      call normalen(1,ntmp)
      ipepa(jnbasp,ip) = nmoyap + sqrt(nmoyap)*ntmp(1)
      ipepa(jnbasp,ip) = max(0,ipepa(jnbasp,ip))
   else
      call fische1(nmoyap, ipepa(jnbasp,ip))
   endif


   if (ipepa(jnbasp,ip).gt.1) then

      paramh = 0.5d0*(2.0d0*rpart+2*rayasp+4.0d0*dcutof)*2.0d0*dcutof     &
           / (rpart+rayasg+rayasp+dcutof)

      nmoyap = 1 + paramh*(2*rayasg-paramh) / rayasp**2 * scovap

      if (nmoyap.gt.600.d0) then
         call normalen(1,ntmp)
         nbasp = nmoyap + sqrt(nmoyap)*ntmp(1)
         nbasp = max(0,nbasp)
      else
         call fische(1, nmoyap, ntmp)
         nbasp = ntmp(1)
      endif
      nbasp = max(1,nbasp)
   else
      nbasp = ipepa(jnbasp,ip)
   endif

   ! Mutliple contacts with large scale asperities?

   nbasp = nbasp * nbasg
   ipepa(jnbasp,ip) = ipepa(jnbasp,ip)*ipepa(jnbasg,ip)

   ! Determination of the minimal distance between the particle and the plate
   dismin = rayasp * min(1.0d0,nbasp*1.0d0)                  &
        + rayasg * min(1.0d0,nbasg*1.0d0)

endif ! End of determination of ipepa(jnbasp,ip) and ipepa(jnbasg,ip)


! Sum of {particle-plane} and {particle-asperity} energies


! Interaction between the sphere and the plate
do np = 1,2
   udlvor(np) = 0.0d0
   distp = dismin + dcutof + step * (3-2*np)

   call vdwsp(distp, rpart, uvdwsp)
   call edlsp(distp, rpart, tempf, uedlsp)

   udlvor(np) = (uvdwsp + uedlsp) * (1 - scovag - scovap)
enddo

fadhes = (udlvor(2) - udlvor(1)) / (2.d0 * step)
adhesion_energ = udlvor(1)

! Interaction between the sphere and small-scale asperities

do np = 1,2

   udlvor(np) = 0.0d0
   distcc =  dcutof + step * (3-2*np) + rpart + rayasp

   call vdwsa(distcc, rpart, rayasp, uvdwss)
   call edlsa(distcc, rpart, rayasp, tempf , uedlss)

   udlvor(np) = uvdwss + uedlss

enddo

fadhes = fadhes + (udlvor(2) - udlvor(1)) / (2.d0 * step)  * nbasp

adhesion_energ = adhesion_energ + udlvor(1)*nbasp

! Interaction between the sphere and large-scale asperities

do np = 1,2
   udlvor(np) = 0.0d0

   if (nbasp.eq.0) then
      distcc =  dcutof + step * (3-2*np) + rpart + rayasg
   elseif (nbasp.gt.0) then
      distcc =  dcutof + rayasp + step * (3-2*np) + rpart + rayasg
   endif

   call vdwsa(distcc, rpart, rayasg, uvdwss)
   call edlsa(distcc, rpart, rayasg, tempf , uedlss)

   udlvor(np) = uvdwss + uedlss
enddo

fadhes = fadhes + (udlvor(2) - udlvor(1)) / (2.0d0 * step) * nbasg
adhesion_energ = adhesion_energ + udlvor(1) * nbasg

! The force is negative when it is attractive

if (fadhes.ge.0.0d0) then
   pepa(jfadh,ip) = 0.0d0
else
   pepa(jfadh,ip) = - fadhes
endif

! The interaction should be negative to prevent reentrainment (attraction)

if (adhesion_energ.ge.0.0d0) then
   adhesion_energ = 0.0d0
else
   adhesion_energ = abs(adhesion_energ)
endif

!
! Calculation of adhesion torques exerted on the particle

call zufall(1,rtmp)
dismom = rtmp(1)

if (nbasg.gt.0) then
   dismom = dismom * sqrt((2.0d0*rpart+rayasg)*rayasg)
elseif (nbasg.eq.0 .and. nbasp.gt.0) then
   dismom = dismom * sqrt((2.0d0*rpart+rayasp)*rayasp)
else

   !in the sphere-plate case, we use the deformation given by the DMT theory,
   !which is close to our approach

   omsurf = cstham / (24.0d0 * pi * dcutof**2)
   dismom = (4.0d0 * pi * omsurf * (rpart**2)/modyeq)**(1.0d0/3.0d0)

endif

pepa(jmfadh,ip) = pepa(jfadh,ip)*dismom


end subroutine lagadh



! =========================================================================
!     vdw interaction between a sphere and a plane
!     following formulas from Czarnecki (large distances)
!                           and Gregory (small distances)
! =========================================================================

subroutine vdwsp (distp, rpart, var)

use cstnum
use lagran

implicit none

double precision distp, rpart, var


if (distp.lt.lambwl/2/pi) then
   var = -cstham*rpart/(6*distp)*(1/                                &
        (1+14*distp/lambwl+5*pi/4.9d0*distp**3/lambwl/rpart**2))
else
   var = cstham*(2.45/60/pi*lambwl*((distp-rpart)/distp**2          &
            -(distp+3*rpart)/(distp+2*rpart)**2)                  &
            -2.17/720/pi**2*lambwl**2*((distp-2*rpart)            &
            /distp**3 -(distp+4*rpart)/(distp+2*rpart)**3)        &
            +0.59/5040/pi**3*lambwl**3*((distp-3*rpart)/          &
            distp**4 -(distp+5*rpart)/(distp+2*rpart)**4))
endif

end subroutine vdwsp


! =========================================================================
!     Vdw interaction between two spheres
!     following the formula from Gregory (1981a)
! =========================================================================

subroutine vdwsa (distcc, rpart1, rpart2, var)

use cstnum
use lagran

implicit none

double precision distcc, rpart1,rpart2, var

var = - cstham*rpart1*rpart2/(6*(distcc-rpart1-rpart2)              &
           *(rpart1+rpart2))*(1-5.32d0*(distcc-rpart1-rpart2)     &
           /lambwl*log(1+lambwl/(distcc-rpart1-rpart2)/5.32d0))

end subroutine vdwsa



! =========================================================================
!     EDL INTERACTION BETWEEN A SPHERE AND A PLANE
!     following the formula from Bell & al (1970)
!     based on the McCartney & Levine method
! =========================================================================

subroutine edlsp &
!     -----------------
     (distp, rpart, tempf , var)

use cstnum
use lagran
use ppthch
use entsor

implicit none

double precision alpha, omega1, omega2, gamma
double precision charge
double precision distp, rpart, var , tau
double precision ldebye, tempf ,  lphi1, lphi2

charge = 1.6d-19

ldebye = ((2.d3 * cstfar**2 * fion)                            &
           /(epseau * epsvid * rr * tempf))**(-0.5)


! Reduced zeta potential
lphi1 =  valen * charge * phi1 /  kboltz / tempf
lphi2 =  valen * charge * phi2 /  kboltz / tempf


!Extended reduced zeta potential
!  (following the work from Ohshima et al, 1982, JCIS, 90, 17-26)

!For the particle
tau = rpart / ldebye

lphi1 = 8.d0 * tanh(lphi1 / 4.d0) /                       &
( 1.d0 + sqrt(1.d0 - (2.d0 * tau + 1.d0) / (tau + 1)**2   &
* tanh(lphi1 / 4.d0)**2) )

! For the plate
lphi2 = 4.d0 * tanh(lphi2 / 4.d0)

!Calculation for the EDL force
alpha = sqrt((distp+rpart)/rpart)+sqrt(rpart/(distp+rpart))
omega1 = lphi1**2 + lphi2**2 + alpha * lphi1 * lphi2
omega2 = lphi1**2 + lphi2**2 - alpha * lphi1 * lphi2
gamma = sqrt(rpart/(distp+rpart))*exp(-1.d0/ldebye*distp)

var = 2 * pi * epseau * epsvid * (kboltz * tempf / charge)**2     &
     * rpart * (distp + rpart) / (distp + 2 * rpart)              &
     * (omega1 * log ( 1 + gamma) + omega2 * log( 1 - gamma))

 end subroutine edlsp


! =========================================================================
!     EDL INTERACTION BETWEEN TWO SPHERES
!     following the formula from Bell & al (1970)
!     based on the McCartney & Levine method
! =========================================================================

subroutine edlsa &
!     -----------------
     (distcc, rpart1, rpart2, tempf, var)

use cstnum
use lagran
use ppthch

implicit none

double precision alpha, omega1, omega2, gamma
double precision charge, tau
double precision distcc, rpart1, rpart2, var
double precision ldebye, tempf , lphi1, lphi2

charge = 1.6e-19

ldebye = ((2.d3 * cstfar**2 * fion)                            &
           /(epseau * epsvid * rr * tempf))**(-0.5)


! Reduced zeta potential
lphi1 =  valen * charge * phi1 /  kboltz / tempf
lphi2 =  valen * charge * phi2 /  kboltz / tempf

!Extended reduced zeta potential
!  (following the work from Ohshima et al, 1982, JCIS, 90, 17-26)

! For the first particle
tau = rpart1 / (1.d0/ldebye)

lphi1 = 8.d0 * tanh(lphi1 / 4.d0) /                       &
( 1.d0 + sqrt(1.d0 - (2.d0 * tau + 1.d0) / (tau + 1)**2   &
* tanh(lphi1 / 4.d0)**2))

! For the second particle
tau = rpart2 / ldebye

lphi2 = 8.d0 * tanh(lphi2 / 4.d0) /                       &
( 1.d0 + sqrt(1.d0 - (2.d0 * tau + 1.d0) / (tau + 1)**2   &
* tanh(lphi2 / 4.d0)**2) )

! Calculation of the EDL force
alpha = sqrt(rpart2*(distcc-rpart2)/(rpart1*(distcc-rpart1)))     &
         +sqrt(rpart1*(distcc-rpart1)/(rpart2*(distcc-rpart2)))
omega1 = lphi1**2 + lphi2**2 + alpha * lphi1 * lphi2
omega2 = lphi1**2 + lphi2**2 - alpha * lphi1 * lphi2
gamma = sqrt(rpart1*rpart2/(distcc-rpart1)/(distcc-rpart2))       &
        *exp(1.d0/ldebye*(rpart1+rpart2-distcc))

var = 2 * pi * epseau * epsvid * (kboltz * tempf / charge)**2     &
     * rpart1 * rpart2 * (distcc - rpart1) * (distcc - rpart2)    &
     /(distcc *( distcc * (rpart1 + rpart2) - rpart1**2 - rpart2**2))       &
     * (omega1 * log(1+gamma) + omega2 *log(1 - gamma))

end subroutine edlsa

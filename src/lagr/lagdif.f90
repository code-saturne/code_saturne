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

subroutine lagdif &
!================

 ( dx     , vvue   , vpart  , marko  , tempf  , depint ,             &
   dtl    , tstruc , tdiffu , ttotal , vstruc ,                      &
   romp   , taup   , kdif   , tlag2  , lvisq  , yplus  ,             &
   unif1  ,unif2   , dintrf , rpart  ,                               &
   kdifcl ,indint  , gnorm  , vnorm  , grpn   , piiln)

!===============================================================================
! Purpose:
! --------
!
!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------
!
!
!   Deposition sub-model:
!
!   Management of the diffusion phases (marko = 2)
!
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
! depint           ! r  ! <-- ! interface location near-wall/core-flow         !
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

!===============================================================================
!     Modules
!===============================================================================

use paramx
use cstnum
use cstphy
use lagran

!===============================================================================

implicit none

! Arguments

integer          marko , indint

double precision tempf
double precision dx, vvue, vpart, depint
double precision dtl, tstruc, tdiffu, ttotal, vstruc
double precision romp,taup,kdif,tlag2,lvisq,yplus
double precision unif1,unif2,dintrf, rpart
double precision kdifcl, gnorm, vnorm, grpn, piiln

! Local variables

double precision vpart0,vvue0, tci , force
double precision aux1 , aux2 , aux3 , aux4 , aux5 , aux6
double precision aux7 , aux8 , aux9 , aux10, aux11
double precision aa , bb , cc , dd , ee
double precision ter1x , ter2x , ter3x , ter4x , ter5x
double precision ter1f , ter2f , ter3f
double precision ter1p , ter2p , ter3p , ter4p, ter5p
double precision gama2 , omegam , omega2
double precision p11 , p21 , p22 , p31 , p32 , p33
double precision grga2 , gagam , gaome
double precision yplusa,dxaux,dtp1
double precision vagaus(4)

!===============================================================================

call normalen(4,vagaus)

vpart0 = vpart

if (marko.eq.12) then
  vvue0 = vagaus(4) * sqrt(kdif**2 * tlag2 / 2.d0)
else
  vvue0 = vvue
endif

tci   = piiln * tlag2 + vnorm
force =  (ro0*grpn/romp + gnorm) * taup

!---> Coefficients and deterministic terms computation
!     ------------------------------------------------

aux1 = exp( -dtl / taup)
aux2 = exp( -dtl / tlag2 )
aux3 = tlag2 / (tlag2 - taup)

aux4 = tlag2 / (tlag2 + taup)
aux5 = tlag2 * (1.d0-aux2)
aux6 = kdif**2 * tlag2

aux7 = tlag2 - taup
aux8 = kdif**2 * aux3**2

!---> terms for the trajectory

aa = taup * (1.d0 - aux1)
bb = (aux5 - aa) * aux3
cc = dtl - aa - bb

ter1x = aa * vpart0
ter2x = bb * vvue0
ter3x = cc * tci
ter4x = (dtl - aa) * force

!---> vu fluid terms
ter1f = vvue0 * aux2
ter2f = tci * (1.d0-aux2)

!---> particles velocity terms

dd = aux3 * (aux2 - aux1)
ee = 1.d0 - aux1

ter1p = vpart0 * aux1
ter2p = vvue0 * dd
ter3p = tci * (ee-dd)
ter4p = force * ee

!---> Coefficients computation for the stochastic integrals:
!---> Integral on the particles position

gama2  = 0.5d0 * (1.d0 - aux2*aux2 )
omegam = 0.5d0 * aux4 * ( aux5 - aux2*aa )                        &
        -0.5d0 * aux2 * bb
omegam = omegam * sqrt(aux6)

omega2 = aux7                                                     &
       * (aux7*dtl - 2.d0 * (tlag2 * aux5-taup*aa))               &
       + 0.5d0 * tlag2 * tlag2 * aux5                             &
       * (1.d0 + aux2)                                            &
       + 0.5d0 * taup * taup * aa * (1.d0+aux1)                   &
       - 2.d0 * aux4 * tlag2 * taup * taup                        &
       * (1.d0 - aux1*aux2)

omega2 = aux8 * omega2

if (abs(gama2).gt.epzero) then

  p21 = omegam / sqrt(gama2)
  p22 = omega2 - p21**2
  p22 = sqrt( max(zero,p22) )

else
  p21 = 0.d0
  p22 = 0.d0
endif

ter5x = p21 * vagaus(1) + p22 * vagaus(2)

!---> vu fluid integral

p11   = sqrt( gama2*aux6 )
ter3f = p11*vagaus(1)

!---> particles velocity Integral

aux9  = 0.5d0 * tlag2 * (1.d0 - aux2*aux2)
aux10 = 0.5d0 * taup  * (1.d0 - aux1*aux1)
aux11 =  taup * tlag2 * (1.d0 - aux1*aux2)                        &
       /(taup + tlag2)

grga2 = (aux9 - 2.d0*aux11 + aux10) * aux8
gagam = (aux9 - aux11) * (aux8 / aux3)
gaome = ( (tlag2 - taup) * (aux5 - aa)                            &
            - tlag2 * aux9 - taup * aux10                         &
            +(tlag2 + taup) * aux11) * aux8

if(p11.gt.epzero) then
  p31 = gagam / p11
else
  p31 = 0.d0
endif
!
if(p22.gt.epzero) then
  p32 = (gaome-p31*p21) / p22
else
  p32 = 0.d0
endif

p33 = grga2 - p31**2 - p32**2
p33 = sqrt( max(zero,p33) )

ter5p =  p31 * vagaus(1)                                           &
       + p32 * vagaus(2)                                           &
       + p33 * vagaus(3)

!===============================================================================
! 3. Writings finalization
!===============================================================================

!---> trajectory

dx = ter1x + ter2x + ter3x + ter4x + ter5x

!---> vu fluid velocity

vvue = ter1f + ter2f + ter3f

!---> particles velocity

vpart = ter1p + ter2p + ter3p + ter4p + ter5p

yplusa = yplus - dx / lvisq

!---------------------------------------------------
!  Dissociation of cases by the arrival position
!---------------------------------------------------
!
if (yplusa.gt.depint) then

  marko = -2

elseif (yplusa.lt.dintrf) then

  marko = 0
  vvue = sqrt((kdifcl)**2*tlag2/2.d0)                             &
        *sqrt(2.d0*pi)*0.5d0

  dx = dx*(dintrf - yplus)/(yplusa - yplus)

  dxaux = dx

  vpart =  (yplus - yplusa)*lvisq / dtl

  dtp1  = dtl*(dintrf-yplusa)/(yplus-yplusa)
  yplus = dintrf

  call lagdcl                                                     &
  !==========
        ( dx     , vvue   , vpart  , marko  ,                     &
          tempf  , depint ,                                       &
          dtp1   , tstruc , tdiffu , ttotal , vstruc ,            &
          romp   , taup   , kdif   , tlag2  , yplus  ,            &
          lvisq  ,                                                &
          unif1  , unif2  ,dintrf  , rpart  ,                     &
          kdifcl , indint , gnorm  , vnorm  , grpn   , piiln )

  dx = dxaux + dx

else

  if (unif1.lt.(dtl/tdiffu)) then

    if (unif2.lt.0.5d0) then

      marko = 1
      vvue  = vstruc + gnorm*taup + vnorm

    else

      marko = 3
      vvue  = -vstruc + gnorm*taup + vnorm

    endif

  else

    marko = 2

  endif

endif

end subroutine

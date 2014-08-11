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

subroutine lagswe &
!================

 ( dx    , vvue   , vpart  , marko  ,                                 &
   tempf , depint ,                                                   &
   dtp   , tstruc , tdiffu , ttotal , vstruc ,                        &
   romp  , taup   , kdif   , tlag2  , lvisq  , yplus ,                &
   unif1 , unif2  , dintrf , rpart  , kdifcl , gnorm ,                &
   vnorm , grpn   , piiln  )

!===============================================================================
!
! Purpose:
! --------
!
!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------
!
!
!   Deposition submodel:
!
!   Management of the sweep coherent structure (marko = 1)
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! dx               ! r  ! <-> ! wall-normal displacement                       !
! vvue             ! r  ! <-> ! wall-normal velocity of the flow seen          !
! vpart            ! r  ! <-> ! particle wall-normal velocity                  !
! marko            ! i  ! --> ! state of the jump process                      !
! tempf            ! r  ! <-- ! temperature of the fluid                       !
! dtp              ! r  ! <-- ! Lagrangian time step                           !
! tstruc           ! r  ! <-- ! coherent structure mean duration               !
! tdiffu           ! r  ! <-- ! diffusion phase mean duration                  !
! ttotal           ! r  ! <-- ! tdiffu + tstruc                                !
! vstruc           ! r  ! <-- ! coherent structure velocity                    !
! romp             ! r  ! <-- ! particle density                               !
! taup             ! r  ! <-- ! particle relaxation time                       !
! kdif             ! r  ! <-- ! diffusion phase diffusion coefficient          !
! tlag2            ! r  ! <-- ! diffusion relaxation timescale                 !
! lvisq            ! r  ! <-- ! wall-unit lengthscale                          !
! yplus            ! r  ! <-- ! wall-normal velocity of the flow seen          !
! unif1            ! r  ! <-- ! random number (uniform law)                    !
! unif2            ! r  ! <-- ! random number (uniform law)                    !
! dintrf           ! r  ! <-- ! extern-intern interface location               !
! rpart            ! r  ! <-- ! particle radius                                !
! kdifcl           ! r  ! <-- ! internal zone diffusion coefficient            !
! gnorm            ! r  ! <-- ! wall-normal gravity component                  !
! vnorm            ! r  ! <-- ! wall-normal fluid (Eulerian) velocity          !
! grpn             ! r  ! <-- ! wall-normal pressure gradient                  !
! depint           ! r  ! <-- ! interface location near-wall/core-flow         !
! piiln            ! r  ! <-- ! SDE integration auxiliary term                 !
!-------------------------------------------------------------------------------
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array

!===============================================================================
!     Module files
!===============================================================================

use cstnum

!===============================================================================

implicit none

! Arguments

integer          marko

double precision tempf,depint
double precision dx, vvue, vpart
double precision dtp, tstruc, tdiffu, ttotal, vstruc
double precision romp,taup,kdif,tlag2,lvisq,yplus
double precision unif1,unif2,dintrf, rpart
double precision kdifcl, gnorm, vnorm, grpn, piiln

! Local variables

integer          indint
double precision vvue0,vpart0,yplusa,ypluss
double precision dxaux,dtp1

!--------------------------------------------------------
!  Computation phase
!--------------------------------------------------------

vvue0  = vvue
vpart0 = vpart

vvue  =  vstruc + gnorm*taup + vnorm

!  Deposition submodel

vpart =  vpart0*exp(-dtp/taup)                                &
      + (1-exp(-dtp/taup))*vvue0

dx    = vvue0*dtp + vvue0*taup                                &
           *(exp(-dtp/taup)-1)                                &
      + vpart0*taup*(1-exp(-dtp/taup))

yplusa = yplus - dx / lvisq

!---------------------------------------------------------
!    Dissociation of cases by the arrival position
!---------------------------------------------------------

if (yplusa.gt.depint) then

  marko = -2

else if (yplusa.lt.dintrf) then

  dtp1  = (dintrf - yplusa)*lvisq / abs(vpart)
  dx    = dx*(dintrf - yplus)/(yplusa - yplus)
  dxaux = dx

  ypluss = yplus
  yplus  = dintrf
  vvue   = - vstruc + gnorm*taup + vnorm

  marko  = 0
  indint = 1

  call lagdcl &
  !==========
( dx    , vvue   , vpart  , marko  ,                      &
  tempf , depint ,                                        &
  dtp1  , tstruc , tdiffu , ttotal , vstruc ,             &
  romp  , taup   , kdif   , tlag2  , yplus  ,             &
  lvisq ,                                                 &
  unif1 , unif2  , dintrf , rpart  ,                      &
  kdifcl, indint , gnorm  , vnorm  , grpn   , piiln )

  indint = 0

  dx = dx + dxaux

  yplusa = ypluss - dx / lvisq

  if (yplusa.gt.dintrf) then

    marko = 3
    vvue  = - vstruc + gnorm*taup + vnorm

    call lageje &
    !==========
  ( marko ,                                                 &
    tempf , depint ,                                        &
    dtp1  , tstruc , vstruc , lvisq ,                       &
    dx    , vvue   , vpart  , taup  , yplus,                &
    unif1 , unif2  , dintrf , gnorm , vnorm )

    dx = dx + dxaux

  endif

else

   if (unif1.lt.(dtp/tstruc)) then
     marko = 12
   else
     marko = 1
   endif

endif

return
end subroutine

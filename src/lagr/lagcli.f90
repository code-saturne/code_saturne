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

subroutine lagcli &
!================

 ( marko  ,                                                       &
   tempf  , lvisq  , tvisq  ,                                     &
   vpart  , vvue , dx    ,                                        &
   diamp  , romp , taup  , yplus  , dintrf , enertur , gnorm ,    &
   vnorm  , grpn , piiln , depint )

!===============================================================================
! Purpose:
! --------
!
!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------
!   Deposition submodel:
!
!   1/ Parameter initialization
!   2/ Call of the different subroutines with respect to the marko indicator
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! marko            ! i  ! <-> ! state of the jump process                      !
! tempf            ! r  ! <-- ! temperature of the fluid                       !
! lvisq            ! r  ! <-- ! wall-unit lenghtscale                          !
! tvisq            ! r  ! <-- ! wall-unit timescale                            !
! vpart            ! r  ! <-- ! particle wall-normal velocity                  !
! vvue             ! r  ! <-- ! wall-normal velocity of the flow seen          !
! dx               ! r  ! <-- ! wall-normal displacement                       !
! diamp            ! r  ! <-- ! particle diameter                              !
! romp             ! r  ! <-- ! particle density                               !
! taup             ! r  ! <-- ! particle relaxation time                       !
! yplus            ! r  ! <-- ! particle wall-normal normalized distance       !
! dintrf           ! r  ! <-- ! extern-intern interface location               !
! enertur          ! r  ! <-- ! turbulent kinetic energy                       !
! gnorm            ! r  ! <-- ! wall-normal gravity component                  !
! vnorm            ! r  ! <-- ! wall-normal fluid (Eulerian) velocity          !
! grpn             ! r  ! <-- ! wall-normal pressure gradient                  !
! piiln            ! r  ! <-- ! SDE integration auxiliary term                 !
! depint           ! r  ! <-- ! interface location near-wall/core-flow         !
!-------------------------------------------------------------------------------
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran

!===============================================================================

implicit none

integer          marko

double precision tempf
double precision vpart , vvue , vnorm, grpn, piiln
double precision diamp , romp , taup ,yplus , dx   , dintrf, gnorm, depint

! Local variables

 integer          indint

 double precision vstruc,tlag2,kdif,kdifcl,rpart
 double precision tstruc,tdiffu,ttotal, lvisq, tvisq
 double precision unif(2),unif1(1), enertur, ectype, paux
 double precision rapkvp

!===============================================================================
! 1. Initializations
!===============================================================================

! Ratio between k and v'

rapkvp = 0.39d0


! The temporal parameters estimated from the DNS computations
! and written in adimensional form

tlag2 = 3.d0 * tvisq
tstruc = 30.d0 * tvisq
tdiffu = 10.d0 * tvisq
ttotal = tstruc + tdiffu

! The velocity Vstruc is estimated as the square-root of turbulent kinetic
! energy times rapkvp, which corresponds roughly to v' in most of the turbulent boundary
! layer

vstruc = sqrt(enertur * rapkvp)

! From Vstruc we estimate the diffusion coefficient Kdif to balance the fluxes.
! Kdif is roughly equal to sqrt(k/(4*pi)) in the core flow
! (which is the theoretical value of the standard Langevin model with a C0 = 2.1)
! such as:     flux_langevin = sig / sqrt(2*pi)       = v' / sqrt(2*pi)
!                   and (v') = k * C0 /( 1 + 3*C0/2 ) = approx. k/2 (see Minier & Pozorski, 1999)

if (ttotal .gt. (sqrt(pi * rapkvp)*tstruc)) then
  kdif = sqrt(enertur / tlag2) * (ttotal - sqrt(pi * rapkvp)*tstruc) /  tdiffu
else
  write(*,*) "valeur des parametres incorrect dans lagcli"
  call csexit(1)
endif

! Ratios computation of the flux to determine the kdifcl value

ectype = sqrt(kdif**2 * tlag2 / 2.d0)
paux = sqrt(pi / 2.d0) * tstruc * vstruc / (ectype * tdiffu)
paux = paux / (1.d0 + paux)

kdifcl = kdif * (tdiffu / ttotal)

call zufall(2,unif)
!==========

indint = 0

!===============================================================================
! 2. Treatment of the 'degenerated' cases (marko = 10, 20, 30)
!===============================================================================

if (marko.eq.10) then

  marko = 0
  vvue  = 0.d0

else if (marko.eq.20) then

  call zufall(1, unif1(1))
  !==========
  if (unif1(1).lt.paux) then
    marko = 1
  else
    marko = 12
  endif

else if (marko.eq.30) then

  call zufall(1, unif1(1))
  !==========
  if (unif1(1).lt.0.5d0) then
    marko = 1
  else
    marko = 3
  endif

endif

!===============================================================================
! 2. Call of different subroutines given the value of marko
!
! marko = 1 --> sweep phase --> call of lagswe
! marko = 2 or 12 --> diffusion phase --> call of lagdif
! marko = 3 --> ejection phase --> call of lageje
! marko = 0 --> inner-zone diffusion
!
!===============================================================================

rpart = diamp * 0.5d0

if (marko.eq.1) then

  call lagswe &
  !==========
( dx      , vvue   , vpart  , marko  ,                   &
  tempf   , depint ,                                     &
  dtp     , tstruc , tdiffu , ttotal , vstruc ,          &
  romp    , taup   , kdif   , tlag2  , lvisq  , yplus,   &
  unif(1) , unif(2), dintrf , rpart  ,                   &
  kdifcl  , gnorm  , vnorm  , grpn   , piiln )

else if (marko.eq.2 .or. marko.eq.12) then

  call lagdif &
  !==========
( dx     , vvue   , vpart  , marko  ,                    &
  tempf  , depint ,                                      &
  dtp    , tstruc , tdiffu , ttotal , vstruc,            &
  romp   , taup   , kdif   , tlag2  , lvisq , yplus ,    &
  unif(1), unif(2), dintrf , rpart  ,                    &
  kdifcl , indint , gnorm  , vnorm  , grpn  , piiln )

elseif (marko.eq.3) then

  call lageje &
  !==========
( marko  ,                                                  &
  tempf  ,  depint,                                         &
  dtp    , tstruc , vstruc , lvisq ,                        &
  dx     , vvue   , vpart  , taup  , yplus ,                &
  unif(1), unif(2), dintrf , gnorm , vnorm )

elseif (marko.eq.0) then

  call lagdcl &
  !==========
( dx    , vvue   , vpart , marko ,                       &
  tempf , depint ,                                       &
  dtp   , tstruc , tdiffu, ttotal, vstruc,               &
  romp  , taup   , kdif  , tlag2 , yplus ,               &
  lvisq  ,                                               &
  unif(1), unif(2), dintrf, rpart ,                      &
  kdifcl , indint , gnorm , vnorm , grpn  , piiln )

endif

return
end subroutine

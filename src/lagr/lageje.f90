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

subroutine lageje &
!================

 ( marko ,                                                        &
   tempf , depint,                                                &
   dtp   , tstruc , vstruc , lvisq ,                              &
   dx    , vvue   , vpart  , taup  , yplus ,                      &
   unif1 , unif2  , dintrf, gnorm, vnorm)

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
!   Management of the ejection coherent structure (marko = 3)
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! marko            ! i  ! --> ! state of the jump process                      !
! tempf            ! r  ! <-- ! temperature of the fluid                       !
! depint           !  r ! <-- ! interface location near-wall/core-flow         !
! rpart            ! r  ! <-- ! particle radius                                !
! kdifcl           ! r  ! <-- ! internal zone diffusion coefficient            !
! dtp              ! r  ! <-- ! Lagrangian timestep                            !
! tstruc           ! r  ! <-- ! coherent structure mean duration               !
! vstruc           ! r  ! <-- ! coherent structure velocity                    !
! lvisq            ! r  ! <-- ! wall-unit lengthscale                          !
! dx               ! r  ! <-> ! wall-normal displacement                       !
! vpart            ! r  ! <-> ! particle wall-normal velocity                  !
! vvue             ! r  ! <-> ! wall-normal velocity of the flow seen          !
! taup             ! r  ! <-- ! particle relaxation time                       !
! yplus            ! r  ! <-- ! particle wall-normal normalized distance       !
! unif1            ! r  ! <-- ! random number (uniform law)                    !
! unif2            ! r  ! <-- ! random number (uniform law)                    !
! dintrf           ! r  ! <-- ! extern-intern interface location               !
! gnorm            ! r  ! <-- ! wall-normal gravity component                  !
! vnorm            ! r  ! <-- ! wall-normal fluid (Eulerian) velocity          !
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

integer marko

double precision tempf
double precision tstruc, vstruc
double precision dtp, lvisq
double precision vpart  , vvue  , dx

double precision unif1 , unif2 , dintrf, depint
double precision taup  , yplus, gnorm, vnorm

! Local variables

double precision vpart0 , vvue0 , ypaux

!===============================================================================

vvue0  = vvue
vpart0 = vpart

! Gravity and ormal fluid velocity added

vvue   =  -vstruc + gnorm*taup + vnorm

vpart  =  vpart0*exp(-dtp/taup)                                 &
        + (1-exp(-dtp/taup))*vvue0

dx     =  vvue0*dtp + vvue0                                     &
        * taup*(exp(-dtp/taup)-1)                               &
        + vpart0*taup*(1-exp(-dtp/taup))

ypaux = yplus - dx / lvisq

!---------------------------------------------------------
!    Dissociation of cases by the arrival position
!---------------------------------------------------------

if (ypaux.gt.depint) then
  marko = -2
elseif (ypaux.lt.dintrf) then
  marko =  0
else
  if (unif1 .lt. (dtp/tstruc) ) then
    marko = 12
  else
    marko = 3
  endif
endif

return
end subroutine

!-------------------------------------------------------------------------------

!VERS

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
! Purpose:
! -------

!> \file cs_user_initialization-unified_combustion.f90
!> \brief Unified combustion example
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time step)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine cs_user_initialization &
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)

! Local variables

!< [loc_var_dec]
integer          iel, ige, mode, icla, icha
integer          ioxy

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision f1mc(ncharm), f2mc(ncharm)
double precision xkent, xeent, d2s3
double precision wmh2o,wmco2,wmn2,wmo2,dmas

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

! Control Print

write(nfecra,9001)

! Local variables initialization

d2s3 = 2.d0/3.d0

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if (isuite.eq.0) then

! --> Initialisation of k and epsilon (exemple)

  xkent = 1.d-10
  xeent = 1.d-10

! ---- TURBULENCE

  if (itytur.eq.2) then

    do iel = 1, ncel
      rtp(iel,ik)  = xkent
      rtp(iel,iep) = xeent
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      rtp(iel,ir11) = d2s3*xkent
      rtp(iel,ir22) = d2s3*xkent
      rtp(iel,ir33) = d2s3*xkent
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
      rtp(iel,iep)  = xeent
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      rtp(iel,ik)   = xkent
      rtp(iel,iep)  = xeent
      rtp(iel,iphi) = d2s3
      rtp(iel,ifb)  = 0.d0
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      rtp(iel,ik)   = xkent
      rtp(iel,iomg) = xeent/cmu/xkent
    enddo

  elseif (iturb.eq.70) then

    do iel = 1, ncel
      rtp(iel,inusa) = cmu*xkent**2/xeent
    enddo

  endif

! --> All the domain is filled with the first oxidizer at TINITK
!                   ============================================

! ---- Computation of H1INIT and H2INIT

  t1init = t0
  t2init = t0

! ------ Transported variables for the solid phase
!         initialy lacking

  do icla = 1, nclacp
    icha = ichcor(icla)

    do iel = 1, ncel

      rtp(iel,isca(ixch(icla))) = zero
      rtp(iel,isca(ixck(icla))) = zero
      rtp(iel,isca(inp(icla) )) = zero
      rtp(iel,isca(ih2(icla)))  = zero
    enddo
  enddo

! ------ Transported variables for the mix (solid+carrying gas)^2

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo

!       Oxidizer are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
!       the composition of the fisrt oxidiser is taken in account

  coefe(io2) = wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ih2o) = wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ico2) = wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
  mode = -1
  call cs_coal_htconvers1(mode, h1init, coefe, f1mc, f2mc, t1init)
 !============================

  do iel = 1, ncel
    rtp(iel,isca(iscalt)) = h1init
  enddo

! ------ Transported variables for the mix (passive scalars, variance)

  do icha = 1, ncharb
    do iel = 1, ncel
      rtp(iel,isca(if1m(icha))) = zero
      rtp(iel,isca(if2m(icha))) = zero
    enddo
  enddo

  do iel = 1, ncel

    if ( noxyd .ge. 2 ) then
      rtp(iel,isca(if4m)) = zero
    endif
    if ( noxyd .ge. 3 ) then
      rtp(iel,isca(if5m)) = zero
    endif

    if ( ippmod(iccoal) .ge. 1 ) then
      rtp(iel,isca(if6m)) = zero
    endif

    rtp(iel,isca(if7m)) = zero

    if ( ihtco2 .eq. 1 ) then
      rtp(iel,isca(if8m)) = zero
    endif

    if ( ihth2o .eq. 1 ) then
      rtp(iel,isca(if9m)) = zero
    endif

    rtp(iel,isca(ifvp2m)) = zero

    if ( ieqco2.ge.1 ) then

      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)

      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas

      rtp(iel,isca(iyco2)) = oxyco2(ioxy)*wmco2/dmas

    endif

    if ( ieqnox .eq. 1 ) then
      rtp(iel,isca(iyhcn)) = 0.d0
      rtp(iel,isca(iyno )) = 0.d0
      rtp(iel,isca(ihox )) = h1init
    endif

  enddo

endif
!< [init]

!--------
! Formats
!--------

 9001 format(                                                   /,&
'  cs_user_initialization : Variables Initialisation for'      ,/,&
'                    Pulverized Coal by User'                  ,/,&
                                                                /)

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_initialization

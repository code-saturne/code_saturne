!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file cs_user_initialization-pulverized_coal.f90
!>
!> \brief Pulverized coal example
!>
!> See \ref cs_user_initialization for examples.
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

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
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, ige, mode, icha
integer          ioxy

double precision t1init, h1init, coefe(ngazem)
double precision t2init
double precision f1mc(ncharm), f2mc(ncharm)
double precision dmas, wmco2, wmh2o, wmn2, wmo2

double precision, dimension(:), pointer :: cvar_scalt, cvar_yco2
double precision, dimension(:), pointer :: cvar_hox
!< [loc_var_dec]

!===============================================================================

!< [init]
! Variables initialization:
!   ONLY when this is not a restarted computation

if (isuite.gt.0) return

! Control Print

write(nfecra,9001)

! All the domain is filled with the first oxidizer at TINITK
! ==========================================================

! Computation of H1INIT and T2INIT

t1init = t0
t2init = t0

! Transported variables for the mix (solid+carrying gas)^2

do ige = 1, ngazem
  coefe(ige) = zero
enddo

! Oxidizers are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
! the composition of the fisrt oxidiser is taken in account

coefe(io2) =   wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
coefe(ih2o) =   wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
coefe(ico2) =   wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

do icha = 1, ncharm
  f1mc(icha) = zero
  f2mc(icha) = zero
enddo

mode = -1
call cs_coal_htconvers1(mode, h1init, coefe, f1mc, f2mc, t1init)

call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

do iel = 1, ncel
  cvar_scalt(iel) = h1init
enddo

! Transported variables for the mix (passive scalars, variance)
! Variables not present here are initialized to 0.

call field_get_val_s(ivarfl(isca(iyco2)), cvar_yco2)
call field_get_val_s(ivarfl(isca(ihox)), cvar_hox)

if (ieqco2 .ge. 1) then

  ioxy   = 1
  wmo2   = wmole(io2)
  wmco2  = wmole(ico2)
  wmh2o  = wmole(ih2o)
  wmn2   = wmole(in2)
  dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
          +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2)
  xco2 = oxyco2(ioxy)*wmco2/dmas

  do iel = 1, ncel
    cvar_yco2(iel) = oxyco2(ioxy)*wmco2/dmas
  enddo

endif

if (ieqnox .eq. 1) then

  do iel = 1, ncel
    cvar_hox(iel) = h1init
  enddo

endif

! Formats
!--------

 9001 format(                                                   /,&
'  cs_user_initialization: settings for pulverized coal',       /,&
                                                                /)
!< [init]

!----
! End
!----

return
end subroutine cs_user_f_initialization

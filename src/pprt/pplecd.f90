!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!  Function:
!  ---------

!> file pplecd.f90
!>
!> \brief Read specific physical model data file

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!

subroutine pplecd

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ihmpre
use cs_coal_incl
use ppcpfu
use radiat

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================

! ---> Diffusion flame - 3-point chemistry
!      Premix flame    - EBU model
!      Premix flame    - LWC model

if (ippmod(icod3p).ge.0 .or. ippmod(icoebu).ge.0                      &
                        .or. ippmod(icolwc).ge.0) then
  call colecd
endif

! ---> Pulverized coal combustion

if (ippmod(iccoal).ge.0 .or. ippmod(icpl3c).ge.0) then
  call uisofu(iirayo, iihmpr, ncharm, ncharb, nclpch, nclacp,         &
              ncpcmx, ichcor, diam20, cch,                            &
              hch, och, nch, sch, ipci, pcich, cp2ch, rho0ch,         &
              thcdch , cck, hck, ock, nck, sck, xashch,               &
              xashsec, xwatch, h0ashc, cpashc,                        &
              iy1ch, y1ch, iy2ch, y2ch, a1ch, a2ch, e1ch, e2ch,       &
              crepn1, crepn2, ahetch, ehetch, iochet, ahetc2,         &
              ehetc2, ioetc2, ahetwt, ehetwt, ioetwt,                 &
              ieqnox, ieqco2, imdnox, irb, ihtco2, ihth2o, qpr, fn,   &
              ckabs1, noxyd, oxyo2, oxyn2, oxyh2o, oxyco2,            &
              repnck, repnle, repnlo)
  call cs_coal_readata
endif

! ---> Fuel flame

if (ippmod(icfuel).ge.0) then
  call cs_fuel_readata
endif

! ---> Joule effect, electric arc, or ionic conduction

if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
  call ellecd(ippmod(ieljou), ippmod(ielarc))
endif

!----
! End
!----

return
end subroutine


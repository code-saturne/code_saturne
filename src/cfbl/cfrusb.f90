!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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

subroutine cfrusb &
!================

 ( ifac   ,                                                       &
   bc_en  , bc_pr  , bc_vel    )

!===============================================================================
! FUNCTION :
! ---------

! Rusanov flux at the boundary for Euler + Energy

! d rho   /dt + div rho u             = 0
! d rho u /dt + div rho u u + grad  P = 0
! d E     /dt + div rho u E + div u P = 0

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ifac             ! i  ! <-- ! face number                                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe, only:rvoid1
use numvar
use optcal
use cstphy
use cstnum
use parall
use entsor
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

integer          ifac

double precision bc_en(nfabor), bc_pr(nfabor), bc_vel(3,nfabor)

! Local variables

integer          iel, l_size
integer          ien
integer          iflmab

double precision und, uni, rund, runi, cd, ci, cd2(1), ci2(1)
double precision rrus, runb
double precision, dimension(:,:), pointer :: cofacv
double precision, dimension(:), pointer :: coface
double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr, cvar_en, cpro_cp, cpro_cv

double precision cpi(1), cvi(1)

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

ien = isca(ienerg)

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_s(ivarfl(ien), cvar_en)

call field_get_coefac_v(ivarfl(iu), cofacv)
call field_get_coefac_s(ivarfl(ien), coface)

iel   = ifabor(ifac)

! Initialize local specific heat values and handle uniform cases

if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
  cpi(1) = cpro_cp(iel)
else
  cpro_cp => rvoid1
  cpi(1) = 0.d0
endif

if (icv.ge.0) then
  call field_get_val_s(icv, cpro_cv)
  cvi(1) = cpro_cv(iel)
else
  cpro_cv => rvoid1
  cvi(1) = 0.d0
endif

!===============================================================================
! 1. COMPUTE VALUES NEEDED FOR RUSANOV SCHEME
!===============================================================================

und   = (bc_vel(1,ifac)*surfbo(1,ifac)                          &
       + bc_vel(2,ifac)*surfbo(2,ifac)                          &
       + bc_vel(3,ifac)*surfbo(3,ifac))/surfbn(ifac)
uni   = (vel(1,iel)*surfbo(1,ifac)                            &
       + vel(2,iel)*surfbo(2,ifac)                            &
       + vel(3,iel)*surfbo(3,ifac))/surfbn(ifac)
rund  = brom(ifac)*und
runi  = crom(iel)     *uni

l_size = 1


call cs_cf_thermo_c_square(cpi, cvi, &
                           bc_pr(ifac:ifac), brom(ifac:ifac), cd2, l_size)
call cs_cf_thermo_c_square(cpi, cvi, &
                           cvar_pr(iel:iel), crom(iel:iel), ci2, l_size)

cd    = sqrt(cd2(1))
ci    = sqrt(ci2(1))
rrus  = max(abs(und)+cd,abs(uni)+ci)

runb  = 0.5d0*(brom(ifac)*und+crom(iel)*uni)          &
      - 0.5d0*rrus*(brom(ifac)-crom(iel))

!===============================================================================
! 2. CONVECTIVE RUSANOV FLUX
!===============================================================================

! Tag the faces where a Rusanov flux is computed
! The tag will be used in bilsc2 to retrieve the faces where a Rusanov flux
! has to be imposed
icvfli(ifac) = 1

! Momentum flux (the centered pressure contribution is directly taken into account
! in the pressure BC)
cofacv(1,ifac) = suffbn(ifac)*                                                  &
                 0.5d0*( rund*bc_vel(1,ifac) + runi*vel(1,iel)                   &
                         -rrus*(brom(ifac)*bc_vel(1,ifac)                        &
                         -crom(iel)*vel(1,iel)) )

cofacv(2,ifac) = suffbn(ifac)*                                                  &
                 0.5d0*( rund*bc_vel(2,ifac) + runi*vel(2,iel)                   &
                         -rrus*( brom(ifac)*bc_vel(2,ifac)                       &
                         -crom(iel)*vel(2,iel)) )

cofacv(3,ifac) = suffbn(ifac)*                                                  &
                 0.5d0*( rund*bc_vel(3,ifac) + runi*vel(3,iel)                   &
                         -rrus*(brom(ifac)*bc_vel(3,ifac)                        &
                         -crom(iel)*vel(3,iel)) )

! BC for the pressure gradient in the momentum balance
bc_pr(ifac) = 0.5d0 * (bc_pr(ifac) + cvar_pr(iel))

! Total energy flux
coface(ifac) = suffbn(ifac)*                                                    &
               0.5d0*( rund*bc_en(ifac) + runi*cvar_en(iel)                  &
                       +und*bc_pr(ifac) + uni*cvar_pr(iel)                   &
                       -rrus*(brom(ifac)*bc_en(ifac)                         &
                       -crom(iel)*cvar_en(iel)) )

return

end subroutine

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

subroutine cfrusb &
!================

 ( nvar   ,                                                       &
   ifac   ,                                                       &
   gammag ,                                                       &
   bval   )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! ifac             ! i  ! <-- ! face number                                    !
! gammag           ! r  ! <-- ! gamma du gaz                                   !
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
use numvar
use optcal
use cstphy
use cstnum
use parall
use pointe
use entsor
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar
integer          ifac

double precision gammag

double precision bval(nfabor,nvar)

! Local variables

integer          iel
integer          ien
integer          iflmab

double precision und    , uni    , rund   , runi   , cd     , ci
double precision rrus   , runb
double precision, dimension(:), pointer :: bmasfl
double precision, dimension(:,:), pointer :: cofacv
double precision, dimension(:), pointer :: coface
double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_pr, cvar_en

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 0. INITIALISATION
!===============================================================================

ien = isca(ienerg)

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_s(ivarfl(ien), cvar_en)

call field_get_key_int(ivarfl(ien), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

call field_get_coefac_v(ivarfl(iu), cofacv)
call field_get_coefac_s(ivarfl(ien), coface)

iel   = ifabor(ifac)

!===============================================================================
! 1. COMPUTE VALUES NEEDED FOR RUSANOV SCHEME
!===============================================================================

und   = (bval(ifac,iu)*surfbo(1,ifac)                          &
       + bval(ifac,iv)*surfbo(2,ifac)                          &
       + bval(ifac,iw)*surfbo(3,ifac))/surfbn(ifac)
uni   = (vel(1,iel)*surfbo(1,ifac)                            &
       + vel(2,iel)*surfbo(2,ifac)                            &
       + vel(3,iel)*surfbo(3,ifac))/surfbn(ifac)
rund  = brom(ifac)*und
runi  = crom(iel)     *uni
cd    = sqrt(gammag*bval(ifac,ipr)/brom(ifac))
ci    = sqrt(gammag*cvar_pr(iel)/crom(iel))
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

! Mass flux
bmasfl(ifac) = runb*surfbn(ifac)

! Momentum flux (the centered pressure contribution is directly taken into account
! in the pressure BC)
cofacv(1,ifac) = surfbn(ifac)*                                                  &
                 0.5d0*( rund*bval(ifac,iu) + runi*vel(1,iel)                   &
                         -rrus*(brom(ifac)*bval(ifac,iu)                        &
                         -crom(iel)*vel(1,iel)) )

cofacv(2,ifac) = surfbn(ifac)*                                                  &
                 0.5d0*( rund*bval(ifac,iv) + runi*vel(2,iel)                   &
                         -rrus*( brom(ifac)*bval(ifac,iv)                       &
                         -crom(iel)*vel(2,iel)) )

cofacv(3,ifac) = surfbn(ifac)*                                                  &
                 0.5d0*( rund*bval(ifac,iw) + runi*vel(3,iel)                   &
                         -rrus*(brom(ifac)*bval(ifac,iw)                        &
                         -crom(iel)*vel(3,iel)) )

! BC for the pressure gradient in the momentum balance
bval(ifac,ipr) = 0.5d0 * (bval(ifac,ipr) + cvar_pr(iel))

! Total energy flux
coface(ifac) = surfbn(ifac)*                                                    &
               0.5d0*( rund*bval(ifac,ien) + runi*cvar_en(iel)                  &
                       +und*bval(ifac,ipr) + uni*cvar_pr(iel)                   &
                       -rrus*(brom(ifac)*bval(ifac,ien)                         &
                       -crom(iel)*cvar_en(iel)) )

return

end subroutine

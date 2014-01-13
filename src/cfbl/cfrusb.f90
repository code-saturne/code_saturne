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

subroutine cfrusb &
!================

 ( nvar   ,                                                       &
   ifac   ,                                                       &
   gammag ,                                                       &
   rtp    , bval   )

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
! rtp(ncelet,*)    ! tr ! <-- ! variables de calcul au centre des cellules     !
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

double precision rtp(ncelet,*)
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

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

ien = isca(ienerg)

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

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
uni   = (rtp(iel,iu)*surfbo(1,ifac)                            &
       + rtp(iel,iv)*surfbo(2,ifac)                            &
       + rtp(iel,iw)*surfbo(3,ifac))/surfbn(ifac)
rund  = brom(ifac)*und
runi  = crom(iel)     *uni
cd    = sqrt(gammag*bval(ifac,ipr)/brom(ifac))
ci    = sqrt(gammag*rtp(iel,ipr)/crom(iel))
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
                 0.5d0*( rund*bval(ifac,iu) + runi*rtp(iel,iu)                  &
                         -rrus*(brom(ifac)*bval(ifac,iu)                        &
                         -crom(iel)*rtp(iel,iu)) )

cofacv(2,ifac) = surfbn(ifac)*                                                  &
                 0.5d0*( rund*bval(ifac,iv) + runi*rtp(iel,iv)                  &
                         -rrus*( brom(ifac)*bval(ifac,iv)                       &
                         -crom(iel)*rtp(iel,iv)) )

cofacv(3,ifac) = surfbn(ifac)*                                                  &
                 0.5d0*( rund*bval(ifac,iw) + runi*rtp(iel,iw)                  &
                         -rrus*(brom(ifac)*bval(ifac,iw)                        &
                         -crom(iel)*rtp(iel,iw)) )

! BC for the pressure gradient in the momentum balance
bval(ifac,ipr) = 0.5d0 * (bval(ifac,ipr) + rtp(iel,ipr))

! Total energy flux
coface(ifac) = surfbn(ifac)*                                                    &
               0.5d0*( rund*bval(ifac,ien) + runi*rtp(iel,ien)                  &
                       +und*bval(ifac,ipr) + uni*rtp(iel,ipr)                   &
                       -rrus*(brom(ifac)*bval(ifac,ien)                         &
                       -crom(iel)*rtp(iel,ien)) )

return

end subroutine

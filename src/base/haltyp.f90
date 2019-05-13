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

subroutine haltyp &
!================

 ( ivoset )

!===============================================================================
! FONCTION :
! ---------

! TEST DE LA NECESSITE DU VOISINAGE ETENDU, POUR ENVOI AU C
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivoset           ! e  ! <-- ! indicateur d'activation du vois. et.           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use dimens
use paramx
use cstphy
use optcal
use ppppar
use ppthch
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ivoset, imrgrl

!===============================================================================

ivoset = 0

imrgrl = abs(imrgra)
imrgrl = modulo(imrgrl,10)

if (     imrgrl.eq.2 .or. imrgrl.eq.3 &
    .or. imrgrl.eq.5 .or. imrgrl.eq.6 &
    .or. imrgrl.eq.8 .or. imrgrl.eq.9) then
  ivoset = 1
endif

if (iturb.eq.41) ivoset = 1

if (ippmod(iaeros).ge.0) ivoset = 1

if (ippmod(iatmos).ge.0) then
  ivoset = max(ivoset, cs_at_opt_interp_is_p1_proj_needed())
endif

return
end subroutine

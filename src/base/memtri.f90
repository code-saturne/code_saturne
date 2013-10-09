!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine memtri &
!================

 ( nvar   ,                                                       &
   nproce ,                                                       &
   idt    , irtp   , irtpa  , ipropc ,                            &
   ifinra )

!===============================================================================
!  FONCTION
!  --------

!  GESTION MEMOIRE VARIABLES NON GEOMETRIQUES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nproce           ! e  ! <-- ! nombre de prop phy aux centres                 !
! idt              ! e  ! --> ! "pointeur" sur dt                              !
! irtp, irtpa      ! e  ! --> ! "pointeur" sur rtp, rtpa                       !
! ipropc           ! e  ! --> ! "pointeur" sur propce                          !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use cstphy
use numvar
use entsor
use pointe
use albase
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use lagpar
use lagdim
use lagran
use ihmpre
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar
integer          nproce
integer          idt
integer          irtp   , irtpa
integer          ipropc
integer          ifinra

! Local variables

!===============================================================================

! Allocate main real arrays
irtp   = 1
irtpa  = irtp   + ncelet *nvar
ipropc = irtpa  + ncelet *nvar
idt    = ipropc + ncelet *nproce
ifinra = idt    + ncelet

return
end subroutine

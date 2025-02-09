!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

subroutine ppprcl()  &
 bind(C, name='cs_f_ppprcl')

!===============================================================================
! FONCTION :
! --------

!    PREPARATION DU REMPLISSAGE DES CONDITIONS AUX LIMITES

!           AIGUILLAGE SPECIFIQUE AUX PHYSIQUES PARTICULIERES


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppthch
use coincl
use ppincl
use atincl
use cs_c_bindings
use mesh

!===============================================================================

implicit none

! Arguments

! Local variables

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine

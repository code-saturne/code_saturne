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

subroutine cs_user_metal_structures_source_terms                  &
 ( nvar   , nscal  ,                                              &
   ncmast ,                                                       &
   ltmast , itypst , izmast ,                                     &
   svcond , tmet  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppincl
use mesh
use field
use cs_c_bindings
use cs_f_interfaces

use cs_tagms

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncmast
integer          izmet, met_znb

integer          ltmast(ncelet)
integer          izmast(ncelet)

integer          itypst(ncelet,nvar)
double precision svcond(ncelet,nvar)
double precision tmet

! Local variables

! INSERT_VARIABLE_DEFINITIONS_HERE

!===============================================================================

! INSERT_ADDITIONAL_INITIALIZATION_CODE_HERE

!===============================================================================
! Select the cells which are associated to the metal structures volume
! with a function, getcel(), already defined by the sources code.
!===============================================================================

! INSERT_MAIN_CODE_HERE

!===============================================================================
! Define here by the user the values to specify arrays used by the modelling of
! the metal structures condensation.
!                 -----------------------------------------
! with :
!         - itypst(:,ivar) to specify the condensation source term type,
!         - svcond(:,ivar) the scalar value to multiply by the sink term array
!                          of the metal structures condensation model.
!===============================================================================

! INSERT_MAIN_CODE_HERE

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_user_metal_structures_source_terms

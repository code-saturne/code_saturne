!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine usalin
!================

!===============================================================================
!  Purpose :
! --------

! --- User subroutine dedicated to the use of ALE (Arbitrary Lagrangian Eulerian)
!     method :
!
!          Here one defines parameters and input data dedicated to the use ALE
!          method
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "optcal.f90"
include "albase.f90"

!===============================================================================

! Arguments

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
!
!     Here are some examples that can be adapted and changed by Code Saturne
!     users.
!
!
! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Activation of ALE (Arbitrary Lagrangian Eulerian) method
iale = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Number of iterations for fluid initialization. Contrary to ntmabs (for example)
!     nalinf is not an absolute iteration number, meaning that in case of
!     restart calculation nalinf corresponds to the number of iterations
!     for fuid initialization beginning from the first current iteration of
!     the calculation restart. In general nalinf = 0 in that case.

  nalinf = 75

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Maximum number of iterations in case of implicit Fluid Structure Coupling with structural
!     calculations (internal and/or external(i.e. using Code_Aster)). NALIMX = 1, in case of
!     explicit FSI algorithm.

nalimx = 15

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Relative precision of sub-cycling Fluid Structure Coupling algorithm.
epalim = 1.d-5

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END
! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Mesh viscosity modeling (cf. usvima)
!     0 : isotropic
!     1 : orthotropic
iortvm = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! FORMATS
!----

!----
! End
!----

return
end subroutine

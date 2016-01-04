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

!===============================================================================
! function:
! ---------

!> \file prodsc.f90
!>
!> \brief Dot product \c vava = \f$ v_a \cdot v_b \f$
!>
!> The flag \c isqrt can be used to compute the square root of the dot product
!> or the normed residual of two extensive vectors:
!>  - if \c isqrt = 0:
!>    \f[
!>       v_a \cdot v_b = \sum_{\celli =1}^{\ncell} v_{a_\celli} v_{b_\celli}
!>    \f]
!>  - if \c isqrt = 1:
!>    \f[
!>      v_a \cdot v_b
!>      = \sqrt{\sum_{\celli =1}^{\ncell} v_{a_\celli} v_{b_\celli}}
!>    \f]
!>  - if \c isqrt = 10:
!>    \f[
!>      v_a \cdot v_b
!>      = \sum_{\celli =1}^{\ncell}
!>        \dfrac{v_{a_\celli} v_{b_\celli}}{\norm{\vol{\celli}}}
!>    \f]
!>  - if \c isqrt = 11:
!>    \f[
!>      v_a \cdot v_b
!>      = \sqrt{\sum_{\celli =1}^{\ncell}
!>              \dfrac{v_{a_\celli} v_{b_\celli}}{\norm{\vol{\celli}}}}
!>    \f]
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncel          number of cells
!> \param[in]     isqrt         flag:
!>                               - 0 to return the canonic scalar product
!>                               - 1 to return the square root
!>                               - 10 to return the scalar product of extensive
!>                                    vectors
!>                               - 11 to return the square root of the scalar
!>                                    product of extensive vectors
!> \param[in]     va            first vector to multiply
!> \param[in]     vb            second vector to multiply
!> \param[out]    vavb          dot product
!_______________________________________________________________________________

subroutine prodsc &
 ( ncel   , isqrt  , va     , vb     , vavb   )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use mesh, only: volume
use parall

!===============================================================================

implicit none

! Arguments

integer          ncel,isqrt
double precision vavb
double precision va(*), vb(*)

! Local variables

double precision csdot, csres
external         csdot, csres

!===============================================================================

if (isqrt.le.1) then
  vavb = csdot(ncel, va, vb)

  if (irangp.ge.0) call parsom (vavb)

  if (isqrt.eq.1) vavb = sqrt(vavb)

! Compute the residual of extensive vectors
else

  vavb = csres(ncel, volume, va, vb)

  if (isqrt.eq.11) vavb = sqrt(vavb)

endif

!----
! End
!----

return

end subroutine


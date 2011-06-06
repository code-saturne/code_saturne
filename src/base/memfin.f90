!-------------------------------------------------------------------------------

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

subroutine memfin
!================

!===============================================================================
! Purpose:
! -------

! Free memory allocated in Fortran code.
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

!===============================================================================
! Module files
!===============================================================================

use paramx
use parall
use ppincl
use pointe
use optcal
use atincl
use cfpoin
use albase
use cplsat
use vorinc

!===============================================================================

implicit none

! Local variables

!===============================================================================

call finalize_fortran_omp
call finalize_aux_arrays

if (ippmod(iatmos).ge.0) then
  call finalize_meteo
endif

if (ippmod(icompf).ge.0) then
  call finalize_compf
endif

if (iale.eq.1.or.imobil.eq.1) then
  call finalize_ale
endif

if (ncpdct.gt.0) then
  call finalize_kpdc
endif

if (nctsmt.gt.0) then
  call finalize_tsma
endif

if (nfpt1d.gt.0) then
  call finalize_pt1d
endif

if (ivrtex.eq.1) then
  call finalize_vortex
endif

return
end subroutine

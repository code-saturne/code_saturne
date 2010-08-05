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

subroutine usjoin &
!================

 ( )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

! Define (conforming or non-conforming) mesh joinings.

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
use entsor
use parall

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iutile, ii, nbjoin
integer          iwarnj
integer          tml, tmb, tcm, icm, maxsf, maxbrk
double precision fract, plane, mtf, pmf, tmr

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! Get the number of joining operations already defined
! (equal to zero at the moment)

call numjoi(nbjoin)
! ==========

! ---------------
! Main parameters
! ---------------

! The initial tolerance radius associated to each
! vertex is equal to the lenght of the shortest
! incident edge, multiplied by this fraction.

fract = 0.10d0

! When subdividing faces, 2 faces are considered
! coplanar and may be joined if angle between their
! unit normals (in degree) does not exceed this parameter.

plane = 25.0

! associated verbosity level (debug level if >= 3)

iwarnj = 1


! Joining definition

call defjoi(ii, '98 or 99', fract, plane, iwarnj)
!==========


! -------------------
! Advanced parameters
! -------------------

! Use advanced parameters in case of problem during the joining step
! or to get a better mesh quality

! Merge tolerance factor
! Used to locally modify the tolerance associated to each
! vertex BEFORE the merge step.
!   = 0 => no vertex merge;
!   < 1 => vertex merge is more strict. It may increase the number
!          of tolerance reduction and so define smaller subset of
!          vertices to merge together but it can drive to loose
!          intersections;
!   = 1 => no change;
!   > 1 => vertex merge is less strict. The subset of vertices able
!          to be merged together is greater.

mtf = 1.00d0

! Pre-merge factor. This parameter is used to define a limit
! under which two vertices are merged before the merge step.
! Tolerance limit for the pre-merge = pmf * fraction

pmf = 0.10d0

! Tolerance computation mode: tcm
!
!   1: (default) tol = min. edge length related to a vertex * fraction
!   2: tolerance is computed like in mode 1 with in addition, the
!      multiplication by a coef. which is equal to the max sin(e1, e2)
!      where e1 and e2 are two edges sharing the same vertex V for which
!      we want to compute the tolerance
!  11: as 1 but taking into account only the selected faces
!  12: as 2 but taking into account only the selected faces

tcm = 1

! Intersection computation mode: icm
!  1: (default) Original algorithm. Try to clip intersection on extremity
!  2: New intersection algorithm. Avoid to clip intersection on extremity

icm = 1

! Maximum number of equivalence breaks which is
! enabled during the merge step

maxbrk = 500

! Maximum number of sub-faces when splitting a selected face

maxsf = 100

! tml, tmb and tmr are parameters of the searching algorithm for
! face intersections between selected faces (octree structure).
! Useful if there is a memory limitation.

! Tree Max Level: deepest level reachable during the tree building

tml = 30

! Tree Max. Boxes: max. number of bounding boxes (BB) which can be
! linked to a leaf of the tree (not necessary true for the deepest level)

tmb = 25

! Tree Max. Ratio: stop to build the tree structure when
! number of bounding boxes > tmr * number of faces to locate
! Efficient parameter to reduce memory consumption.

tmr = 5.0

! Advanced parameters setup

iutile = 0
if (iutile.eq.1) then
  call setajp(ii, mtf, pmf, tcm, icm, maxbrk, maxsf, tml, tmb, tmr)
  !==========
endif

return
end subroutine

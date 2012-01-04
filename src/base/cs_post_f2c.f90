!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine psteva &
!================

 ( nummai , nomvar , dimvar , ientla , ivarpr , ntcabs , ttcabs , &
   varcel , varfac , varfbo )

!===============================================================================
! Purpose:
! --------

! Write a cell of face located field based on data provided by the
! Fortran layer: encapsulation so as to provide character string lengths.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nummai           ! a  ! <-- ! numero du maillage associe                     !
! nomvar           ! e  ! <-- ! nom de la variable associee                    !
! dimvar           ! e  ! <-- ! 1 pour scalaire, 3 pour vecteur                !
! ientla           ! e  ! <-- ! si vecteur, 1 si valeurs entrelacees           !
!                  !    !     ! (x1, y1, z1, x2, y2, ..., yn, zn),             !
!                  !    !     ! 0 sinon (x1, x2, ...xn, y1, y2, ...)           !
! ivarpr           ! e  ! <-- ! 1 si variable definie sur maillage             !
!                  !    !     ! "parent", 0 si variable restreinte             !
!                  !    !     ! au maillage nummai                             !
! ntcabs           ! e  ! <-- ! numero de pas de temps (-1 pour une            !
!                  !    !     ! variable independante du temps)                !
! ttcabs           ! r  ! <-- ! temps physique associe                         !
! varcel(*)        ! r  ! <-- ! valeurs aux cellules associees                 !
! varfac(*)        ! r  ! <-- ! valeurs aux faces internes associees           !
! varfbo(*)        ! r  ! <-- ! valeurs aux faces de bord associees            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*32     nomvar
integer          nummai, dimvar, ientla, ivarpr, ntcabs
double precision ttcabs, varcel(*), varfac(*), varfbo(*)

! Local variables

integer          lnmvar

!===============================================================================

lnmvar = len(nomvar)

call pstev1 (nummai, nomvar, lnmvar, dimvar, ientla, ivarpr,      &
!==========
             ntcabs, ttcabs, varcel, varfac, varfbo)

return

end subroutine

subroutine pstsnv &
!================

 ( nomvar , nomva2 , nomva3 )

!===============================================================================
! Purpose:
! --------

! Remove character X, x, or 1 from a Fortran character string if the
! compared strings are identical except for the last character, respectively
! Y, y, or 2 and Z, z, or 3.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nomvar           ! s  ! <-- ! name of the first associated variable          !
! nomva2           ! s  ! <-- ! name of the second associated variable         !
! nomva3           ! s  ! <-- ! name of the third associated variable          !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*32     nomvar, nomva2, nomva3

! Local variables

integer          ii, jj
integer          lnmvar, lnmva2, lnmva3

!===============================================================================

lnmvar = len(nomvar)
lnmva2 = len(nomva2)
lnmva3 = len(nomva3)

if ((lnmvar .eq. lnmva2) .and. (lnmvar .eq. lnmva3)) then

  do 10 ii = lnmvar, 1, -1
    if (     nomvar(ii:ii) .ne. ' '                               &
        .or. nomva2(ii:ii) .ne. ' '                               &
        .or. nomva3(ii:ii) .ne. ' ') then
      goto 20
    endif
 10     continue

 20     continue

  if (ii .gt. 1) then

    jj = ii

    ! Handle the case where the next-to-last character changes, such
    ! as with VelocityX1, VelocityX2, ... in case of a calculation
    ! with multiple phases.

    if (      (ii .gt. 2)                                         &
        .and. (nomvar(ii:ii) .eq. nomva2(ii:ii))                  &
        .and. (nomvar(ii:ii) .eq. nomva3(ii:ii))) then
      ii = jj-1
    endif

    ! Remove the character related to the spatial axis

    if (      nomvar(ii:ii) .eq. 'X'                              &
        .and. nomva2(ii:ii) .eq. 'Y'                              &
        .and. nomva3(ii:ii) .eq. 'Z') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. 'x'                         &
             .and. nomva2(ii:ii) .eq. 'y'                         &
             .and. nomva3(ii:ii) .eq. 'z') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. 'U'                         &
             .and. nomva2(ii:ii) .eq. 'V'                         &
             .and. nomva3(ii:ii) .eq. 'W') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. 'u'                         &
             .and. nomva2(ii:ii) .eq. 'v'                         &
             .and. nomva3(ii:ii) .eq. 'w') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. '1'                         &
             .and. nomva2(ii:ii) .eq. '2'                         &
             .and. nomva3(ii:ii) .eq. '3') then
      nomvar(ii:ii) = ' '
    endif

    ! If the next-to last character was removed, the last one must be shifted.

    if (ii .eq. jj+1) then
      nomvar(ii:ii) = nomvar(jj:jj)
      nomvar(jj:jj) = ' '
    endif

  endif

endif

return

end subroutine

subroutine pstmom &
!================

 ( imom   , dtcm )

!===============================================================================
! Purpose:
! --------

! Get cumulative moment from dtcmom.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! imom             ! i  ! <-- ! moment number                                  !
! dtcm             ! r  ! --> ! cumulative moment                              !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use optcal

!===============================================================================

implicit none

! Arguments

integer          imom
double precision dtcm

!===============================================================================

dtcm = dtcmom(imom)

return

end subroutine


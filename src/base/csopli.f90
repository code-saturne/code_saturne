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

subroutine csopli &
!================

 (irkpar, nrkpar, ilogr0, ilogrp)

!===============================================================================
! Purpose:
! -------

!    Initialize log files using Fortran IO.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! irkpar           ! i  ! <-- ! rank if parallel; -1 if sequential             !
! nrkpar           ! i  ! <-- ! number of parallel ranks                       !
! ilogr0           ! i  ! <-- ! log output option for rank 0                   !
!                  !    !     !   0: not redirected                            !
!                  !    !     !   1: redirected to "listing" file              !
! ilogrp           ! i  ! <-- ! log output option for ranks > 0                !
!                  !    !     !   0: not redirected (for debugging)            !
!                  !    !     !   1: redirected to "listing_n*" files          !
!                  !    !     !   2: redirected to /dev/null (suppressed)      !
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
include "entsor.f90"

!===============================================================================

! Arguments

integer          irkpar, nrkpar, ilogr0, ilogrp

! Local variables

character        name*30

!===============================================================================

nfecra = 6  ! default value for Fortran "stdout"

if (irkpar .le. 0) then
  if (ilogr0 .eq. 1) then
    nfecra = 9
    name = 'listing'
  endif
else
  if (ilogrp .eq. 1) then
    nfecra = 9
    if (nrkpar .ge. 10000) then
      write (name,'(a9,i7.4)') 'listing_n', irkpar + 1
    else
      write (name,'(a9,i4.4)') 'listing_n', irkpar + 1
    endif
  else if (ilogrp.eq.2) then
    nfecra = 9
    name = '/dev/null'
  endif
endif

if (nfecra.eq.9) then
   open (file=name, unit=nfecra, form='formatted', status='unknown', err=900)
endif

goto 950

 900  write (0, 999) name
call csexit (1)

 950  continue

#if defined(_CS_LANG_FR)

 999  format(/,                                    &
'Code_Saturne : Erreur d''initialisation :',/,     &
'Impossible d''ouvrir le fichier : ', a, /)

#else

 999  format(/,                                    &
'Code_Saturne: Initialization error:',/,           &
'Impossible to open the file: ', a, /)

#endif

return
end subroutine

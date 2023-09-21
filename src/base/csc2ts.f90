!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file csc2ts.f90
!> \brief Code-code coupling with source terms.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     n_elts        number of local (distant in reverse mode)
!>                              overlapped cells
!> \param[in]     elt_ids       localisation of coupling cells
!> \param[in]     f_id          field index
!> \param[in]     f_dim         field dimension
!> \param[in]     reverse       reverse mode if 1
!> \param[in]     cw1_received  distant variable array
!>                              or explicit part for reverse mode
!> \param[in]     cw2_received  implicit part for reverse mode
!> \param[in,out] crvexp        explicit source term
!> \param[in,out] crvimp        working table for implicit part
!______________________________________________________________________________

subroutine csc2ts &
 ( n_elts , elt_ids , f_id , f_dim , reverse , cw1_received , cw2_received, &
   crvexp , crvimp)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          n_elts
integer          f_id, f_dim
integer          reverse
double precision crvexp(f_dim, ncelet)
double precision crvimp(f_dim, f_dim, ncelet)
double precision cw1_received(f_dim,n_elts)
double precision cw2_received(f_dim,f_dim,n_elts)
integer elt_ids(n_elts)

! Local variables
integer          isou   , jsou
integer          ipt    , iel
double precision xdis   , xloc   , xtau   , rovtau
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cvara_s
double precision, dimension(:,:), pointer :: cvara_v

!----------------------------------------------------------------------------------


!get the array of cells density
call field_get_val_s(icrom, crom) ! Why icrom and not f_id?

xtau = 100.d0*dtref
if (reverse.eq.0) then
  ! For scalars
  if (f_dim.eq.1) then
  ! get the previous scalar array values
    call field_get_val_prev_s(f_id, cvara_s)

    do ipt = 1, n_elts

      ! Get the true localization in the mesh of a given coupled cell
      iel = elt_ids(ipt)

      rovtau = cell_f_vol(iel)*crom(iel)/xtau
      xdis = cw1_received(1, ipt)

      crvexp(1,iel) = crvexp(1,iel) + rovtau*xdis;
      crvimp(1,1,iel) = crvimp(1,1,iel) - rovtau;
    enddo
    ! For vectors and tensors
  else
    call field_get_val_prev_v(f_id, cvara_v)

    do ipt = 1, n_elts

      iel = elt_ids(ipt)
      rovtau = cell_f_vol(iel)*crom(iel)/xtau

      do isou = 1, f_dim
        xdis = cw1_received(isou,ipt)
        crvexp(isou,iel) = crvexp(isou,iel) + rovtau*xdis
        crvimp(isou,isou,iel) = crvimp(isou, isou,iel) - rovtau
      enddo

    enddo
  endif
endif

if (reverse.eq.1) then
  ! For scalars
  if (f_dim.eq.1) then
  !get the previous scalar array values
    call field_get_val_prev_s(f_id, cvara_s)

    do ipt = 1, n_elts

      !get the true localization in the mesh of a given coupled cell
      iel = elt_ids(ipt)

      crvexp(1,iel) = crvexp(1,iel) + cw1_received(1,ipt)/xtau
      crvimp(1,1,iel) = crvimp(1,1,iel) - cw2_received(1,1,ipt)/xtau

    enddo


    ! For vectors and tensors
  else
    call field_get_val_prev_v(f_id, cvara_v)

    do ipt = 1, n_elts

      iel = elt_ids(ipt)

      do isou = 1, f_dim
        crvexp(isou,iel) = crvexp(isou,iel) + cw1_received(isou,ipt)/xtau
        do jsou = 1, f_dim
          crvimp(jsou,isou,iel) = crvimp(jsou,isou,iel) &
                                - cw2_received(jsou,isou,ipt)/xtau
        enddo
      enddo

    enddo
  endif
endif

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine

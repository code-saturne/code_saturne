!-------------------------------------------------------------------------------

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

!> \file post_util.f90

!===============================================================================
! Function:
! ---------

!> \brief Compute thermal flux at boundary.

!> If working with enthalpy, compute an enthalpy flux.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    bflux         boundary heat flux at selected faces
!_______________________________________________________________________________

subroutine post_boundary_thermal_flux &
 ( nfbrps , lstfbr ,                                              &
   bflux )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use optcal
use numvar
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                :: nfbrps
integer, dimension(nfbrps), intent(in)             :: lstfbr
double precision, dimension(nfbrps), intent(out)   :: bflux

! Local variables

integer ::         f_id
integer ::         iloc, ivar

character(len=80) :: f_name
integer(c_int), dimension(:), allocatable :: c_faces

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_post_boundary_flux(scalar_name, n_loc_b_faces, b_face_ids,   &
                                   b_face_flux)                              &
    bind(C, name='cs_post_boundary_flux')
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=1), dimension(*), intent(in) :: scalar_name
    integer(c_int), value :: n_loc_b_faces
    integer(c_int), dimension(*) :: b_face_ids
    real(kind=c_double), dimension(*) :: b_face_flux
  end subroutine cs_post_boundary_flux

end interface

!===============================================================================

! Initialize variables to avoid compiler warnings

if (iscalt.gt.0) then

  ivar = isca(iscalt)
  f_id = ivarfl(ivar)

  call field_get_name (f_id, f_name)

  allocate(c_faces(nfbrps))

  do iloc = 1, nfbrps
    c_faces(iloc) = lstfbr(iloc) - 1
  enddo

  call cs_post_boundary_flux(trim(f_name)//c_null_char, nfbrps, c_faces, bflux)

  deallocate(c_faces)

else ! if thermal variable is not available

  do iloc = 1, nfbrps
    bflux(iloc) = 0.d0
  enddo

endif

!----
! End
!----

return
end subroutine post_boundary_thermal_flux

!===============================================================================
! Function:
! ---------

!> \brief Compute stress at boundary.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    stress        stress at selected faces
!_______________________________________________________________________________

subroutine post_stress &
 ( nfbrps , lstfbr ,                                                    &
   stress )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use optcal
use numvar
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                 :: nfbrps
integer, dimension(nfbrps), intent(in)              :: lstfbr
double precision, dimension(3, nfbrps), intent(out) :: stress

! Local variables

integer          :: ifac  , iloc
double precision :: srfbn
double precision, dimension(:,:), pointer :: forbr

!===============================================================================

call field_get_val_v(iforbr, forbr)

do iloc = 1, nfbrps
  ifac = lstfbr(iloc)
  srfbn = surfbn(ifac)
  stress(1,iloc) = forbr(1,ifac)/srfbn
  stress(2,iloc) = forbr(2,ifac)/srfbn
  stress(3,iloc) = forbr(3,ifac)/srfbn
enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_stress

!===============================================================================
! Function:
! ---------

!> \brief Extract stress normal to the boundary.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    effnrm        stress normal to wall at selected faces
!_______________________________________________________________________________

subroutine post_stress_normal &
 ( nfbrps , lstfbr ,                                              &
   effnrm )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use optcal
use numvar
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                 :: nfbrps
integer, dimension(nfbrps), intent(in)              :: lstfbr
double precision, dimension(nfbrps), intent(out)    :: effnrm

! Local variables

integer                        :: ifac  , iloc
double precision               :: srfbn
double precision, dimension(3) :: srfnor
double precision, dimension(:,:), pointer :: forbr

!===============================================================================

call field_get_val_v(iforbr, forbr)

do iloc = 1, nfbrps
  ifac = lstfbr(iloc)
  srfbn = surfbn(ifac)
  srfnor(1) = surfbo(1,ifac) / srfbn
  srfnor(2) = surfbo(2,ifac) / srfbn
  srfnor(3) = surfbo(3,ifac) / srfbn
  effnrm(iloc) =  (  forbr(1,ifac)*srfnor(1)                                 &
                   + forbr(2,ifac)*srfnor(2)                                 &
                   + forbr(3,ifac)*srfnor(3)) / srfbn
enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_stress_normal

!===============================================================================
! Function:
! ---------

!> \brief Compute tangential stress at boundary.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    stress        stress at selected faces
!_______________________________________________________________________________

subroutine post_stress_tangential &
 ( nfbrps , lstfbr ,                                              &
   stress )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use optcal
use numvar
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                 :: nfbrps
integer, dimension(nfbrps), intent(in)              :: lstfbr
double precision, dimension(3, nfbrps), intent(out) :: stress

! Local variables

integer                        :: ifac  , iloc
double precision               :: srfbn, fornor
double precision, dimension(3) :: srfnor
double precision, dimension(:,:), pointer :: forbr

!===============================================================================

call field_get_val_v(iforbr, forbr)

do iloc = 1, nfbrps
  ifac = lstfbr(iloc)
  srfbn = surfbn(ifac)
  srfnor(1) = surfbo(1,ifac) / srfbn
  srfnor(2) = surfbo(2,ifac) / srfbn
  srfnor(3) = surfbo(3,ifac) / srfbn
  fornor =    forbr(1,ifac)*srfnor(1)                                 &
            + forbr(2,ifac)*srfnor(2)                                 &
            + forbr(3,ifac)*srfnor(3)
  stress(1,iloc) = (forbr(1,ifac) - fornor*srfnor(1)) / srfbn
  stress(2,iloc) = (forbr(2,ifac) - fornor*srfnor(2)) / srfbn
  stress(3,iloc) = (forbr(3,ifac) - fornor*srfnor(3)) / srfbn
enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_stress_tangential

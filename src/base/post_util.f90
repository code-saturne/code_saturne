!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

integer ::         inc, iccocg, f_id
integer ::         ifac, iloc, ivar
integer ::         iel
integer ::         iflmab

double precision :: cpp   , srfbn , heq
double precision :: flumab, diipbx, diipby, diipbz

double precision, allocatable, dimension(:) :: theipb
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: bmasfl, tscalp, hextp, hintp, dist_theipb

type(var_cal_opt) :: vcopt

logical(c_bool), dimension(:), pointer ::  cpl_faces

!===============================================================================

! Initialize variables to avoid compiler warnings

if (iscalt.gt.0) then

  ivar   = isca(iscalt)

  f_id = ivarfl(ivar)

  ! Boundary condition pointers for gradients and advection

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)

  ! Boundary condition pointers for diffusion

  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! Boundary condition pointers for diffusion with coupling

  call field_get_hext(f_id, hextp)
  call field_get_hint(f_id, hintp)

  ! Pointers to fields and properties

  call field_get_val_prev_s(ivarfl(ivar), tscalp)

  if (iscacp(iscalt).eq.1 .and. icp.ge.0) then
    call field_get_val_s(icp, cpro_cp)
  endif

  call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
  call field_get_val_s(iflmab, bmasfl)

  ! Compute variable values at boundary faces

  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

  allocate(theipb(nfabor))
  theipb = 0.d0

  do iloc = 1, nfbrps
    ifac = lstfbr(iloc)
    iel = ifabor(ifac)
    theipb(ifac) = tscalp(iel)
  enddo

  ! Reconstructed fluxes
  if (vcopt%ircflu .gt. 0 .and. itbrrb.eq.1) then

    ! Compute gradient of temperature / enthalpy

    allocate(grad(3,ncelet))

    inc = 1
    iccocg = 1

    call field_gradient_scalar &
      (ivarfl(ivar), 0, imrgra, inc,         &
       iccocg,                               &
       grad)

    ! Compute reconstructed temperature

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)

      theipb(ifac) = theipb(ifac)       &
                   + diipbx*grad(1,iel) &
                   + diipby*grad(2,iel) &
                   + diipbz*grad(3,iel)
    enddo

    deallocate(grad)
  endif

  if (vcopt%icoupl.gt.0) then
    call field_get_coupled_faces(ivarfl(ivar), cpl_faces)
    allocate(dist_theipb(nfabor))
    call cs_ic_field_dist_data_by_face_id(ivarfl(ivar), 1, theipb, dist_theipb)
  endif

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    iel = ifabor(ifac)

    if (iscacp(iscalt).eq.1) then
      if (icp.ge.0) then
        cpp = cpro_cp(iel)
      else
        cpp = cp0
      endif
    else
      cpp = 1.d0
    endif

    srfbn = max(surfbn(ifac), epzero**2)
    flumab = bmasfl(ifac)

    bflux(iloc) = (cofafp(ifac) + cofbfp(ifac)*theipb(ifac))   &
                - flumab*cpp/srfbn*(coefap(ifac) + coefbp(ifac)*theipb(ifac))
    ! here bflux = 0 if current face is coupled

    if (vcopt%icoupl.gt.0) then
      if (cpl_faces(ifac)) then
        heq = hextp(ifac) * hintp(ifac) / ((hextp(ifac) + hintp(ifac))*srfbn)
        bflux(iloc) = heq*(theipb(ifac)-dist_theipb(ifac))
      endif
    endif

  enddo

  if (vcopt%icoupl.gt.0) then
    deallocate(dist_theipb)
  endif

  deallocate(theipb)

else ! if thermal variable is not available

  do iloc = 1, nfbrps
    bflux(iloc) = 0.d0
  enddo

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_boundary_thermal_flux

!===============================================================================
! Function:
! ---------

!> \brief Compute Nusselt number near boundary.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    bnussl        Nusselt near boundary
!_______________________________________________________________________________

subroutine post_boundary_nusselt &
 ( nfbrps , lstfbr ,                                              &
   bnussl )

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
double precision, dimension(nfbrps), intent(out)   :: bnussl

! Local variables

integer ::         inc, iccocg
integer ::         iel, ifac, iloc, ivar
integer ::         ifcvsl, itplus, itstar

double precision :: xvsl  , srfbn , heq
double precision :: diipbx, diipby, diipbz
double precision :: numer, denom

double precision, allocatable, dimension(:) :: theipb
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer :: cviscl
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: tplusp, tstarp
double precision, dimension(:), pointer :: tscalp, hextp, hintp, dist_theipb

type(var_cal_opt) :: vcopt

logical(c_bool), dimension(:), pointer ::  cpl_faces

!===============================================================================

! pointers to T+ and T* if saved

call field_get_id_try('tplus', itplus)
call field_get_id_try('tstar', itstar)

if (itstar.ge.0 .and. itplus.ge.0) then

  ivar = isca(iscalt)

  call field_get_val_prev_s(ivarfl(ivar), tscalp)

  call field_get_val_s(itplus, tplusp)
  call field_get_val_s(itstar, tstarp)

  ! Boundary condition pointers for diffusion

  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! Boundary condition pointers for diffusion with coupling

  call field_get_hext(ivarfl(ivar), hextp)
  call field_get_hint(ivarfl(ivar), hintp)

  ! Compute variable values at boundary faces

  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

  allocate(theipb(nfabor))
  theipb = 0.d0

  do iloc = 1, nfbrps
    ifac = lstfbr(iloc)
    iel = ifabor(ifac)
    theipb(ifac) = tscalp(iel)
  enddo

  ! Reconstructed fluxes
  if (vcopt%ircflu .gt. 0 .and. itbrrb.eq.1) then

    ! Compute gradient of temperature / enthalpy

    allocate(grad(3,ncelet))

    inc = 1
    iccocg = 1

    call field_gradient_scalar &
      (ivarfl(ivar), 0, imrgra, inc,         &
       iccocg,                               &
       grad)

    ! Compute reconstructed temperature

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)

      theipb(ifac) = theipb(ifac)                                               &
                   + diipbx*grad(1,iel) + diipby*grad(2,iel) + diipbz*grad(3,iel)
    enddo

    deallocate(grad)
  endif

  if (vcopt%icoupl.gt.0) then
    call field_get_coupled_faces(ivarfl(ivar), cpl_faces)
    allocate(dist_theipb(nfabor))
    call cs_ic_field_dist_data_by_face_id(ivarfl(ivar), 1, theipb, dist_theipb)
  endif

  ! Physical property pointers

  call field_get_key_int (ivarfl(ivar), kivisl, ifcvsl)
  if (ifcvsl .ge. 0) then
    call field_get_val_s(ifcvsl, cviscl)
  endif

  ! Boundary condition pointers for gradients and advection

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)

  ! Compute using reconstructed temperature value in boundary cells

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    iel = ifabor(ifac)

    if (ifcvsl.ge.0) then
      xvsl = cviscl(iel)
    else
      xvsl = visls0(iscalt)
    endif

    srfbn = surfbn(ifac)

    numer = (cofafp(ifac) + cofbfp(ifac)*theipb(ifac)) * distb(ifac)
    ! here numer = 0 if current face is coupled

    if (vcopt%icoupl.gt.0) then
      if (cpl_faces(ifac)) then
        heq = hextp(ifac) * hintp(ifac) / ((hextp(ifac) + hintp(ifac))*srfbn)
        numer = heq*(theipb(ifac)-dist_theipb(ifac)) * distb(ifac)
      endif
    endif

    denom = xvsl * tplusp(ifac)*tstarp(ifac)

    if (abs(denom).gt.1e-30) then
      bnussl(iloc) = numer / denom
    else
      bnussl(iloc) = 0.d0
    endif

  enddo

  if (vcopt%icoupl.gt.0) then
    deallocate(dist_theipb)
  endif

  deallocate(theipb)

else ! default if not computable

  do iloc = 1, nfbrps
    bnussl(iloc) = -1.d0
  enddo

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_boundary_nusselt

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

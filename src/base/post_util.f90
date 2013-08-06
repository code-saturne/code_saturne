!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
!> \param[in]     rtp           calculated variables at cell centers
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[out]    bflux         boundary heat flux at selected faces
!_______________________________________________________________________________

subroutine post_boundary_thermal_flux &
 ( nfbrps , lstfbr ,                                              &
   rtp    , propce , propfb ,                                     &
   bflux )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                        :: nfbrps
integer, dimension(nfbrps), intent(in)                     :: lstfbr
double precision, dimension(ncelet, *), intent(in), target :: rtp, propce
double precision, dimension(ndimfb, *), intent(in)         :: propfb
double precision, dimension(nfbrps), intent(out)           :: bflux

! Local variables

integer ::         inc, iccocg, nswrgp, imligp, iwarnp
integer ::         ifac, iloc, ivar
integer ::         iel
integer ::         ipcvsl, ipcvst, iflmab

double precision :: xcp   , xvsl  , srfbn
double precision :: visct , flumab, diipbx, diipby, diipbz, tcel
double precision :: epsrgp, climgp, extrap

double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer :: bmasfl

!===============================================================================

! Initialize variables to avoid compiler warnings

if (iscalt.gt.0) then

  ivar   = isca(iscalt)

  ! Boundary condition pointers for gradients and advection

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)

  ! Boundary condition pointers for diffusion

  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! Pointers to properties

  if (ivisls(iscalt).gt.0) then
    ipcvsl = ipproc(ivisls(iscalt))
  else
    ipcvsl = 0
  endif
  ipcvst = ipproc(ivisct)
  call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
  call field_get_val_s(iflmab, bmasfl)

  ! Compute variable values at boundary faces

  if (ircflu(ivar) .gt. 0) then

    ! Compute gradient of temperature / enthalpy

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
    endif

    ! Allocate a temporary array for the gradient calculation
    allocate(grad(ncelet,3))

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel &
    !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap , rtp(1,ivar) , coefap , coefbp ,     &
   grad   )

    ! Compute diffusive and convective flux using reconstructed temperature

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)

      tcel =   rtp(iel,ivar)                                                  &
             + diipbx*grad(iel,1) + diipby*grad(iel,2) + diipbz*grad(iel,3)

      if (ipcvsl.gt.0) then
        xvsl = propce(iel,ipcvsl)
      else
        xvsl = visls0(iscalt)
      endif
      srfbn = max(surfbn(ifac), epzero**2)
      visct  = propce(iel,ipcvst)
      flumab = bmasfl(ifac)

      bflux(iloc) =                (cofafp(ifac) + cofbfp(ifac)*tcel)   &
                    - flumab/srfbn*(coefap(ifac) + coefbp(ifac)*tcel)

    enddo

    deallocate(grad)

  else ! If flux is not reconstructed

    ! Compute diffusive and convective flux using non-reconstructed temperature

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      tcel = rtp(iel,ivar)

      if (ipcvsl.gt.0) then
        xvsl = propce(iel,ipcvsl)
      else
        xvsl = visls0(iscalt)
      endif
      srfbn = max(surfbn(ifac), epzero**2)
      visct  = propce(iel,ipcvst)
      flumab = bmasfl(ifac)

      bflux(iloc) =                (cofafp(ifac) + cofbfp(ifac)*tcel)   &
                    - flumab/srfbn*(coefap(ifac) + coefbp(ifac)*tcel)

    enddo

  endif ! test on reconstruction

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

!> \brief Compute temperature at boundary.

!> If working with enthalpy, compute an enthalpy.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[in]     rtp           calculated variables at cell centers
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[out]    btemp         boundary temperature at selected faces
!_______________________________________________________________________________

subroutine post_boundary_temperature &
 ( nfbrps , lstfbr ,                                              &
   rtp    , propce , propfb ,                                     &
   btemp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                        :: nfbrps
integer, dimension(nfbrps), intent(in)                     :: lstfbr
double precision, dimension(ncelet, *), intent(in), target :: rtp, propce
double precision, dimension(ndimfb, *), intent(in)         :: propfb
double precision, dimension(nfbrps), intent(out)           :: btemp

! Local variables

integer ::         inc, iccocg, nswrgp, imligp, iwarnp
integer ::         iel, ifac, iloc, ivar
integer ::         itplus, itstar

double precision :: diipbx, diipby, diipbz
double precision :: epsrgp, climgp, extrap, tcel

double precision, dimension(:), pointer :: coefap, coefbp
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer :: tplusp, tstarp

!===============================================================================

! pointers to T+ and T* if saved

call field_get_id('tplus', itplus)
call field_get_id('tstar', itstar)

if (itstar.ge.0 .and. itplus.ge.0) then

  call field_get_val_s (itplus, tplusp)
  call field_get_val_s (itstar, tstarp)

  ivar   = isca(iscalt)

  ! Compute variable values at boundary faces

  if (ircflu(ivar) .gt. 0) then

    ! Boundary condition pointers for gradients and advection

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)

    ! Compute gradient of temperature / enthalpy

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
    endif

    ! Allocate a temporary array for the gradient calculation
    allocate(grad(ncelet,3))

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel &
    !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap , rtp(1,ivar) , coefap , coefbp ,     &
   grad   )

    ! Compute reconstructed value in boundary cells

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      tcel =   rtp(iel,ivar)                                                 &
             + diipbx*grad(iel,1) + diipby*grad(iel,2) + diipbz*grad(iel,3)

      btemp(iloc) = tcel - tplusp(ifac)*tstarp(ifac)

    enddo

    deallocate(grad)

  else ! If flux is not reconstructed

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      tcel = rtp(iel,ivar)

      btemp(iloc) = tcel - tplusp(ifac)*tstarp(ifac)

    enddo

  endif

else ! default if not computable

  do iloc = 1, nfbrps
    btemp(iloc) = -1.d0
  enddo

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_boundary_temperature

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
!> \param[in]     rtp           calculated variables at cell centers
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[out]    bnussl        Nusselt near boundary
!_______________________________________________________________________________

subroutine post_boundary_nusselt &
 ( nfbrps , lstfbr ,                                              &
   rtp    , propce , propfb ,                                     &
   bnussl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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

!===============================================================================

implicit none

! Arguments

integer, intent(in)                                        :: nfbrps
integer, dimension(nfbrps), intent(in)                     :: lstfbr
double precision, dimension(ncelet, *), intent(in), target :: rtp, propce
double precision, dimension(ndimfb, *), intent(in)         :: propfb
double precision, dimension(nfbrps), intent(out)           :: bnussl

! Local variables

integer ::         inc, iccocg, nswrgp, imligp, iwarnp
integer ::         iel, ifac, iloc, ivar
integer ::         ipcvsl, ipcvst, itplus, itstar

double precision :: xvsl  , srfbn
double precision :: visct , flumab, diipbx, diipby, diipbz
double precision :: epsrgp, climgp, extrap, numer, denom, tcel

double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:), pointer :: tplusp, tstarp

!===============================================================================

! pointers to T+ and T* if saved

call field_get_id('tplus', itplus)
call field_get_id('tstar', itstar)

if (itstar.ge.0 .and. itplus.ge.0) then

  call field_get_val_s (itplus, tplusp)
  call field_get_val_s (itstar, tstarp)

  ivar   = isca(iscalt)

  ! Boundary condition pointers for diffusion

  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! Physical property pointers

  if (ivisls(iscalt).gt.0) then
    ipcvsl = ipproc(ivisls(iscalt))
  else
    ipcvsl = 0
  endif
  ipcvst = ipproc(ivisct)

  ! Compute variable values at boundary faces

  if (ircflu(ivar) .gt. 0) then

    ! Boundary condition pointers for gradients and advection

    call field_get_coefa_s(ivarfl(ivar), coefap)
    call field_get_coefb_s(ivarfl(ivar), coefbp)

    ! Compute gradient of temperature / enthalpy

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
    endif

    ! Allocate a temporary array for the gradient calculation
    allocate(grad(ncelet,3))

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    call grdcel &
    !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap , rtp(1,ivar) , coefap , coefbp ,     &
   grad   )

    ! Compute using reconstructed temperature value in boundary cells

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      tcel =   rtp(iel,ivar)                                                 &
             + diipbx*grad(iel,1) + diipby*grad(iel,2) + diipbz*grad(iel,3)

      if (ipcvsl.gt.0) then
        xvsl = propce(iel,ipcvsl)
      else
        xvsl = visls0(iscalt)
      endif
      srfbn = max(surfbn(ifac), epzero**2)
      visct = propce(iel,ipcvst)

      numer = (cofafp(ifac) + cofbfp(ifac)*tcel) * distb(ifac)
      denom = xvsl * tplusp(ifac)*tstarp(ifac)

      if (abs(denom).gt.1e-30) then
        bnussl(iloc) = numer / denom
      else
        bnussl(iloc) = 0.d0
      endif

    enddo

    deallocate(grad)

  else ! If flux is not reconstructed

    ! Compute using non-reconstructed temperature value in boundary cells

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)
      iel = ifabor(ifac)

      tcel =   rtp(iel,ivar)

      if (ipcvsl.gt.0) then
        xvsl = propce(iel,ipcvsl)
      else
        xvsl = visls0(iscalt)
      endif
      srfbn = max(surfbn(ifac), epzero**2)
      visct = propce(iel,ipcvst)

      numer = (cofafp(ifac) + cofbfp(ifac)*tcel) * distb(ifac)
      denom = xvsl * tplusp(ifac)*tstarp(ifac)

      if (abs(denom).gt.1e-30) then
        bnussl(iloc) = numer / denom
      else
        bnussl(iloc) = 0.d0
      endif

    enddo

  endif

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

!> \brief Compute efforts at boundary     .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    effort        efforts at selected faces
!_______________________________________________________________________________

subroutine post_efforts &
 ( nfbrps , lstfbr ,                                                    &
   effort )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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
double precision, dimension(3, nfbrps), intent(out) :: effort

! Local variables

integer          :: ifac  , iloc
double precision :: srfbn

!===============================================================================

do iloc = 1, nfbrps
  ifac = lstfbr(iloc)
  srfbn = surfbn(ifac)
  effort(1,iloc) = forbr(1,ifac)/srfbn
  effort(2,iloc) = forbr(2,ifac)/srfbn
  effort(3,iloc) = forbr(3,ifac)/srfbn
enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_efforts

!===============================================================================
! Function:
! ---------

!> \brief Extract efforts normal to the boundary     .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    effnrm        efforts normal to wall at selected faces
!_______________________________________________________________________________

subroutine post_efforts_normal &
 ( nfbrps , lstfbr ,                                              &
   effnrm )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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

!===============================================================================

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
end subroutine post_efforts_normal

!===============================================================================
! Function:
! ---------

!> \brief Compute tangential efforts at boundary     .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfbrps        number of boundary faces to postprocess
!> \param[in]     lstfbr        list of boundary faces to postprocess
!> \param[out]    effort        efforts at selected faces
!_______________________________________________________________________________

subroutine post_efforts_tangential &
 ( nfbrps , lstfbr ,                                              &
   effort )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
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
double precision, dimension(3, nfbrps), intent(out) :: effort

! Local variables

integer                        :: ifac  , iloc
double precision               :: srfbn, fornor
double precision, dimension(3) :: srfnor

!===============================================================================

do iloc = 1, nfbrps
  ifac = lstfbr(iloc)
  srfbn = surfbn(ifac)
  srfnor(1) = surfbo(1,ifac) / srfbn
  srfnor(2) = surfbo(2,ifac) / srfbn
  srfnor(3) = surfbo(3,ifac) / srfbn
  fornor =    forbr(1,ifac)*srfnor(1)                                 &
            + forbr(2,ifac)*srfnor(2)                                 &
            + forbr(3,ifac)*srfnor(3)
  effort(1,iloc) = (forbr(1,ifac) - fornor*srfnor(1)) / srfbn
  effort(2,iloc) = (forbr(2,ifac) - fornor*srfnor(2)) / srfbn
  effort(3,iloc) = (forbr(3,ifac) - fornor*srfnor(3)) / srfbn
enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine post_efforts_tangential


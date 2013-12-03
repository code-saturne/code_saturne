!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

! This is an example of cs_user_extra_operations.f90 which
! performs an energy balance.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     nbpmax        max. number of particles allowed
!> \param[in]     nvp           number of particle-defined variables
!> \param[in]     nvep          number of real particle properties
!> \param[in]     nivep         number of integer particle properties
!> \param[in]     ntersl        number of return coupling source terms
!> \param[in]     nvlsta        number of Lagrangian statistical variables
!> \param[in]     nvisbr        number of boundary statistics
!> \param[in]     itepa         integer particle attributes
!>                                (containing cell, ...)
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     propfa        physical properties at interior face centers
!> \param[in]     propfb        physical properties at boundary face centers
!> \param[in]     ettp, ettpa   particle-defined variables
!> \param[in]                    (at current and previous time steps)
!> \param[in]     tepa          real particle properties
!> \param[in]                    (statistical weight, ...
!> \param[in]     statis        statistic means
!> \param[in]     stativ        accumulator for variance of volume statisitics
!> \param[in]     tslagr        Lagrangian return coupling term
!> \param[in]                    on carrier phase
!> \param[in]     parbor        particle interaction properties
!> \param[in]                    on boundary faces
!_______________________________________________________________________________


subroutine cs_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   itepa  ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagpar
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta), stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)


! Local variables

!< [loc_var_dec]
integer          iel    , ifac   , ivar
integer          iel1   , iel2   , ieltsm
integer          iortho
integer          inc    , iccocg
integer          nswrgp , imligp , iwarnp
integer          ipcrom , ipcvst , iflmas , iflmab , ipccp, ipcvsl
integer          iscal
integer          ncesmp
integer          ilelt  , nlelt

double precision xrtpa  , xrtp
double precision xbilan , xbilvl , xbilpa , xbilpt
double precision xbilsy , xbilen , xbilso , xbildv
double precision xbilmi , xbilma
double precision epsrgp , climgp , extrap
double precision xfluxf , xgamma
double precision diipbx, diipby, diipbz, distbr
double precision visct, flumab , xcp , xvsl, ctb1, ctb2

integer, allocatable, dimension(:) :: lstelt

double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: treco
!< [loc_var_dec]

!===============================================================================
! Initialization
!===============================================================================

!< [init]
! Allocate a temporary array for cells or interior/boundary faces selection
allocate(lstelt(max(ncel,nfac,nfabor)))
!< [init]

!===============================================================================

!< [example_1]
! The balance is not valid if inpdt0=1
if (inpdt0.eq.0) then

  ! 2.1 Initialization
  ! ==================

  ! --> Local variables
  !     ---------------

  ! xbilvl: volume contribution of unsteady terms
  ! xbildv: volume contribution due to to term in div(rho u)
  ! xbilpa: contribution from adiabatic walls
  ! xbilpt: contribution from walls with fixed temperature
  ! xbilsy: contribution from symmetry boundaries
  ! xbilen: contribution from inlets
  ! xbilso: contribution from outlets
  ! xbilmi: contribution from mass injections
  ! xbilma: constribution from mass suctions
  ! xbilan: total balance

  xbilvl = 0.d0
  xbildv = 0.d0
  xbilpa = 0.d0
  xbilpt = 0.d0
  xbilsy = 0.d0
  xbilen = 0.d0
  xbilso = 0.d0
  xbilmi = 0.d0
  xbilma = 0.d0
  xbilan = 0.d0

  iscal = iscalt         ! temperature scalar number
  ivar =  isca(iscal)    ! temperature variable number

  ! Physical quantity numbers
  ipcrom = ipproc(irom)
  ipcvst = ipproc(ivisct)
  iflmas = ipprof(ifluma(ivar))
  iflmab = ipprob(ifluma(ivar))

  ! We save in ipccp a flag allowing to determine if the specific heat is
  ! constant (= cp0) or variable. It will be used to compute balances
  ! (xbilvl is in Joules).
  if (icp.gt.0) then
    ipccp  = ipproc(icp   )
  else
    ipccp  = 0
  endif

  ! We save in ipcvsl a flag allowing to determine if the diffusivity is
  ! constant (= visls0) or variable. It will be used for diffusive terms.
  if (ivisls(iscal).gt.0) then
    ipcvsl = ipproc(ivisls(iscal))
  else
    ipcvsl = 0
  endif

  ! Boundary condition pointers for gradients and advection

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)

  ! Boundary condition pointers for diffusion

  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! --> Synchronization of Cp and Dt
  !     ----------------------------

  ! To compute fluxes at interior faces, it is necessary to have access
  ! to variables at neighboring cells. Notably, it is necessary to know
  ! the specific heat and the time step value. For this,

  ! - in parallel calculations, it is necessary on faces at sub-domain
  !   boundaries to know the value of these variables in cells from the
  !   neighboring subdomain.
  ! - in periodic calculations, it is necessary at periodic faces to know
  !   the value of these variables in matching periodic cells.

  ! To ensure that these values are up to date, it is necessary to use
  ! the synchronization routines to update parallel and periodic ghost
  ! values for Cp and Dt before computing the gradient.

  ! If the calculation is neither parallel nor periodic, the calls may be
  ! kept, as tests on iperio and irangp ensure generality).

  ! Parallel and periodic update

  if (irangp.ge.0.or.iperio.eq.1) then

    ! update Dt
    call synsca(dt)
    !==========

    ! update Cp if variable (otherwise cp0 is used)
    if (ipccp.gt.0) then
      call synsca(propce(1,ipccp))
      !==========
    endif

  endif

  ! --> Compute value reconstructed at I' for boundary faces

  allocate(treco(nfabor))

  ! For non-orthogonal meshes, it must be equal to the value at the
  ! cell center, which is computed in:
  ! treco(ifac) (with ifac=1, nfabor)

  ! For orthogonal meshes, it is sufficient to assign:
  ! rtp(iel, ivar) to treco(ifac), with iel=ifabor(ifac)
  ! (this option corresponds to the second branch of the test below,
  ! with iortho different from 0).

  iortho = 0

  ! --> General case (for non-orthogonal meshes)

  if (iortho.eq.0) then

    ! Allocate a work array for the gradient calculation
    allocate(grad(ncelet,3))

    ! --- Compute temperature gradient

    ! To compute the temperature gradient in a given cell, it is necessary
    ! to have access to values at neighboring cells.  For this,

    ! - in parallel calculations, it is necessary at cells on sub-domain
    !   boundaries to know the value of these variables in cells from the
    !   neighboring subdomain.
    ! - in periodic calculations, it is necessary at cells on periodic
    !   boundaries to know the value of these variables in matching
    !   periodic cells.

    ! To ensure that these values are up to date, it is necessary to use
    ! the synchronization routines to update parallel and periodic ghost
    ! values for the temperature before computing the gradient.

    ! If the calculation is neither parallel nor periodic, the calls may be
    ! kept, as tests on iperio and irangp ensure generality).

    ! - Parallel and periodic update

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
      !==========
    endif


    ! - Compute gradient

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
        epsrgp , climgp , extrap ,                                     &
        rtp(1,ivar) , coefap , coefbp ,                                &
        grad   )

    ! - Compute reconstructed value in boundary cells

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      treco(ifac) =   rtp(iel,ivar)            &
                    + diipbx*grad(iel,1)  &
                    + diipby*grad(iel,2)  &
                    + diipbz*grad(iel,3)
    enddo

    ! Free memory
    deallocate(grad)

  ! --> Case of orthogonal meshes

  else

    ! Compute reconstructed value
    ! (here, we assign the non-reconstructed value)

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      treco(ifac) = rtp(iel,ivar)
    enddo

  endif

  ! 2.1 Compute the balance at time step n
  ! ======================================

  ! --> Balance on interior volumes
  !     ---------------------------

  ! If it is variable, the density 'rom' has been computed at the beginning
  ! of the time step using the temperature from the previous time step.

  if (ipccp.gt.0) then
    do iel = 1, ncel
      xrtpa = rtpa(iel,ivar)
      xrtp  = rtp (iel,ivar)
      xbilvl =   xbilvl                                                &
               + volume(iel) * propce(iel,ipccp) * propce(iel,ipcrom)  &
                                                 * (xrtpa - xrtp)
    enddo
  else
    do iel = 1, ncel
      xrtpa = rtpa(iel,ivar)
      xrtp  = rtp (iel,ivar)
      xbilvl =   xbilvl  &
               + volume(iel) * cp0 * propce(iel,ipcrom) * (xrtpa - xrtp)
    enddo
  endif

  ! --> Balance on all faces (interior and boundary), for div(rho u)
  !     ------------------------------------------------------------

  ! Caution: values of Cp and Dt in cells adjacent to interior faces are
  !          used, which implies having synchronized these values for
  !          parallelism and periodicity.

  ! Note that if Cp is variable, writing a balance on the temperature
  ! equation is not absolutely correct.

  if (ipccp.gt.0) then
    do ifac = 1, nfac

      iel1 = ifacel(1,ifac)
      if (iel1.le.ncel) then
        ctb1 = propfa(ifac,iflmas)*propce(iel1,ipccp)*rtp(iel1,ivar)
      else
        ctb1 = 0d0
      endif

      iel2 = ifacel(2,ifac)
      if (iel2.le.ncel) then
        ctb2 = propfa(ifac,iflmas)*propce(iel2,ipccp)*rtp(iel2,ivar)
      else
        ctb2 = 0d0
      endif

      xbildv =  xbildv + (dt(iel1)*ctb1 - dt(iel2)*ctb2)
    enddo

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      xbildv = xbildv + dt(iel) * propce(iel,ipccp)    &
                                * propfb(ifac,iflmab)  &
                                * rtp(iel,ivar)
    enddo

  ! --- if Cp is constant

  else
    do ifac = 1, nfac

      iel1 = ifacel(1,ifac)
      if (iel1.le.ncel) then
        ctb1 = propfa(ifac,iflmas)*cp0*rtp(iel1,ivar)
      else
        ctb1 = 0d0
      endif

      iel2 = ifacel(2,ifac)
      if (iel2.le.ncel) then
        ctb2 = propfa(ifac,iflmas)*cp0*rtp(iel2,ivar)
      else
        ctb2 = 0d0
      endif

      xbildv = xbildv + (dt(iel1) + dt(iel2))*0.5d0*(ctb1 - ctb2)
    enddo

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      xbildv = xbildv + dt(iel) * cp0                  &
                                * propfb(ifac,iflmab)  &
                                * rtp(iel,ivar)
    enddo
  endif

  ! In case of a mass source term, add contribution from Gamma*Tn+1

  ncesmp = ncetsm
  if (ncesmp.gt.0) then
    do ieltsm = 1, ncesmp
      iel = icetsm(ieltsm)
      xrtp  = rtp (iel,ivar)
      xgamma = smacel(ieltsm,ipr)
      if (ipccp.gt.0) then
        xbildv =   xbildv                                     &
                 - volume(iel) * propce(iel,ipccp) * dt(iel)  &
                               * xgamma * xrtp
      else
        xbildv =   xbildv  &
                 - volume(iel) * cp0 * dt(iel) * xgamma * xrtp
      endif
    enddo
  endif

  ! --> Balance on boundary faces
  !     -------------------------

  ! We handle different types of boundary faces separately to better
  ! analyze the information, but this is not mandatory.

  ! - Compute the contribution from walls with colors 2, 4, and 7
  !   (adiabatic here, so flux should be 0)

  call getfbr('2 or 4 or 7', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Geometric variables

    distbr = distb(ifac)

    ! Physical variables

    visct  = propce(iel,ipcvst)
    flumab = propfb(ifac,iflmab)

    if (ipccp.gt.0) then
      xcp = propce(iel,ipccp)
    else
      xcp    = cp0
    endif

    if (ipcvsl.gt.0) then
      xvsl = propce(iel,ipcvsl)
    else
      xvsl = visls0(iscal)
    endif

    ! Contribution to flux from the current face
    ! (diffusion and convection flux, negative if incoming)

    xfluxf = - surfbn(ifac) * dt(iel) * xcp                 &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp                         &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilpa = xbilpa + xfluxf

  enddo

  ! Contribution from walls with color 6
  ! (here at fixed temperature; the convective flux should be 0)

  call getfbr('6', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Geometric variables

    distbr = distb(ifac)

    ! Physical variables

    visct  = propce(iel,ipcvst)
    flumab = propfb(ifac,iflmab)

    if (ipccp.gt.0) then
      xcp = propce(iel,ipccp)
    else
      xcp    = cp0
    endif

    if (ipcvsl.gt.0) then
      xvsl = propce(iel,ipcvsl)
    else
      xvsl = visls0(iscal)
    endif

    ! Contribution to flux from the current face
    ! (diffusion and convection flux, negative if incoming)

    xfluxf = - surfbn(ifac) * dt(iel) * xcp                 &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp                         &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilpt = xbilpt + xfluxf

  enddo

  ! Contribution from symmetries (should be 0).

  call getfbr('1', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Geometric variables

    distbr = distb(ifac)

    ! Physical variables

    visct  = propce(iel,ipcvst)
    flumab = propfb(ifac,iflmab)

    if (ipccp.gt.0) then
      xcp = propce(iel,ipccp)
    else
      xcp    = cp0
    endif

    if (ipcvsl.gt.0) then
      xvsl = propce(iel,ipcvsl)
    else
      xvsl = visls0(iscal)
    endif

    ! Contribution to flux from the current face
    ! (diffusion and convection flux, negative if incoming)

    xfluxf = - surfbn(ifac) * dt(iel) * xcp                 &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp                         &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilsy = xbilsy + xfluxf

  enddo

  ! Contribution from inlet (color 3, diffusion and convection flux)

  call getfbr('3', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Geometric variables

    distbr = distb(ifac)

    ! Physical variables

    visct  = propce(iel,ipcvst)
    flumab = propfb(ifac,iflmab)

    if (ipccp.gt.0) then
      xcp = propce(iel,ipccp)
    else
      xcp    = cp0
    endif

    if (ipcvsl.gt.0) then
      xvsl = propce(iel,ipcvsl)
    else
      xvsl = visls0(iscal)
    endif

    ! Contribution to flux from the current face
    ! (diffusion and convection flux, negative if incoming)

    xfluxf = - surfbn(ifac) * dt(iel) * xcp                 &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp                         &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilen = xbilen + xfluxf

  enddo

  ! Contribution from outlet (color 5, diffusion and convection flux)

  call getfbr('5', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Geometric variables

    distbr = distb(ifac)

    ! Physical variables

    visct  = propce(iel,ipcvst)
    flumab = propfb(ifac,iflmab)

    if (ipccp.gt.0) then
      xcp = propce(iel,ipccp)
    else
      xcp    = cp0
    endif

    if (ipcvsl.gt.0) then
      xvsl = propce(iel,ipcvsl)
    else
      xvsl = visls0(iscal)
    endif

    ! Contribution to flux from the current face
    ! (diffusion and convection flux, negative if incoming)

    xfluxf = - surfbn(ifac) * dt(iel) * xcp                 &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp                         &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilso = xbilso + xfluxf

  enddo

  ! Now the work array for the temperature can be freed
  deallocate(treco)


  ! --> Balance on mass source terms
  !     ----------------------------

  ! We separate mass injections from suctions for better generality

  ncesmp = ncetsm
  if (ncesmp.gt.0) then
    do ieltsm = 1, ncesmp
      ! depending on the type of injection we use the 'smacell' value
      ! or the ambient temperature
      iel = icetsm(ieltsm)
      xgamma = smacel(ieltsm,ipr)
      if (itypsm(ieltsm,ivar).eq.0 .or. xgamma.lt.0.d0) then
        xrtp = rtp (iel,ivar)
      else
        xrtp = smacel(ieltsm,ivar)
      endif
      if (ipccp.gt.0) then
        if (xgamma.lt.0.d0) then
          xbilma =   xbilma  &
                   + volume(iel) * propce(iel,ipccp) * dt(iel) * xgamma * xrtp
        else
          xbilmi =   xbilmi  &
                   + volume(iel) * propce(iel,ipccp) * dt(iel) * xgamma * xrtp
        endif
      else
        if (xgamma.lt.0.d0) then
          xbilma =   xbilma  &
                   + volume(iel) * cp0 * dt(iel) * xgamma * xrtp
        else
          xbilmi =   xbilmi  &
                   + volume(iel) * cp0 * dt(iel) * xgamma * xrtp
        endif
      endif
    enddo
  endif

  ! Sum of values on all ranks (parallel calculations)

  if (irangp.ge.0) then
    call parsom(xbilvl)
    call parsom(xbildv)
    call parsom(xbilpa)
    call parsom(xbilpt)
    call parsom(xbilsy)
    call parsom(xbilen)
    call parsom(xbilso)
    call parsom(xbilmi)
    call parsom(xbilma)
  endif

  ! --> Total balance
  !     -------------

  ! We add the different contributions calculated above.

  xbilan =   xbilvl + xbildv + xbilpa + xbilpt + xbilsy + xbilen   &
           + xbilso + xbilmi + xbilma

  ! 2.3 Write the balance at time step n
  ! ====================================

  write (nfecra, 2000)                                               &
    ntcabs, xbilvl, xbildv, xbilpa, xbilpt, xbilsy, xbilen, xbilso,  &
    xbilmi, xbilma, xbilan

2000 format                                                           &
  (/,                                                                 &
   3X,'** Thermal balance **', /,                                     &
   3X,'   ---------------', /,                                        &
   '---', '------',                                                   &
   '------------------------------------------------------------', /, &
   'bt ','  Iter',                                                    &
   '   Volume     Divergence  Adia Wall   Fixed_T Wall  Symmetry',    &
   '      Inlet       Outlet  Inj. Mass.  Suc. Mass.  Total', /,      &
   'bt ', i6, 10e12.4, /,                                             &
   '---','------',                                                    &
   '------------------------------------------------------------')

endif ! End of test on inpdt0
!< [example_1]

!< [finalize]
! Deallocate the temporary array
deallocate(lstelt)
!< [finalize]

return
end subroutine cs_user_extra_operations

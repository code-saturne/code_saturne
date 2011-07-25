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

subroutine usproj &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   itepa  ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Called at end of each time step, very general purpose
!    (i.e. anything that does not have another dedicated user subroutine)


! Several examples are given here:

!  - compute a thermal balance
!    (if needed, see note  below on adapting this to any scalar)

!  - compute global efforts on a subset of faces

!  - arbitrarily modify a calculation variable

!  - extract a 1 d profile

!  - print a moment

!  - examples on using parallel utility functions

! These examples are valid when using periodicity (iperio .gt. 0)
! and in parallel (irangp .ge. 0).

! The thermal balance compution also illustates a few other features,
! including the required precautions in parallel or with periodicity):
! - gradient calculation
! - computation of a value depending on cells adjacent to a face
!   (see synchronization of Dt and Cp)
! - computation of a global sum in parallel (parsom)


! Cells, boundary faces and interior faces identification
! =======================================================

! Cells, boundary faces and interior faces may be identified using
! the subroutines 'getcel', 'getfbr' and 'getfac' (respectively).

!  getfbr(string, nelts, eltlst):
!  - string is a user-supplied character string containing selection criteria;
!  - nelts is set by the subroutine. It is an integer value corresponding to
!    the number of boundary faces verifying the selection criteria;
!  - lstelt is set by the subroutine. It is an integer array of size nelts
!    containing the list of boundary faces verifying the selection criteria.

!  string may contain:
!  - references to colors (ex.: 1, 8, 26, ...)
!  - references to groups (ex.: inlet, group1, ...)
!  - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!  These criteria may be combined using logical operators ('and', 'or') and
!  parentheses.
!  Example: '1 and (group2 or group3) and y < 1' will select boundary faces
!  of color 1, belonging to groups 'group2' or 'group3' and with face center
!  coordinate y less than 1.

! Similarly, interior faces and cells can be identified using the 'getfac'
! and 'getcel' subroutines (respectively). Their syntax are identical to
! 'getfbr' syntax.

! For a more thorough description of the criteria syntax, it can be referred
! to the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! max. number of particles allowed               !
! nvp              ! i  ! <-- ! number of particle-defined variables           !
! nvep             ! i  ! <-- ! number of real particle properties             !
! nivep            ! i  ! <-- ! number of integer particle properties          !
! ntersl           ! i  ! <-- ! number of return coupling source terms         !
! nvlsta           ! i  ! <-- ! number of Lagrangian statistical variables     !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! itepa            ! ia ! <-- ! integer particle attributes                    !
!  (nbpmax, nivep) !    !     !   (containing cell, ...)                       !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp, ettpa      ! ra ! <-- ! particle-defined variables                     !
!  (nbpmax, nvp)   !    !     !  (at current and previous time steps)          !
! tepa             ! ra ! <-- ! real particle properties                       !
!  (nbpmax, nvep)  !    !     !  (statistical weight, ...                      !
! statis           ! ra ! <-- ! statistic means                                !
!  (ncelet, nvlsta)!    !     !                                                !
! stativ(ncelet,   ! ra ! <-- ! accumulator for variance of volume statisitics !
!        nvlsta -1)!    !     !                                                !
! tslagr           ! ra ! <-- ! Lagrangian return coupling term                !
!  (ncelet, ntersl)!    !     !  on carrier phase                              !
! parbor           ! ra ! <-- ! particle interaction properties                !
!  (nfabor, nvisbr)!    !     !  on boundary faces                             !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta), stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)


! Local variables

integer          iel    , ielg   , ifac   , ifacg  , ivar
integer          iel1   , iel2   , ieltsm
integer          iortho , impout
integer          inc    , iccocg
integer          nswrgp , imligp , iwarnp
integer          iutile , iclvar , iii
integer          ipcrom , ipcvst , iflmas , iflmab , ipccp, ipcvsl
integer          iscal
integer          ii     , nbr    , irangv , irang1 , npoint
integer          imom   , ipcmom , idtcm
integer          itab(3), iun
integer          ncesmp
integer          ilelt  , nlelt

double precision xrtpa  , xrtp
double precision xbilan , xbilvl , xbilpa , xbilpt
double precision xbilsy , xbilen , xbilso , xbildv
double precision xbilmi , xbilma
double precision epsrgp , climgp , extrap
double precision xfluxf , xgamma
double precision diipbx, diipby, diipbz, distbr
double precision visct, flumab , xcp , xvsl, rrr
double precision xfor(3), xyz(3), xabs, xu, xv, xw, xk, xeps

integer, allocatable, dimension(:) :: lstelt

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: treco

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Initialization
!===============================================================================

! Allocate a temporary array for cells or interior/boundary faces selection
allocate(lstelt(max(ncel,nfac,nfabor)))

! Memory management


!===============================================================================
! 2. Example: compute energy balance relative to temperature
!    -------------------------------------------------------

! We assume that we want to compute balances  (convective and diffusive)
! at the boundaries of the calculation domain represented below
! (with boundaries marked by colors).

! The scalar considered if the temperature. We will also use the
! specific heat (to btain balances in Joules)


! Domain and associated boundary colors
! -------------------------------------
!                  6
!      --------------------------
!      |                        |
!      |                        |
!   7  |           1            | 5
!      |     ^                  |
!      |     |                  |
!      --------------------------

!         2  3             4

! 2, 4, 7 : adiabatic walls
! 6       : wall with fixed temperature
! 3       : inlet
! 5       : outlet
! 1       : symmetry

!-------------------------------------------------------------------------------

! To ensure calculations have physical meaning, it is best to use
! a spatially uniform time step (idtvar = 0 or 1).
! In addition, when restarting a calculation, the balance is
! incorrect if inpdt0 = 1 (visct not initialized and t(n-1) not known)

!-------------------------------------------------------------------------------

! Temperature variable: ivar = isca(iscalt) (use rtp(iel, ivar))

!-------------------------------------------------------------------------------

! The balance at time step n is equal to:

!        n        iel=ncelet           n-1
! balance  =   sum { volume(iel)*cp*rom(iel)*(rtpa(iel,ivar)-rtp(iel,ivar)) }
!                 iel=1

!                 ifac=nfabor
!            + sum {
!                 ifac=1

!                     surfbn(ifac)*dt(ifabor(ifac))*cp
!                   * [visls0(iscalt) + visct(ifabor(ifac))/sigmas(iscalt) ]
!                   / distbr(ifac)
!                   * [  coefa(ifac,iclvar)
!                      + (coefb(ifac,iclvar)-1.d0)*rtp(ifabor(ifac,ivar))]
!                  }

!                 ifac=nfabor
!            + sum {
!                 ifac=1
!                     dt(ifabor(ifac))*cp
!                   * rtp(ifabor(ifac,ivar))*(-flumab(ifac))
!                  }

! The first term is negative if the amount of energy in the volume
! has decreased (it is 0 in a steady regime).

! The other terms (convection, diffusion) are positive if the amount
! of energy in the volume has increased due to boundary conditions.

! In a steady regime, a positive balance thus indicates an energy gain.

!-------------------------------------------------------------------------------

! With 'rom' calculated using the density law from the usphyv subroutine,
! for example:

!    n-1
! rom(iel) = p0 / [rr * (rtpa(iel,ivar) + tkelv)]

!-------------------------------------------------------------------------------

! Cp and lambda/Cp may be variable

!-------------------------------------------------------------------------------

! Adaptation to an arbitrary scalar
! ---------------------------------

! The approach may be used for the balance of any other scalar (but the
! balances are not in Joules and the specific heat is not used)

! In this case:

! - replace iscalt by the number iscal of the required scalar,
!   iscal having an allowed range of 1 to nscal.

! - set ipccp to 0 independently of the value of icp and use 1 instead of cp0

!===============================================================================

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
  ivar =  isca(iscal)           ! temperature variable number
  iclvar = iclrtp(ivar,icoef)   ! boundary condition number

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

    call grdcel                                                     &
    !==========
      ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
        iwarnp , nfecra ,                                              &
        epsrgp , climgp , extrap ,                                     &
        rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
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
      iel2 = ifacel(2,ifac)
      xbildv =   xbildv + propfa(ifac,iflmas)                 &
               * (dt(iel1)*propce(iel1,ipccp)*rtp(iel1,ivar)  &
               - dt(iel2)*propce(iel2,ipccp)*rtp(iel2,ivar))
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
      iel2 = ifacel(2,ifac)
      xbildv = xbildv +   (dt(iel1)+ dt(iel2))*0.5d0         &
                        * cp0                                &
                        * propfa(ifac,iflmas)                &
                        * (rtp(iel1,ivar) - rtp(iel2,ivar))
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

    xfluxf =      surfbn(ifac) * dt(iel) * xcp                      &
                * (xvsl+visct/sigmas(iscal)) / distbr               &
                * (  coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)   &
                   * treco(ifac))                                   &
              -   flumab * dt(iel) * xcp                            &
                * (  coefa(ifac,iclvar) + coefb(ifac,iclvar)        &
                   * treco(ifac))

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

    xfluxf =    surfbn(ifac) * dt(iel) * xcp                      &
              * (xvsl+visct/sigmas(iscal)) / distbr               &
              * (  coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)   &
                 * treco(ifac))                                   &
              -   flumab * dt(iel) * xcp                          &
                * (  coefa(ifac,iclvar) + coefb(ifac,iclvar)       &
                   * treco(ifac))

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

    xfluxf =          surfbn(ifac) * dt(iel) * xcp *              &
     (xvsl+visct/sigmas(iscal))/distbr *                          &
     (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)*treco(ifac))   &
                    - flumab * dt(iel) * xcp *                    &
     (coefa(ifac,iclvar)+ coefb(ifac,iclvar)*treco(ifac))

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

    xfluxf =    surfbn(ifac) * dt(iel) * xcp                      &
              * (xvsl+visct/sigmas(iscal))/distbr                 &
              * (  coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)   &
                 * treco(ifac))                                   &
              -   flumab * dt(iel) * xcp                          &
                * (  coefa(ifac,iclvar)+ coefb(ifac,iclvar)       &
                   * treco(ifac))

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

    xfluxf =     surfbn(ifac) * dt(iel) * xcp                      &
               * (xvsl+visct/sigmas(iscal))/distbr                 &
               * (  coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)   &
                  * treco(ifac))                                   &
             -   flumab * dt(iel) * xcp                            &
               * (  coefa(ifac,iclvar)+ coefb(ifac,iclvar)         &
                  * treco(ifac))

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

!===============================================================================
! 3. Example: compute global efforts on a subset of faces
!    ----------------------------------------------------
!===============================================================================

! ----------------------------------------------

! The test below allows checking that the following example compiles
! while disabling it by default.

iutile = 0

if (iutile.eq.0) return

! ----------------------------------------------

! If efforts have been calculated correctly:

if (ineedf.eq.1) then

  do ii = 1, ndim
    xfor(ii) = 0.d0
  enddo

  call getfbr('2 or 3', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    do ii = 1, ndim
      xfor(ii) = xfor(ii) + forbr(ii, ifac)
    enddo

  enddo

  if (irangp.ge.0) then
    call parrsm(ndim,xfor)
  endif

endif

!===============================================================================
! 4. Example: set temperature to 20 in a given region starting at t = 12s
!    --------------------------------------------------------------------

! Do this with precaution...
! The user is responsible for the validity of results.
!===============================================================================

! ----------------------------------------------

! The test below allows checking that the following example compiles
! while disabling it by default.

iutile = 0

if (iutile.eq.0) return

! ----------------------------------------------

iscal = iscalt

if (ttcabs .ge. 12.d0) then

  if (iscal.gt.0 .and. iscal.le.nscal) then
    do iel = 1, ncel
      rtp(iel,isca(iscal)) = 20.d0
    enddo
  endif

  write(nfecra,3000)

endif

 3000 format                                                       &
  (/,                                                              &
   ' User modification of variables at the end of the time step',  &
   /)

!===============================================================================
! 5. Example: extraction of a 1D profile
!    -----------------------------------

! We seek here to extract the profile of U, V, W, k and epsilon on an
! arbitrary 1D curve based on a curvilear abscissa.
! The profile is described in the 'profile.dat' file (do not forget to
! define it as user data in the run script).

! - the curve used here is the segment: [(0;0;0),(0;0.1;0)], but the
!   generalization to an arbitrary curve is simple.
! - the routine handles parallelism an periodicity, as well as the different
!   turbulence models.
! - the 1D curve is discretized into 'npoint' points. For each of these
!   points, we search for the closest cell center and we output the variable
!   values at this cell center. For better consistency, the coordinate
!   which is output is that of the cell center (instead of the initial point).
! - we avoid using the same cell multiple times (in case several points
!   an the curve are associated with the same cell).
!===============================================================================

! ----------------------------------------------

! The test below allows checking that the following example compiles
! while disabling it by default.

iutile = 0

if (iutile.eq.0) return

! ----------------------------------------------

if (ntcabs.eq.ntmabs) then

  ! Only process of rank 0 (parallel) or -1 (scalar) writes to this file.
  ! We use 'user' Fortran units.
  impout = impusr(1)
  if (irangp.le.0) then
    open(impout,file='profile.dat')
    write(impout,*)  &
         '# z(m) U(m/s) V(m/s) W(m/s) k(m2/s2) eps(m2/s3)'
  endif

  npoint = 200
  iel1   = -999
  irang1 = -999
  do ii = 1, npoint

    xyz(1) = 0.d0
    xyz(2) = float(ii-1)/float(npoint-1)*0.1d0
    xyz(3) = 0.d0

    call findpt(ncelet, ncel, xyzcen, xyz(1), xyz(2), xyz(3), iel, irangv)
    !==========

    if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
      iel1   = iel
      irang1 = irangv

      ! Set temporary variables xu, xv, ... for the process containing
      ! the point and then send it to other processes.
      if (irangp.eq.irangv) then
        xabs = xyzcen(2,iel)
        xu   = rtp(iel,iu)
        xv   = rtp(iel,iv)
        xw   = rtp(iel,iw)
        xk   = 0.d0
        xeps = 0.d0
        if (     itytur.eq.2 .or. iturb.eq.50    &
            .or. iturb.eq.60) then
          xk = rtp(iel,ik)
        elseif (itytur.eq.3) then
          xk = (  rtp(iel,ir11) + rtp(iel,ir22)  &
                + rtp(iel,ir33)) / 2.d0
        endif
        if (     itytur.eq.2 .or. itytur.eq.3    &
            .or. iturb.eq.50) then
          xeps = rtp(iel,iep)
        elseif (iturb.eq.60) then
          xeps = cmu*rtp(iel,ik)*rtp(iel,iomg)
        endif
      else
        xabs = 0.d0
        xu   = 0.d0
        xv   = 0.d0
        xw   = 0.d0
        xk   = 0.d0
        xeps = 0.d0
      endif

      ! Broadcast to other ranks in parallel
      if (irangp.ge.0) then
        iun = 1
        call parbcr(irangv, iun, xabs)
        call parbcr(irangv, iun, xu)
        call parbcr(irangv, iun, xv)
        call parbcr(irangv, iun, xw)
        call parbcr(irangv, iun, xk)
        call parbcr(irangv, iun, xeps)
      endif

      if (irangp.le.0) write(impout,99) xabs, xu, xv, xw, xk, xeps

99    format(6g17.9)

    endif

  enddo

  if (irangp.le.0) close(impout)

endif

!===============================================================================
! 6. Example: print first calculated statistical moment
!===============================================================================

! ----------------------------------------------

! The test below allows checking that the following example compiles
! while disabling it by default.

iutile = 0

if (iutile.eq.0) return

! ----------------------------------------------

if (nbmomt.gt.0) then

  imom = 1 ! Moment number

  ! Position in 'propce' of the array of temporal accumulation for moments,
  ! propce(iel,ipcmom)

  ipcmom = ipproc(icmome(imom))

  ! The temporal accumulation for moments must be divided by the accumulated
  ! time, which id an array of size ncel or a single real number:
  ! - array of size ncel if idtmom(imom) > 0 : propce(iel, idtcm)
  ! - or simple real     if idtmom(imom) < 0 : dtcmom(idtcm)

  if (idtmom(imom).gt.0) then
    idtcm = ipproc(icdtmo(idtmom(imom)))
    do iel = 1, ncel
      write(nfecra, 4000)  &
           iel, propce(iel, ipcmom)/max(propce(iel, idtcm), epzero)
    enddo
  elseif (idtmom(imom).lt.0) then
    idtcm = -idtmom(imom)
    do iel = 1, ncel
      write(nfecra,4000)  &
           iel, propce(iel, ipcmom)/max(dtcmom(idtcm), epzero)
    enddo
  endif

endif

4000 format(' Cell ',i10,'   First moment ',e14.5)

!===============================================================================
! 6. Example: use of parallel utility functions for several operations
!===============================================================================

! This example demonstrates the parallel utility functions that may be used
! to simplify parallel operations.

! CAUTION: these routines modify their input

! ----------------------------------------------

! The test below allows checking that the following example compiles
! while disabling it by default.

iutile = 0

if (iutile.eq.0) return

! ----------------------------------------------

! Sum of an integer counter 'ii', here the number of cells

! local value
ii = ncel
! global sum
if (irangp.ge.0) then
  call parcpt(ii)
endif
! print the global sum
write(nfecra,5020)ii
 5020 format(' usproj: total number of cells = ', i10)

! Maximum of an integer counter 'ii', here the number of cells

! local value
ii = ncel
! global maximum
if (irangp.ge.0) then
  call parcmx(ii)
endif
! print the global maximum value
write(nfecra,5010)ii
 5010 format(' usproj: max. number of cells per process = ', i10)

! Sum of a real 'rrr', here the volume

! local value
rrr = 0.d0
do iel = 1, ncel
  rrr = rrr + volume(iel)
enddo
! global sum
if (irangp.ge.0) then
  call parsom(rrr)
endif
! print the global sum
write(nfecra,5030)rrr
 5030 format(' usproj: total domain volume = ', e14.5)

! Maximum of a real 'rrr', here the volume

! local value
rrr = 0.d0
do iel = 1, ncel
  if (volume(iel).gt.rrr) rrr = volume(iel)
enddo
! global maximum
if (irangp.ge.0) then
  call parmax(rrr)
endif
! print the global maximum
write(nfecra,5040)rrr
 5040 format(' usproj: max volume per process = ', e14.5)

! Minimum of a real 'rrr', here the volume

! local value
rrr = grand
do iel = 1, ncel
  if (volume(iel).lt.rrr) rrr = volume(iel)
enddo
! global minimum
if (irangp.ge.0) then
  call parmin(rrr)
endif
! print the global minimum
write(nfecra,5050)rrr
 5050 format(' usproj: min volume per process = ', e14.5)

! Maximum of a real and associated real values;
! here the volume and its location (3 coordinates)

nbr = 3
rrr  = -1.d0
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  if (rrr.lt.volume(iel)) then
    rrr = volume(iel)
    xyz(1) = xyzcen(1,iel)
    xyz(2) = xyzcen(2,iel)
    xyz(3) = xyzcen(3,iel)
  endif
enddo
! global maximum and associated location
if (irangp.ge.0) then
  call parmxl(nbr, rrr, xyz)
endif
! print the global maximum and its associated values
write(nfecra,5060) rrr, xyz(1), xyz(2), xyz(3)
 5060 format(' Usproj: Max. volume =      ', e14.5, /,  &
             '         Location (x,y,z) = ', 3e14.5)

! Minimum of a real and associated real values;
! here the volume and its location (3 coordinates)

nbr = 3
rrr  = 1.d+30
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  if (rrr.gt.volume(iel)) then
    rrr = volume(iel)
    xyz(1) = xyzcen(1,iel)
    xyz(2) = xyzcen(2,iel)
    xyz(3) = xyzcen(3,iel)
  endif
enddo
! global minimum and associated location
if (irangp.ge.0) then
  call parmnl(nbr,rrr,xyz)
endif
! print the global minimum and its associated values
write(nfecra,5070) rrr, xyz(1), xyz(2), xyz(3)
 5070 format(' Usproj: Min. volume =      ', e14.5, /,  &
             '         Location (x,y,z) = ', 3e14.5)

! Sum of an array of integers;
! here, the number of cells, faces, and boundary faces

! local values; note that to avoid counting interior faces on
! parallel boundaries twice, we check if 'ifacel(1,ifac) .le. ncel',
! as on a parallel boundary, this is always true for one domain
! and false for the other.

nbr = 3
itab(1) = ncel
itab(2) = 0
itab(3) = nfabor
do ifac = 1, nfac
  if (ifacel(1, ifac).le.ncel) itab(2) = itab(2) + 1
enddo
! global sum
if (irangp.ge.0) then
  call parism(nbr, itab)
endif
! print the global sums
write(nfecra,5080) itab(1), itab(2), itab(3)
 5080 format(' usproj: Number of cells =          ', i10, /,  &
             '         Number of interior faces = ', i10, /,  &
             '         Number of boundary faces = ', i10)

! Maxima from an array of integers;
! here, the number of cells, faces, and boundary faces

! local values
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
! global maxima
if (irangp.ge.0) then
  call parimx(nbr, itab)
endif
! print the global maxima
write(nfecra,5090) itab(1), itab(2), itab(3)
 5090 format(' usproj: Max. number of cells per proc. =          ', i10, /,  &
             '         Max. number of interior faces per proc. = ', i10, /,  &
             '         Max. number of boundary faces per proc. = ', i10)

! Minima from an array of integers;
! here, the number of cells, faces, and boundary faces

! local values
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
! global minima
if (irangp.ge.0) then
  call parimn(nbr, itab)
endif
! print the global minima
write(nfecra,5100) itab(1), itab(2), itab(3)
 5100 format(' usproj: Min. number of cells per proc. =          ', i10, /,  &
             '         Min. number of interior faces per proc. = ', i10, /,  &
             '         Min. number of boundary faces per proc. = ', i10)

! Sum of an array of reals;
! here, the 3 velocity components (so as to compute a mean for example)

! local values
nbr = 3
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  xyz(1) = xyz(1)+rtp(iel,iu)
  xyz(2) = xyz(2)+rtp(iel,iv)
  xyz(3) = xyz(3)+rtp(iel,iw)
enddo
! global sum
if (irangp.ge.0) then
  call parrsm(nbr, xyz)
endif
! print the global sums
write(nfecra,5110) xyz(1), xyz(2), xyz(3)
 5110 format(' usproj: Sum of U on the domain = ', e14.5, /,   &
             '         Sum of V on the domain = ', e14.5, /,   &
             '         Sum of V on the domain = ', e14.5)

! Maximum of an array of reals;
! here, the 3 velocity components

! local values
nbr = 3
xyz(1) = rtp(1,iu)
xyz(2) = rtp(1,iv)
xyz(3) = rtp(1,iw)
do iel = 1, ncel
  xyz(1) = max(xyz(1),rtp(iel,iu))
  xyz(2) = max(xyz(2),rtp(iel,iv))
  xyz(3) = max(xyz(3),rtp(iel,iw))
enddo
! global maximum
if (irangp.ge.0) then
  call parrmx(nbr, xyz)
endif
! print the global maxima
write(nfecra,5120) xyz(1), xyz(2), xyz(3)
 5120 format(' usproj: Maximum of U on the domain = ', e14.5, /,   &
             '         Maximum of V on the domain = ', e14.5, /,   &
             '         Maximum of V on the domain = ', e14.5)

! Maximum of an array of reals;
! here, the 3 velocity components

! local values
nbr = 3
xyz(1) = rtp(1,iu)
xyz(2) = rtp(1,iv)
xyz(3) = rtp(1,iw)
do iel = 1, ncel
  xyz(1) = min(xyz(1),rtp(iel,iu))
  xyz(2) = min(xyz(2),rtp(iel,iv))
  xyz(3) = min(xyz(3),rtp(iel,iw))
enddo
! global minimum
if (irangp.ge.0) then
  call parrmn(nbr, xyz)
endif
! print the global maxima
write(nfecra,5130) xyz(1), xyz(2), xyz(3)
 5130 format(' usproj: Minimum of U on the domain = ', e14.5, /,   &
             '         Minimum of V on the domain = ', e14.5, /,   &
             '         Minimum of V on the domain = ', e14.5)

! Broadcast an array of local integers to other ranks;
! in this example, we use the number of cells, interior faces, and boundary
! faces from process rank 0 (irangv).

! local values
irangv = 0
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
! broadcast from rank irangv to all others
if (irangp.ge.0) then
  call parbci(irangv, nbr, itab)
endif
! print values broadcast and received from rank 'irangv'
write(nfecra,5140) irangv, itab(1), itab(2), itab(3)
 5140 format(' usproj: On rank ', i10 , /,                     &
             '         Number of cells          = ', i10, /,   &
             '         Number of interior faces = ', i10, /,   &
             '         Number of boundary faces = ', i10)

! Broadcast an array of local reals to other ranks;
! in this example, we use 3 velocity values from process rank 0 (irangv).

! local values
irangv = 0
nbr = 3
xyz(1) = rtp(1,iu)
xyz(2) = rtp(1,iv)
xyz(3) = rtp(1,iw)
! broadcast from rank irangv to all others
if (irangp.ge.0) then
  call parbcr(irangv, nbr, xyz)
endif
! print values broadcast and received from rank 'irangv'
write(nfecra,5150) irangv, xyz(1), xyz(2), xyz(3)
 5150 format(' usproj: On rank ', i10 , /,                       &
             '         Velocity U in first cell = ', e14.5, /,   &
             '         Velocity V in first cell = ', e14.5, /,   &
             '         Velocity W in first cell = ', e14.5)

! All ranks obtain the global cell number from a given cell on a given rank;
! in this example, cell 'iel' prom rank 'irangv'.
iel = 1
irangv = 0
if (irangp.ge.0) then
  call parcel(iel, irangv, ielg)
else
  ielg = -1
endif
! print global cell number for cell 'iel' from rank 'irangv'
write(nfecra,5160) iel, irangv, ielg
 5160 format(' usproj: local cell iel =         ', i10, /,   &
             '         on rank irangv =         ', i10, /,   &
             '         has global number ielg = ', i10)

! Get the global number of a local cel 'iel';
! in this exemple, we use iel = 1, as all ranks should have at least one cell
iel = 1
call parclg(iel, irangp, ielg)
! each rank prints the global cell number of its cell number iel;
! (or 0 if it has less than iel cells)
write(nfecra,5170) iel, irangp, ielg
 5170 format(' usproj: Local cell number iel =  ', i10, /,    &
             '         on rank irangp =         ', i10, /,    &
             '         has global number ielg = ', i10)

! Get the global number of a local interior face 'ifac';
ifac = 1
call parfig(ifac, irangp, ifacg)
! each rank prints the global face number of its face number ifac;
! (or 0 if it has less than ifac faces)
write(nfecra,5180) ifac, irangp, ifacg
 5180 format(' usproj: Local face number ifac =  ', i10, /,  &
             '         on rank irangp =          ', i10, /,  &
             '         has global number ifacg = ', i10)

! Get the global number of a local boundary face 'ifac';
ifac = 1
call parfbg(ifac, irangp, ifacg)
! each rank prints the global face number of its boundary face number ifac;
! (or 0 if it has less than ifac boundary faces, which may occur)
write(nfecra,5190) ifac, irangp, ifacg
 5190 format(' usproj: Local boundary face number ifac =  ', i10, /,  &
             '         on rank irangp =                   ', i10, /,  &
             '         has global number ifacg =          ', i10)

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine

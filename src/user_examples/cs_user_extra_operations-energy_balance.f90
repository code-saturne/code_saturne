!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_user_extra_operations-energy_balance.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!> This is an example of cs_user_extra_operations.f90 which
!> performs an energy balance.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

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
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use field_operator
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel    , ifac   , ivar
integer          iel1   , iel2   , ieltsm
integer          iortho
integer          inc    , iccocg
integer          iflmas , iflmab , ipccp
integer          iscal
integer          ilelt  , nlelt

double precision xvara  , xvar
double precision xbilan , xbilvl , xbilpa , xbilpt
double precision xbilsy , xbilen , xbilso , xbildv
double precision xbilmi , xbilma
double precision xfluxf , xgamma
double precision flumab , ctb1, ctb2

integer, allocatable, dimension(:) :: lstelt

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: treco
double precision, allocatable, dimension(:) :: xcp
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cpro_visct, cpro_cp
double precision, dimension(:), pointer :: cvar_scal, cvara_scal
!< [loc_var_dec]

!===============================================================================
! Initialization
!===============================================================================

!< [init]
! Allocate a temporary array for cells or interior/boundary faces selection
allocate(lstelt(max(ncel,nfac,nfabor)))
allocate(xcp(ncel))
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

  call field_get_val_s(ivarfl(ivar), cvar_scal)
  call field_get_val_prev_s(ivarfl(ivar), cvara_scal)

  ! Physical quantity numbers
  call field_get_val_s(icrom, crom)
  call field_get_val_s(ivisct, cpro_visct)

  ! Pointers to the mass fluxes
  call field_get_key_int(ivarfl(ivar), kimasf, iflmas)
  call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  ! The balance is in Joule, so store the specific heat when dealing
  ! with the Temperature.

  if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

  ! If it is a temperature
  if (itherm.eq.1) then
    if (icp.ge.0) then
      do iel = 1, ncel
        xcp(iel) = cpro_cp(iel)
      enddo
    else
      do iel = 1, ncel
        xcp(iel) = cp0
      enddo
    endif
  else
    do iel = 1, ncel
      xcp(iel) = 1.d0
    enddo
  endif

  ! Boundary condition pointers for gradients and advection

  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)

  ! Boundary condition pointers for diffusion

  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  ! --> Compute value reconstructed at I' for boundary faces

  allocate(treco(nfabor))

  ! For non-orthogonal meshes, it must be equal to the value at the
  ! cell center, which is computed in:
  ! treco(ifac) (with ifac=1, nfabor)

  ! For orthogonal meshes, it is sufficient to assign:
  ! cvar_scal(iel) to treco(ifac), with iel=ifabor(ifac)
  ! (this option corresponds to the second branch of the test below,
  ! with iortho different from 0).

  iortho = 0

  ! --> General case (for non-orthogonal meshes)

  if (iortho.eq.0) then

    ! Allocate a work array for the gradient calculation
    allocate(grad(3,ncelet))

    ! --- Compute temperature gradient

    inc = 1
    iccocg = 1

    call field_gradient_scalar &
    !=========================
      (ivarfl(ivar), 0, imrgra, inc, iccocg,                               &
       grad)

    ! - Compute reconstructed value in boundary cells

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      treco(ifac) =   cvar_scal(iel)       &
                    + diipb(1,ifac)*grad(1,iel)  &
                    + diipb(2,ifac)*grad(2,iel)  &
                    + diipb(3,ifac)*grad(3,iel)
    enddo

    ! Free memory
    deallocate(grad)

  ! --> Case of orthogonal meshes

  else

    ! Compute reconstructed value
    ! (here, we assign the non-reconstructed value)

    do ifac = 1, nfabor
      iel = ifabor(ifac)
      treco(ifac) = cvar_scal(iel)
    enddo

  endif

  ! 2.1 Compute the balance at time step n
  ! ======================================

  ! --> Balance on interior volumes
  !     ---------------------------

  ! If it is variable, the density 'rom' has been computed at the beginning
  ! of the time step using the temperature from the previous time step.

  do iel = 1, ncel
    xvara = cvara_scal(iel)
    xvar  = cvar_scal(iel)
    xbilvl =   xbilvl                              &
             + volume(iel) * xcp(iel) * crom(iel)  &
                           * (xvara - xvar)
  enddo

  ! --> Balance on all faces (interior and boundary), for div(rho u)
  !     ------------------------------------------------------------

  do ifac = 1, nfac

    iel1 = ifacel(1,ifac)
    if (iel1.le.ncel) then
      ctb1 = imasfl(ifac)*xcp(iel1)*cvar_scal(iel1)*dt(iel1)
    else
      ctb1 = 0d0
    endif

    iel2 = ifacel(2,ifac)
    if (iel2.le.ncel) then
      ctb2 = imasfl(ifac)*xcp(iel2)*cvar_scal(iel2)*dt(iel2)
    else
      ctb2 = 0d0
    endif

    xbildv =  xbildv + (ctb1 - ctb2)
  enddo

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    xbildv = xbildv + dt(iel) * xcp(iel)             &
                              * bmasfl(ifac)         &
                              * cvar_scal(iel)
  enddo


  ! In case of a mass source term, add contribution from Gamma*Tn+1

  if (ncetsm.gt.0) then
    do ieltsm = 1, ncetsm
      iel = icetsm(ieltsm)
      xvar  = cvar_scal(iel)
      xgamma = smacel(ieltsm,ipr)
      xbildv =   xbildv                            &
               - volume(iel) * xcp(iel) * dt(iel)  &
                             * xgamma * xvar
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

    ! Physical variables

    flumab = bmasfl(ifac)

    ! Contribution to flux from the current face
    ! (diffusion and convection fluxes)

    xfluxf = - surfbn(ifac) * dt(iel)                       &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp(iel)                    &
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

    ! Physical variables

    flumab = bmasfl(ifac)

    ! Contribution to flux from the current face
    ! (diffusion and convection fluxes)

    xfluxf = - surfbn(ifac) * dt(iel)                       &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp(iel)                    &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilpt = xbilpt + xfluxf

  enddo

  ! Contribution from symmetries (should be 0).

  call getfbr('1', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Physical variables

    flumab = bmasfl(ifac)

    ! Contribution to flux from the current face
    ! (diffusion and convection fluxes)

    xfluxf = - surfbn(ifac) * dt(iel)                       &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp(iel)                    &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilsy = xbilsy + xfluxf

  enddo

  ! Contribution from inlet (color 3, diffusion and convection flux)

  call getfbr('3', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Physical variables

    flumab = bmasfl(ifac)

    ! Contribution to flux from the current face
    ! (diffusion and convection fluxes)

    xfluxf = - surfbn(ifac) * dt(iel)                       &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp(iel)                    &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilen = xbilen + xfluxf

  enddo

  ! Contribution from outlet (color 5, diffusion and convection flux)

  call getfbr('5', nlelt, lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    iel  = ifabor(ifac)   ! associated boundary cell

    ! Physical variables

    flumab = bmasfl(ifac)

    ! Contribution to flux from the current face
    ! (diffusion and convection fluxes)

    xfluxf = - surfbn(ifac) * dt(iel)                       &
             * (cofafp(ifac) + cofbfp(ifac)*treco(ifac))    &
           - flumab * dt(iel) * xcp(iel)                    &
             * (coefap(ifac) + coefbp(ifac)*treco(ifac))

    xbilso = xbilso + xfluxf

  enddo

  ! Now the work array for the temperature can be freed
  deallocate(treco)


  ! --> Balance on mass source terms
  !     ----------------------------

  ! We separate mass injections from suctions for better generality

  if (ncetsm.gt.0) then
    do ieltsm = 1, ncetsm
      ! depending on the type of injection we use the 'smacell' value
      ! or the ambient temperature
      iel = icetsm(ieltsm)
      xgamma = smacel(ieltsm,ipr)
      if (itypsm(ieltsm,ivar).eq.0 .or. xgamma.lt.0.d0) then
        xvar = cvar_scal(iel)
      else
        xvar = smacel(ieltsm,ivar)
      endif
      if (icp.ge.0) then
        if (xgamma.lt.0.d0) then
          xbilma =   xbilma  &
                   + volume(iel) * cpro_cp(iel) * dt(iel) * xgamma * xvar
        else
          xbilmi =   xbilmi  &
                   + volume(iel) * cpro_cp(iel) * dt(iel) * xgamma * xvar
        endif
      else
        if (xgamma.lt.0.d0) then
          xbilma =   xbilma  &
                   + volume(iel) * cp0 * dt(iel) * xgamma * xvar
        else
          xbilmi =   xbilmi  &
                   + volume(iel) * cp0 * dt(iel) * xgamma * xvar
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
deallocate(xcp)
!< [finalize]

return
end subroutine cs_f_user_extra_operations

!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_extra_operations-nusselt_calculation.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!> See \subpage cs_user_extra_operations_examples and
!> \subpage cs_user_extra_operations-nusselt_calculation for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
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

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)
!< [loc_var_f_user]
! Local variables

integer iscal, ivar, iortho, iprev
integer inc, iccocg, ilelt, irangv
integer ifac, iel, npoint, iel1, irang1, ii, iun, impout, nlelt, neltg

double precision diipbx, diipby, diipbz
double precision xtbulk, xubulk, xyz(3), xtb, xub, tfac, lambda, xab

character*19 namfil

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: coefap, coefbp, cofafp
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: treco,treglo,treloc
double precision, allocatable, dimension(:) :: xnusselt
double precision, allocatable, dimension(:) :: xabs, xabsg
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar, viscl

double precision :: height, prandtl, qwall
parameter (height = 1.0d0)
parameter (prandtl = 1.0d0)
parameter (qwall = 1.0d0)
!< [loc_var_f_user]
!***********************************************************************


! Calculation of the Nusselt number
!< [nusselt_number]
if (ntcabs.eq.ntmabs) then

  allocate(lstelt(max(ncel,nfac,nfabor)))
  call field_get_val_v(ivarfl(iu), vel)

  iscal = iscalt         ! temperature scalar number
  ivar =  isca(iscal)    ! temperature variable number
  call field_get_val_s(iviscl, viscl)
  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)

  call field_get_val_s(ivarfl(ivar), cvar)
!< [nusselt_number]
  ! --> Compute value reconstructed at I' for boundary faces
!< [compute_nusselt]
  allocate(treco(nfabor))

  iortho = 0

!< [compute_nusselt]
  ! --> General case (for non-orthogonal meshes)
!< [gen_nusselt]
  if (iortho.eq.0) then
!< [gen_nusselt]
    ! Allocate a work array for the gradient calculation
!< [allocate_nusselt]
    allocate(grad(3,ncelet))
!< [allocate_nusselt]
    ! - Compute gradient
!< [gradient_nusselt]
    iprev = 0
    inc = 1
    iccocg = 1

    call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,    &
                               iccocg,                              &
                               grad)
!< [gradient_nusselt]
    ! - Compute reconstructed value in boundary cells
!< [value_nusselt]
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      treco(ifac) =   coefap(ifac) + coefbp(ifac)*(cvar(iel)       &
           + diipbx*grad(1,iel)  &
           + diipby*grad(2,iel)  &
           + diipbz*grad(3,iel))
    enddo
!< [value_nusselt]
    ! Free memory
!< [free_nusselt]
    deallocate(grad)
!< [free_nusselt]
    ! --> Case of orthogonal meshes
!< [else_nusselt]
  else
!< [else_nusselt]
    ! Compute reconstructed value
    ! (here, we assign the non-reconstructed value)
!< [value_ortho_nusselt]
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      treco(ifac) = coefap(ifac) + coefbp(ifac)*cvar(iel)
    enddo

  endif

  impout = impusr(1)
  if (irangp.le.0) then
    open(file="Nusselt.dat",unit=impout)
  endif

  call getfbr('normal[0,-1,0,0.1] and y < 0.01',nlelt,lstelt)

  neltg = nlelt
  if (irangp.ge.0) then
    call parcpt(neltg)
  endif

  allocate(xabs(nlelt))
  allocate(treloc(nlelt))
  allocate(xnusselt(neltg))
  if (irangp.ge.0) then
    allocate(xabsg(neltg))
    allocate(treglo(neltg))
  endif


  do ilelt = 1, nlelt
    ifac = lstelt(ilelt)
    xabs(ilelt) = cdgfbo(1,ifac)
    treloc(ilelt) = treco(ifac)
  enddo

  if (irangp.ge.0) then
    call paragv(nlelt, neltg, xabs, xabsg)
    call paragv(nlelt, neltg, treloc, treglo)
  endif

  do ilelt = 1,neltg
!< [value_ortho_nusselt]
    ! Calculation of the bulk temperature
!< [bulk_nusselt]
    if (irangp.ge.0) then
      xab = xabsg(ilelt)
    else
      xab = xabs(ilelt)
    endif

    xtbulk = 0.0d0
    xubulk = 0.0d0

    npoint = 200
    iel1   = -999
    do ii = 1, npoint
      xyz(1) = xab
      xyz(2) = float(ii-1)/float(npoint-1)
      xyz(3) = 0.d0
      call findpt(ncelet, ncel, xyzcen, xyz(1), xyz(2), xyz(3), iel, irangv)
      !==========
      if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
        iel1   = iel
        irang1 = irangv
        if (irangp.eq.irangv) then
          xtb = volume(iel)*cvar(iel)*vel(1,iel)
          xub = volume(iel)*vel(1,iel)
          xtbulk = xtbulk + xtb
          xubulk = xubulk + xub
          lambda = cp0 * viscl(iel) / prandtl
        endif
      endif
    enddo

    if (irangp.ge.0) then
      iun = 1
      call parall_bcast_r(irangv, lambda)
      call parsom(xtbulk)
      call parsom(xubulk)
    endif

    xtbulk = xtbulk/xubulk
    if (irangp.ge.0) then
      tfac = treglo(ilelt)
    else
      tfac = treloc(ilelt)
    endif

    xnusselt(ilelt) = qwall * 2.d0 * height / lambda / (tfac - xtbulk)

  enddo

  if (irangp.eq.-1) then
    do ii = 1, neltg
      write(impout,'(2E17.9)') xabs(ii)*10, xnusselt(ii)/(0.023d0*30000.d0**(0.8d0)*0.71d0**(0.4d0))
    enddo
  endif

  if (irangp.eq.0) then
    call sortc2(xabsg, xnusselt, neltg)
    do ii = 1, neltg
      write(impout,'(2E17.9)') xabsg(ii)*10, xnusselt(ii)/(0.023d0*30000.d0**(0.8d0)*0.71d0**(0.4d0))
    enddo
  endif

  close(impout)

  deallocate(treco)
  deallocate(xabs)
  deallocate(xnusselt)
  deallocate(lstelt)
  deallocate(treloc)
  if (irangp.ge.0) then
    deallocate(xabsg)
    deallocate(treglo)
  endif
endif

return
end subroutine cs_f_user_extra_operations
!< [bulk_nusselt]
subroutine sortc2(tab1, tab2, n)
!================
!< [loc_var_sortc2]
implicit none
integer n
double precision tab1(n), tab2(n)

integer ns, ii, jj, kk
double precision tabmin, t2
!< [loc_var_sortc2]
!< [body_sortc2]
ns = 1
do ii = 1, n-1
  tabmin = 1.d20
  kk = 0
  do jj = ns,n
    if (tabmin.gt.tab1(jj)) then
      tabmin = tab1(jj)
      kk = jj
    endif
  enddo
  t2 = tab2(kk)
  tab1(kk) = tab1(ns)
  tab2(kk) = tab2(ns)
  tab2(ns) = t2
  tab1(ns) = tabmin
  ns = ns + 1
enddo

return
end subroutine sortc2
!< [body_sortc2]

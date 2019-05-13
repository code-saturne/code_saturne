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

!> \file cs_user_initialization-pulverized_coal.f90
!>
!> \brief Pulverized coal example
!>
!> See \subpage cs_user_initialization for examples.
!>
!
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


subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, ige, mode, icla, icha
integer          ioxy

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision f1mc(ncharm), f2mc(ncharm)
double precision xkent, xeent, d2s3
double precision wmh2o,wmco2,wmn2,wmo2,dmas

integer, allocatable, dimension(:) :: lstelt

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi, cvar_fb
double precision, dimension(:), pointer :: cvar_omg, cvar_nusa
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvar_xch, cvar_xck, cvar_np
double precision, dimension(:), pointer :: cvar_h2, cvar_scalt
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m, cvar_f3m
double precision, dimension(:), pointer :: cvar_f6m, cvar_f7m
double precision, dimension(:), pointer :: cvar_yco2, cvar_yhcn, cvar_yno
double precision, dimension(:), pointer :: cvar_taire
!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

! Control Print

write(nfecra,9001)

! Local variables initialization

d2s3 = 2.d0/3.d0

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if (isuite.eq.0) then

! --> Initialisation of k and epsilon (exemple)

  xkent = 1.d-10
  xeent = 1.d-10

! ---- TURBULENCE

  if (itytur.eq.2) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)

    do iel = 1, ncel
      cvar_k(iel)  = xkent
      cvar_ep(iel) = xeent
    enddo

  elseif (itytur.eq.3) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(ir12), cvar_r12)
    call field_get_val_s(ivarfl(ir13), cvar_r13)
    call field_get_val_s(ivarfl(ir23), cvar_r23)
    call field_get_val_s(ivarfl(iep), cvar_ep)

    do iel = 1, ncel
      cvar_r11(iel) = d2s3*xkent
      cvar_r22(iel) = d2s3*xkent
      cvar_r33(iel) = d2s3*xkent
      cvar_r12(iel) = 0.d0
      cvar_r13(iel) = 0.d0
      cvar_r23(iel) = 0.d0
      cvar_ep(iel)  = xeent
    enddo

  elseif (iturb.eq.50) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)
    call field_get_val_s(ivarfl(iphi), cvar_phi)
    call field_get_val_s(ivarfl(ifb), cvar_fb)

    do iel = 1, ncel
      cvar_k(iel)   = xkent
      cvar_ep(iel)  = xeent
      cvar_phi(iel) = d2s3
      cvar_fb(iel)  = 0.d0
    enddo

  elseif (iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)

    do iel = 1, ncel
      cvar_k(iel)   = xkent
      cvar_omg(iel) = xeent/cmu/xkent
    enddo

  elseif (iturb.eq.70) then
    call field_get_val_s(ivarfl(inusa), cvar_nusa)

    do iel = 1, ncel
      cvar_nusa(iel) = cmu*xkent**2/xeent
    enddo

  endif

! --> All the domain is filled with the first oxidizer at TINITK
!                   ============================================

! ---- Computation of H1INIT and H2INIT

  t1init = t0
  t2init = t0

! ------ Transported variables for the solid phase
!         initialy lacking

  do icla = 1, nclacp
    icha = ichcor(icla)

    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xch)
    call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xck)
    call field_get_val_s(ivarfl(isca(inp(icla))), cvar_np)
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2)

    do iel = 1, ncel

      cvar_xch(iel) = zero
      cvar_xck(iel) = zero
      cvar_np(iel) = zero
      cvar_h2(iel)  = zero
    enddo
  enddo

! ------ Transported variables for the mix (solid+carrying gas)^2

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo

!       Oxidizer are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
!       the composition of the fisrt oxidiser is taken in account

  coefe(io2) = wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ih2o) = wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ico2) = wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
  mode = -1
  call cs_coal_htconvers1(mode, h1init, coefe, f1mc, f2mc, t1init)
 !============================

  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

  do iel = 1, ncel
    cvar_scalt(iel) = h1init
  enddo

! ------ Transported variables for the mix (passive scalars, variance)

  do icha = 1, ncharb
    call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
    call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
    do iel = 1, ncel
      cvar_f1m(iel) = zero
      cvar_f2m(iel) = zero
    enddo
  enddo

  call field_get_val_s(ivarfl(isca(if3m)), cvar_f3m)
  call field_get_val_s(ivarfl(isca(if6m)), cvar_f6m)
  call field_get_val_s(ivarfl(isca(if7m)), cvar_f7m)
  call field_get_val_s(ivarfl(isca(iyco2)), cvar_yco2)
  call field_get_val_s(ivarfl(isca(iyhcn)), cvar_yhcn)
  call field_get_val_s(ivarfl(isca(iyno)), cvar_yno)
  call field_get_val_s(ivarfl(isca(itaire)), cvar_taire)

  do iel = 1, ncel

    cvar_f3m(iel) = zero

    if (noxyd .ge. 2) then
      cvar_f6m(iel) = zero
    endif

    if (noxyd .eq. 3) then
      cvar_f7m(iel) = zero
    endif

    if (ieqco2.ge.1) then

      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)

      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      xco2 = oxyco2(ioxy)*wmco2/dmas

      cvar_yco2(iel) = oxyco2(ioxy)*wmco2/dmas

    endif

    if (ieqnox .eq. 1) then
      cvar_yhcn(iel) = 0.d0
      cvar_yno(iel) = 0.d0
      cvar_taire(iel) = 293.d0
    endif

  enddo

endif
!< [init]

!--------
! Formats
!--------

 9001 format(                                                   /,&
'  cs_user_initialization : Variables Initialisation for'      ,/,&
'                    Pulverized Coal by User'                  ,/,&
                                                                /)

!----
! End
!----

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_f_initialization

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

!> \file cs_user_initialization-fuel.f90
!>
!> \brief Fuel example
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
integer          iel, ige, mode, icla
integer          ioxy

double precision t1init, h1init, coefe(ngazem)
double precision t2init, h2init
double precision xkent, xeent, d2s3
double precision dmas , wmco2 , wmh2o , wmn2 , wmo2

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_fb, cvar_omg
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvar_yfol, cvar_ng, cvar_h2
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_fvap, cvar_f7m, cvar_fvp2m
double precision, dimension(:), pointer :: cvar_yco2, cvar_yhcn, cvar_yno
double precision, dimension(:), pointer :: cvar_hox

!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

write(nfecra,9001)


!< [init]
d2s3 = 2.d0/3.d0

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if ( isuite.eq.0 ) then

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

  endif

! --> All the computation domain is initialized with air at TINITK
!             ====================================================

! ---- Computation of H1INIT and  H2INIT

  t1init = 1000.d0
  t2init = 1000.d0

! ------ Transported variables for droplets

  h2init = h02fol +  cp2fol*(t2init-trefth)

  do icla = 1, nclafu
    call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfol)
    call field_get_val_s(ivarfl(isca(ing(icla))), cvar_ng)
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2)
    do iel = 1, ncel
      cvar_yfol(iel) = zero
      cvar_ng(iel)  = zero
      cvar_h2(iel)  = zero
    enddo
  enddo

! ------ Transported variables for the mix (droplets and carrying gases)

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo
!  On considere l'oxydant 1
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

  mode = -1
  call cs_fuel_htconvers1(mode,h1init,coefe,t1init)
 !============================

  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
  do iel = 1, ncel
    cvar_scalt(iel) = h1init
  enddo

! ------ Transported variables for gaseous mixture
!        (passive scalars, variance, reactive species)

  call field_get_val_s(ivarfl(isca(ifvap)), cvar_fvap)
  call field_get_val_s(ivarfl(isca(if7m)), cvar_f7m)
  call field_get_val_s(ivarfl(isca(ifvp2m)), cvar_fvp2m)
  call field_get_val_s(ivarfl(isca(iyco2)), cvar_yco2)
  call field_get_val_s(ivarfl(isca(iyhcn)), cvar_yhcn)
  call field_get_val_s(ivarfl(isca(iyno)), cvar_yno)
  call field_get_val_s(ivarfl(isca(ihox)), cvar_hox)

  do iel = 1, ncel
    cvar_fvap(iel) = 0.d0
    cvar_f7m(iel) = zero
    cvar_fvp2m(iel) = zero
    if ( ieqco2 .ge. 1 ) then
      ioxy   =  1
      wmo2   = wmole(io2)
      wmco2  = wmole(ico2)
      wmh2o  = wmole(ih2o)
      wmn2   = wmole(in2)
      dmas = ( oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2               &
              +oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2 )
      cvar_yco2(iel) = oxyco2(ioxy)*wmco2/dmas
    endif
    if ( ieqnox .eq. 1 ) then
      cvar_yhcn(iel) = zero
      cvar_yno(iel) = zero
      cvar_hox(iel) = h1init
    endif
  enddo

endif
!< [init]


!--------
! Formats
!--------

 9001 format(                                                   /,&
'  Variables Initialisation for Fuel by the User'              ,/,&
                                                                /)

!----
! End
!----

return
end subroutine cs_user_f_initialization

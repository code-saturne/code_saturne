!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine diffst &
!================

 ( nscal, iterns )

!===============================================================================
! Function :
! --------

! Weakly compressible algorithm (semi-analytic):
!  Computation of scalar diffusion terms

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! iterns           ! i  ! <-- ! Navier-Stokes iteration number                 !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nscal , iterns

! Local variables

integer          ivar  , iel   , ifac  , iscal, f_id0
integer          nswrgp, imligp, iwarnp
integer          iccocg, inc, imrgrp
integer          iconvp, idiffp, ircflp
integer          ischcp, isstpp
integer          ifcvsl, iflmas, iflmab
integer          imucpp, idftnp, imasac
integer          key_buoyant_id, is_buoyant_fld
double precision epsrgp, climgp, extrap
double precision blencp, relaxp, thetex
double precision turb_schmidt

integer          icvflb
integer          ivoid(1)

double precision rvoid(1)

double precision, allocatable, dimension(:) :: vistot, viscf, viscb
double precision, allocatable, dimension(:) :: xcpp
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: visct, cpro_cp, cpro_viscls
double precision, dimension(:), pointer :: cvar_scal
double precision, dimension(:), pointer :: cpro_tsscal

type(var_cal_opt) :: vcopt

!===============================================================================

! Memory allocation
allocate(vistot(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(xcpp(ncelet))

if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

do iscal = 1, nscal

  ! Index for variable
  ivar = isca(iscal)

  call field_get_val_s(ivarfl(isca(iscal)), cvar_scal)

  ! Key id for buoyant field (inside the Navier Stokes loop)
  call field_get_key_id("is_buoyant", key_buoyant_id)
  call field_get_key_int(ivarfl(isca(iscal)), key_buoyant_id, is_buoyant_fld)

  ! If the scalar is buoyant, it is inside the Navier Stokes loop, and so iterns >=1
  ! otherwise it is outside of the loop and iterns = -1.
  if (  (is_buoyant_fld.eq. 1 .and. iterns.eq.-1) &
    .or.(is_buoyant_fld.eq. 0 .and. iterns.ne.-1)) cycle

  imucpp = 0
  if (iscavr(iscal).gt.0) then
    if (abs(iscacp(iscavr(iscal))).eq.1) then
      imucpp = 1
    endif
  else
    if (abs(iscacp(iscal)).eq.1) then
      imucpp = 1
    endif
  endif

  if (imucpp.eq.0) then
    do iel = 1, ncel
      xcpp(iel) = 1.d0
    enddo
  elseif (imucpp.eq.1) then
    if (icp.ge.0) then
      do iel = 1, ncel
        xcpp(iel) = cpro_cp(iel)
      enddo
   else
      do iel = 1, ncel
        xcpp(iel) = cp0
      enddo
    endif
  endif

  ! Handle parallelism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(xcpp)
  endif

  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

  f_id0  = -1
  iconvp = 0 ! diffusion term only
  ischcp = 1
  isstpp = 1
  blencp = 0.d0
  imasac = 0
  idiffp = 1
  imrgrp = vcopt%imrgra
  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  ircflp = vcopt%ircflu
  inc    = 1
  iccocg = 1
  idftnp = ISOTROPIC_DIFFUSION !FIXME when activating GGDH
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag
  relaxp = vcopt%relaxv
  thetex = 1.d0
  ! all boundary convective flux with upwind
  icvflb = 0

  ! Pointers to the mass fluxes
  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  ! Diffusion velocity

  ! Index for molecular diffusivity
  call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_viscls)
  endif

  ! Index for turbulent diffusivity
  call field_get_val_s(ivisct, visct)

  if (vcopt%idiff.ge.1) then

    ! Only the positive part of mu_t is considered (MAX(mu_t,0)),
    ! Dynamic LES can cause negative mu_t (clipping on (mu+mu_t))
    ! The positive part of (K+K_t) would have been considered
    ! but should allow negative K_t that is considered non physical here

    call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

    if(ifcvsl.lt.0)then
      do iel = 1, ncel
        vistot(iel) = visls0(iscal)                                     &
           + vcopt%idifft*xcpp(iel)*max(visct(iel),zero)/turb_schmidt
      enddo
    else
      do iel = 1, ncel
        vistot(iel) = cpro_viscls(iel)                                  &
           + vcopt%idifft*xcpp(iel)*max(visct(iel),zero)/turb_schmidt
      enddo
    endif

    call viscfa ( imvisf , vistot , viscf , viscb )
    !==========

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo
    do iel = 1, ncel
      vistot(iel) = 0.d0
    enddo

  endif

  ! Source term
  call field_get_val_s(iustdy(iscal), cpro_tsscal)

  ! Diffusion term calculation
  call field_get_coefa_s(ivarfl(ivar), coefap)
  call field_get_coefb_s(ivarfl(ivar), coefbp)
  call field_get_coefaf_s(ivarfl(ivar), cofafp)
  call field_get_coefbf_s(ivarfl(ivar), cofbfp)

  call bilsca &
  !==========
  ( idtvar , f_id0  , iconvp , idiffp , nswrgp , imligp , ircflp , &
    ischcp , isstpp , inc    , imrgrp , iccocg ,                   &
    iwarnp , imucpp , idftnp , imasac ,                            &
    blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
    cvar_scal       , cvar_scal       ,                            &
    coefap , coefbp , cofafp , cofbfp ,                            &
    imasfl , bmasfl ,                                              &
    viscf  , viscb  , rvoid  , xcpp   ,                            &
    rvoid  , rvoid  ,                                              &
    icvflb , ivoid  ,                                              &
    cpro_tsscal     )

enddo

! Free memory
deallocate(viscf, viscb)
deallocate(vistot)
deallocate(xcpp)

!----
! Formats
!----

!----
! End
!----

return
end subroutine

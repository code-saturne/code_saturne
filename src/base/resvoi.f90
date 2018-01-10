!-------------------------------------------------------------------------------

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
! Function:
! ---------

!> \file resvoi.f90
!>
!> \brief Solving the void fraction \f$ \alpha \f$ for the Volume of Fluid
!>        method (and hence for cavitating flows).
!>
!> This function solves the equation:
!> \f[
!> \dfrac{\alpha^n - \alpha^{n-1}}{\Delta t}
!>     + \divs \left( \alpha^n \vect{u}^n \right)
!>     = \dfrac{\Gamma_V \left( \alpha^{n-1}, p^n \right)}{\rho_v}
!> \f]
!> with \f$ \Gamma_V \f$ the eventual vaporization source term (Merkle model) in
!> case the cavitation model is enabled and \f$ \rho_v \f$ the reference gas
!> density.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     dt            time step (per cell)
!> \param[in]     iterns        Navier-Stokes iteration number
!_______________________________________________________________________________

subroutine resvoi &
 ( dt     , iterns )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use dimens
use paramx
use numvar
use entsor
use optcal
use pointe, only: gamcav, dgdpca
use mesh
use field
use cavitation
use vof
use parall
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision dt(ncelet)
integer          iterns

! Local variables

integer          ivar  , iel, ifac, f_id
integer          init
integer          nswrgp, imligp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          iflmas, iflmab
integer          iwarnp
integer          imucpp, idftnp, iswdyp

integer          icvflb
integer          ivoid(1)
integer          kscmin, kscmax, iclmin(1), iclmax(1)

double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp, epsrsp, thetap

double precision rvoid(1)
double precision vmin(1), vmax(1)
double precision dtmaxl, dtmaxg
double precision scmaxp, scminp
double precision thets, thetv, tsexp

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt
double precision, allocatable, dimension(:) :: dpvar, divu
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: c_st_voidf
double precision, dimension(:), pointer :: cvar_pr, cvara_pr
double precision, dimension(:), pointer :: cvar_voidf, cvara_voidf

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

ivar = ivolf2

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

call field_get_val_s(ivarfl(ivolf2), cvar_voidf)
call field_get_val_prev_s(ivarfl(ivolf2), cvara_voidf)

! implicitation in pressure of the vaporization/condensation model (cavitation)
if (icavit.ge.0.and.itscvi.eq.1) then
  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_prev_s(ivarfl(ipr), cvara_pr)
endif

! Allocate temporary arrays
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrs(ncelet),rovsdt(ncelet))

! Allocate work arrays
allocate(dpvar(ncelet))
allocate(divu(ncelet))

! --- Boundary conditions

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! --- Physical quantities

call field_get_key_int(ivarfl(ivar), kimasf, iflmas)
call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Key id for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

! Theta related to explicit source terms
thets  = thetsn
if (isno2t.gt.0) then
  call field_get_key_int(ivarfl(ivolf2), kstprv, f_id)
  call field_get_val_s(f_id, c_st_voidf)
else
  c_st_voidf => null()
endif

! Theta related to void fraction
thetv = vcopt%thetav

! --- Initialization

do iel = 1, ncel
   smbrs(iel) = 0.d0
enddo

do iel = 1, ncel
   rovsdt(iel) = 0.d0
enddo

! Arbitrary initialization (no diffusion for void fraction)
do ifac = 1, nfac
   viscf(ifac) = 1.d0
enddo
do ifac = 1, nfabor
   viscb(ifac) = 1.d0
enddo

!===============================================================================
! 2. Preliminary computations
!===============================================================================

! Update the cavitation source term with pressure increment
!   if it has been implicited in pressure at correction step,
!   in order to ensure mass conservation.

if (icavit.ge.0.and.itscvi.eq.1) then
  do iel = 1, ncel
    gamcav(iel) = gamcav(iel) + dgdpca(iel)*(cvar_pr(iel)-cvara_pr(iel))
  enddo
endif

! Compute the limiting time step to satisfy min/max principle.
!   Only if a source term is accounted for.

dtmaxl = 1.d15
dtmaxg = 1.d15

if (icavit.gt.0) then
  do iel = 1, ncel
    if (gamcav(iel).lt.0.d0) then
      dtmaxl = -rho2*cvara_voidf(iel)/gamcav(iel)
    else
      dtmaxl = rho1*(1.d0-cvara_voidf(iel))/gamcav(iel)
    endif
    dtmaxg = min(dtmaxl,dtmaxg)
  enddo
  if (irangp.ge.0) call parmin(dtmaxg)

  if (dt(1).gt.dtmaxg)  write(nfecra,1000) dt(1), dtmaxg
endif

!===============================================================================
! 3. Construct the system to solve
!===============================================================================

! Source terms
!-------------

! Cavitation source term (explicit)
if (icavit.ge.0) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + cell_f_vol(iel)*gamcav(iel)/rho2
  enddo
endif

! Source term linked with the non-conservative form of convection term
! in codits (always implicited)

init = 1
call divmas (init,imasfl,bmasfl,divu)

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) - divu(iel)
enddo

! Source terms assembly for codits

! If source terms are extrapolated over time
if (isno2t.gt.0) then
  do iel = 1, ncel
    tsexp = c_st_voidf(iel)
    c_st_voidf(iel) = smbrs(iel)
    smbrs(iel) = -thets*tsexp + (1.d0+thets)*c_st_voidf(iel) &
                 + rovsdt(iel)*cvara_voidf(iel)
    rovsdt(iel) = -thetv*rovsdt(iel)
  enddo
! If source terms are not extrapolated over time
else
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + rovsdt(iel)*cvara_voidf(iel)
    rovsdt(iel) = -rovsdt(iel)
  enddo
endif

! Unteady term
!-------------

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) + vcopt%istat*cell_f_vol(iel)/dt(iel)
enddo

!===============================================================================
! 3. Solving
!===============================================================================

! Solving void fraction
iconvp = vcopt%iconv
idiffp = vcopt%idiff
ndircp = ndircl(ivar)
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
imucpp = 0
idftnp = vcopt%idften
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
relaxp = vcopt%relaxv
thetap = vcopt%thetav
! all boundary convective flux with upwind
icvflb = 0

call codits &
!==========
 ( idtvar , iterns , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_voidf     , cvara_voidf     ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbrs  , cvar_voidf      , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 3. Clipping: only if min/max principle is not satisfied for cavitation
!===============================================================================

iclmax(1) = 0
iclmin(1) = 0

if ((icavit.gt.0.and.dt(1).gt.dtmaxg).or.icavit.le.0) then

  ! --- Calcul du min et max
  vmin(1) = cvar_voidf(1)
  vmax(1) = cvar_voidf(1)
  do iel = 1, ncel
    vmin(1) = min(vmin(1),cvar_voidf(iel))
    vmax(1) = max(vmax(1),cvar_voidf(iel))
  enddo

  ! Get the min and max clipping
  call field_get_key_double(ivarfl(ivar), kscmin, scminp)
  call field_get_key_double(ivarfl(ivar), kscmax, scmaxp)

  if(scmaxp.gt.scminp) then
    do iel = 1, ncel
      if(cvar_voidf(iel).gt.scmaxp)then
        iclmax(1) = iclmax(1) + 1
        cvar_voidf(iel) = scmaxp
      endif
      if(cvar_voidf(iel).lt.scminp)then
        iclmin(1) = iclmin(1) + 1
        cvar_voidf(iel) = scminp
      endif
    enddo
  endif

endif

call log_iteration_clipping_field(ivarfl(ivar), iclmin(1), iclmax(1), &
                                  vmin, vmax, iclmin(1), iclmax(1))

! Free memory
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)
deallocate(dpvar)
deallocate(divu)

!--------
! Formats
!--------

 1000   format(                                                 /,&
'@',                                                            /,&
'@ @@ WARNING: Void fraction resolution',                       /,&
'@    ========',                                                /,&
'@  The current time step is too large to ensure the min/max',  /,&
'@     principle on void fraction.',                            /,&
'@'                                                             /,&
'@  The current time step is', E13.5,' while',                  /,&
'@     the maximum admissible value is', E13.5,                 /,&
'@'                                                             /,&
'@  Clipping on void fraction should occur and',                /,&
'@     mass conservation is lost.',                             /,&
'@ ',                                                           /)

!----
! End
!----

return

end subroutine

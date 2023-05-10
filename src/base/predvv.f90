!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file predvv.f90
!>
!> \brief This subroutine performs the velocity prediction step of the Navier
!> Stokes equations for incompressible or slightly compressible flows for
!> the coupled velocity components solver.
!>
!> - at the first call, the predicted velocities are computed and also
!>   an estimator on the predicted velocity is computed.
!>
!> - at the second call, a global estimator on Navier Stokes is computed.
!>   This second call is done after the correction step
!>   (\ref cs_pressure_correction).
!>
!> Please refer to the
!> <a href="../../theory.pdf#predvv"><b>predvv</b></b></a> section
!> of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iappel        call number (1 or 2)
!> \param[in]     iterns        index of the iteration on Navier-Stokes
!> \param[in]     dt            time step (per cell)
!> \param[in]     vel           velocity
!> \param[in]     vela          velocity at the previous time step
!> \param[in]     velk          velocity at the previous sub iteration (or vela)
!> \param[in,out] da_uu         velocity matrix
!> \param[in]     coefav        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (implicit part)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     trava         working array for the velocity-pressure coupling
!> \param[out]    dfrcxt        variation of the external forces
!                               making the hydrostatic pressure
!> \param[in]     grdphd        hydrostatic pressure gradient to handle the
!>                              imbalance between the pressure gradient and
!>                              gravity source term
!> \param[in]     tpucou        non scalar time step in case of
!>                              velocity pressure coupling
!> \param[in]     trav          right hand side for the normalizing
!>                              the residual
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     viscfi        same as viscf for increments
!> \param[in]     viscbi        same as viscb for increments
!> \param[in]     secvif        secondary viscosity at interior faces
!> \param[in]     secvib        secondary viscosity at boundary faces
!_______________________________________________________________________________

subroutine predvv &
 ( iappel , iterns ,                                              &
   dt     , vel    , vela   , velk   , da_uu  ,                   &
   coefav , coefbv , cofafv , cofbfv ,                            &
   frcxt  , grdphd ,                                              &
   trava  ,                   dfrcxt , tpucou , trav   ,          &
   viscf  , viscb  , viscfi , viscbi , secvif , secvib )          &
  bind(C, name='cs_velocity_prediction')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use cstphy
use cstnum
use optcal
use parall
use period
use lagran
use ppppar
use ppthch
use ppincl
use cplsat
use mesh
use rotation
use turbomachinery
use cs_f_interfaces
use cs_c_bindings
use cfpoin
use field
use field_operator
use cavitation
use vof
use atincl, only: kopint, iatmst, ps
use pointe

!===============================================================================

implicit none

! Arguments

integer, value :: iappel, iterns

double precision dt(ncelet)
double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
double precision grdphd(3, ncelet)
double precision trava(ndim,ncelet)
double precision tpucou(6, ncelet)
double precision trav(3,ncelet)
double precision viscf(*), viscb(nfabor)
double precision viscfi(*), viscbi(nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision coefav(3  ,nfabor)
double precision cofafv(3  ,nfabor)
double precision coefbv(3, 3, nfabor)
double precision cofbfv(3, 3, nfabor)

double precision vel   (3, ncelet)
double precision velk  (3, ncelet)
double precision vela  (3, ncelet)
double precision da_uu (6, ncelet)

! Local variables

integer          f_id  , iel   , ielpdc, ifac  , isou  , itypfl, n_fans
integer          inc   , iprev , init  , ii    , jj
integer          imrgrp, nswrgp, imligp, iwarnp
integer          idftnp
integer          nswrsp, imvisp
integer          iescap, iespre, iescor, iestot
integer          iflmb0
integer          idtva0, icvflb, f_oi_id
integer          jsou  , ivisep, imasac
integer          f_dim , iflwgr
integer          iflmas, iflmab
integer          ivoid(1)

double precision rnorm , vitnor
double precision romvom, drom  , rom
double precision epsrgp, climgp
double precision vit1  , vit2  , vit3, xkb, pip, pfac
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision d2s3  , thetap, thetp1, thets
double precision diipbx, diipby, diipbz
double precision dvol
double precision n(3)
double precision tensor(6)
double precision rscp, tref

double precision rvoid(1)

! Working arrays
double precision, allocatable, dimension(:,:) :: eswork
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:), target :: grad
double precision, allocatable, dimension(:,:), target :: hl_exp
double precision, dimension(:,:), allocatable :: smbr
double precision, dimension(:,:,:), allocatable :: fimp
double precision, dimension(:,:,:), allocatable :: fimpcp
double precision, dimension(:,:), allocatable :: gavinj
double precision, allocatable, dimension(:,:), target :: loctsexp
double precision, allocatable, dimension(:,:,:), target :: loctsimp
double precision, dimension(:,:), pointer :: tsexp
double precision, dimension(:,:,:), pointer :: tsimp
double precision, allocatable, dimension(:,:) :: viscce
double precision, dimension(:,:), allocatable :: vect
double precision, dimension(:), pointer :: crom, croma, cromaa, pcrom
double precision, dimension(:), pointer :: brom_eos, crom_eos, brom, broma
double precision, dimension(:), allocatable, target :: cproa_rho_tc
double precision, dimension(:), pointer :: coefa_k, coefb_k
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:,:), allocatable :: rij
double precision, dimension(:,:), pointer :: coefap
double precision, dimension(:,:,:), pointer :: coefbp
double precision, dimension(:,:), allocatable :: coefat
double precision, dimension(:,:,:), allocatable :: coefbt
double precision, dimension(:,:), allocatable :: tflmas, tflmab
double precision, allocatable, dimension(:,:), target :: divt
double precision, dimension(:,:), pointer :: forbr, c_st_vel
double precision, dimension(:), pointer :: cvar_pr, cvara_k
double precision, dimension(:,:), pointer :: cvara_rij
double precision, dimension(:), pointer :: viscl, visct, c_estim
double precision, dimension(:,:), pointer :: lapla, lagr_st_vel
double precision, dimension(:,:), pointer :: cpro_gradp
double precision, dimension(:,:), pointer :: cpro_divr
double precision, dimension(:,:), pointer :: cpro_pred_vel
double precision, dimension(:), pointer :: cpro_wgrec_s, wgrec_crom
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: imasfl_prev, bmasfl_prev
double precision, dimension(:), pointer :: cpro_beta, cvar_t
double precision, dimension(:), allocatable, target :: cpro_rho_tc
double precision, dimension(:), pointer :: cpro_rho_mass
double precision, dimension(:,:), allocatable, target :: stf
double precision, dimension(:), pointer :: ipro_rusanov
double precision, dimension(:), pointer :: bpro_rusanov

type(var_cal_opt) :: vcopt_p, vcopt_u, vcopt
type(var_cal_opt), target :: vcopt_loc
type(var_cal_opt), pointer :: p_k_value
type(c_ptr) :: c_k_value

!===============================================================================
! 1. Initialization
!===============================================================================

! Id of the mass flux
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)

! Pointers to the mass fluxes
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Pointers to properties
! Density at time n+1,iteration iterns+1
call field_get_val_s(icrom, crom_eos)
call field_get_val_s(ibrom, brom_eos)

! Density at time (n)
if (irovar.eq.1) then
  call field_get_val_prev_s(icrom, croma)
  call field_get_val_prev_s(ibrom, broma)
else
  croma => crom_eos
  broma => brom_eos
endif

! Density at time (n-1) if needed
if (idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3.or.irovar.eq.1) then
  call field_get_val_prev2_s(icrom, cromaa)
endif

! Add Rusanov
if (irijnu.eq.2) then
  call field_get_val_s_by_name('i_rusanov_diff', ipro_rusanov)
endif

! Density for the unsteady term (at time n)
! Compressible algorithm (mass equation is already solved)
! or Low Mach compressible algos with mass flux prediction
if ((    (ippmod(icompf).ge.0.and.ippmod(icompf).ne.3)   &
     .or.(idilat.gt.1.and.ipredfl.eq.1))                 &
    .and.irovar.eq.1) then
  pcrom => croma

! VOF algorithm and Low Mach compressible algos: density at time n-1
else if ((idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3) &
         .and.irovar.eq.1) then
  ! Time-segregated solver
  if (itpcol.eq.0.and.iterns.eq.1) then
    pcrom => cromaa
  else
    pcrom => croma
  endif

! Weakly variable density algo. (idilat <=1) or constant density
else
  pcrom => crom_eos
endif

! Allocate temporary arrays
allocate(smbr(3,ncelet))
allocate(fimp(3,3,ncelet))
allocate(fimpcp(3,3,ncelet))
allocate(w1(ncelet))
call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

! Density for other terms such as buoyancy term
! 2nd order in time
if (vcopt_u%thetav .lt. 1.d0) then
  ! map the density pointer:
  ! 1/4(n-1) + 1/2(n) + 1/4(n+1)
  ! here replaced by (n)
  crom => croma
  brom => broma
! 1st order in time
else
  crom => crom_eos
  brom => brom_eos
endif

! Interpolation of rho^n-1/2 (stored in pcrom)
! Interpolation of the mass flux at (n+1/2)
! NB: the mass flux (n+1) is overwritten because not used after.
! The mass flux for (n->n+1) will be recomputed in cs_pressure_correction
! FIXME irovar=1 and if dt varies, use theta(rho) = theta(u)*...
if (vcopt_u%thetav .lt. 1.d0 .and. iappel.eq.1 .and. iterns.gt.1   &
    .and. itpcol .eq. 0) then
  allocate(cproa_rho_tc(ncelet))

  ! Pointer to the previous mass fluxes
  call field_get_val_prev_s(iflmas, imasfl_prev)
  call field_get_val_prev_s(iflmab, bmasfl_prev)

  if (irovar.eq.1) then
    ! remap the density pointer: n-1/2
    do iel = 1, ncelet
      cproa_rho_tc(iel) =  vcopt_u%thetav * croma(iel)              &
                         + (1.d0 - vcopt_u%thetav) * cromaa(iel)
    enddo
    pcrom => cproa_rho_tc
  endif

  ! Inner mass flux interpolation: n-1/2->n+1/2
  do ifac = 1, nfac
    imasfl(ifac) =  vcopt_u%thetav * imasfl(ifac)                   &
                  + (1.d0 - vcopt_u%thetav) * imasfl_prev(ifac)
  enddo

  ! Boundary mass flux interpolation: n-1/2->n+1/2
  do ifac = 1, nfabor
    bmasfl(ifac) =  vcopt_u%thetav * bmasfl(ifac)                   &
                  + (1.d0 - vcopt_u%thetav) * bmasfl_prev(ifac)
  enddo

endif

idftnp = vcopt_u%idften
imvisp = vcopt_u%imvisf

if (iand(idftnp, ANISOTROPIC_LEFT_DIFFUSION).ne.0) allocate(viscce(6,ncelet))

call field_get_id_try("est_error_pre_2", iespre)
call field_get_id_try("est_error_cor_2", iescor)
call field_get_id_try("est_error_tot_2", iestot)

! Allocate a temporary array for the prediction-stage error estimator
if (iespre.ge.0) then
  allocate(eswork(3,ncelet))
endif

if (iappel.eq.2) then
  if (iforbr.ge.0 .and. iterns.eq.1 .or. ivofmt.gt.0) then
    call field_get_val_s(ivarfl(ipr), cvar_pr)
  endif
  if(iforbr.ge.0 .and. iterns.eq.1                                          &
     .and. (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) .and. igrhok.eq.1) then
    call field_get_val_s(ivarfl(ik), cvara_k)
  endif
  if (itytur.eq.3.and.iterns.eq.1) then
    call field_get_val_v(ivarfl(irij), cvara_rij)
  endif
else
  if (iforbr.ge.0 .and. iterns.eq.1 .or. ivofmt.gt.0) then
    call field_get_val_s(ivarfl(ipr), cvar_pr)
  endif
  if(iforbr.ge.0 .and. iterns.eq.1                                          &
     .and. (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) .and. igrhok.eq.1) then
      call field_get_val_prev_s(ivarfl(ik), cvara_k)
  endif
  if (itytur.eq.3.and.iterns.eq.1) then
    call field_get_val_v(ivarfl(irij), cvara_rij)
  endif
endif

if (iforbr.ge.0 .and. iterns.eq.1) call field_get_val_v(iforbr, forbr)

! Theta relatif aux termes sources explicites
thets  = thetsn
if (isno2t.gt.0) then
  call field_get_key_int(ivarfl(iu), kstprv, f_id)
  call field_get_val_v(f_id, c_st_vel)
else
  c_st_vel => null()
endif

!-------------------------------------------------------------------------------
! ---> Get user source terms

call field_get_id_try("velocity_source_term_exp", f_id)
if (f_id.ge.0) then
  call field_get_val_v(f_id, tsexp)
else
  allocate(loctsexp(3,ncelet))
  tsexp => loctsexp
endif

call field_get_id_try("velocity_source_term_imp", f_id)
if (f_id.ge.0) then
  call field_get_val_t(f_id, tsimp)
else
  allocate(loctsimp(3,3,ncelet))
  tsimp => loctsimp
endif

do iel = 1, ncel
  do isou = 1, 3
    tsexp(isou,iel) = 0.d0
    do jsou = 1, 3
      tsimp(jsou,isou,iel) = 0.d0
    enddo
  enddo
enddo

! The computation of explicit and implicit source terms is performed
! at the first iteration only.
! If iphydr=1 or if we have buoyant scalars then we need to update source terms
call uitsnv (vel, tsexp, tsimp)

! C version
call user_source_terms(ivarfl(iu), tsexp, tsimp)

! Boundary contributions for all immersed boundaries
if (iporos.eq.3) then
  call cs_immersed_boundary_wall_functions(ivarfl(iu), tsexp, tsimp)
endif

n_fans = cs_fan_n_fans()
if (n_fans .gt. 0) then
  if (ntcabs.eq.ntpabs+1) then
    call debvtl(imasfl, bmasfl, crom, brom)
  endif
  call tsvvtl(tsexp)
endif

! Skip first time step after restart if previous values have not been read.
if (vcopt_u%ibdtso.lt.0) then
  vcopt_u%ibdtso = iabs(vcopt_u%ibdtso)
  call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
endif

! Nudging towards optimal interpolation for velocity
if (ippmod(iatmos).ge.0) then
  call field_get_key_int(ivarfl(iu), kopint, f_oi_id)
  if (f_oi_id.ge.0) then
    call cs_at_data_assim_source_term(ivarfl(iu), tsexp, tsimp)
  endif
  if (iatmst.ge.1) then
    call cs_at_source_term_for_inlet(tsexp)
  endif
endif

! Coupling between two code_saturne
if (nbrcpl.gt.0) then
  !vectorial interleaved exchange
  call csccel(ivarfl(iu), tsexp)
endif

if (vcopt_u%ibdtso.gt.1.and.ntcabs.gt.ntinit &
  .and.(idtvar.eq.0.or.idtvar.eq.1)) then
  ! TODO: remove test on ntcabs and implemente a "proper" condition for
  ! initialization.
  f_id = ivarfl(iu)
  call cs_backward_differentiation_in_time(f_id, tsexp, tsimp)
endif

!===============================================================================
! 2. Potential forces (pressure gradient and gravity)
!===============================================================================

! Pressure gradient

call field_get_id_try("pressure_gradient", f_id)
if (f_id.ge.0) then
  call field_get_val_v(f_id, cpro_gradp)
else
  allocate(grad(3,ncelet))
  cpro_gradp => grad
endif

inc    = 1

! Take the latest pressure field
iprev = 0

! Namely for the VOF algorithm: consistency of the gradient
! with the diffusive flux scheme of the correction step
if (vcopt_p%iwgrec.eq.1) then

  ! retrieve density used in diffusive flux scheme (correction step)
  if (irovar.eq.1.and.(idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3)) then
    call field_get_id("density_mass", f_id)
    call field_get_val_s(f_id, cpro_rho_mass)

    ! Time interpolated density
    if (vcopt_u%thetav.lt.1.d0.and.iterns.gt.1) then
      allocate(cpro_rho_tc(ncelet))

      do iel = 1, ncelet
        cpro_rho_tc(iel) =           vcopt_u%thetav * cpro_rho_mass(iel) &
                          + (1.d0 - vcopt_u%thetav) * croma(iel)
      enddo

      wgrec_crom => cpro_rho_tc
    else
      wgrec_crom => cpro_rho_mass
    endif

  ! Weakly variable density algo. (idilat <=1) or constant density
  else
    wgrec_crom => crom_eos
  endif

  ! Id weighting field for gradient
  call field_get_key_int(ivarfl(ipr), kwgrec, iflwgr)
  call field_get_dim(iflwgr, f_dim)
  if (f_dim.gt.1) then
    call field_get_val_v(iflwgr, cpro_wgrec_v)
    do iel = 1, ncel
      ! FIXME should take head losses into account,
      ! not compatible either with ipucou=1...
      cpro_wgrec_v(1,iel) = dt(iel) / wgrec_crom(iel)
      cpro_wgrec_v(2,iel) = dt(iel) / wgrec_crom(iel)
      cpro_wgrec_v(3,iel) = dt(iel) / wgrec_crom(iel)
      cpro_wgrec_v(4,iel) = 0.d0
      cpro_wgrec_v(5,iel) = 0.d0
      cpro_wgrec_v(6,iel) = 0.d0
    enddo
    call syntis(cpro_wgrec_v)

  else
    call field_get_val_s(iflwgr, cpro_wgrec_s)
    do iel = 1, ncel
      cpro_wgrec_s(iel) = dt(iel) / wgrec_crom(iel)
    enddo
    call synsca(cpro_wgrec_s)
  endif
  if (allocated(cpro_rho_tc)) deallocate(cpro_rho_tc)
endif

call grdpor(inc)

! Pressure gradient
if (iprcdo.eq.0) then
  call field_gradient_potential(ivarfl(ipr), iprev, inc,      &
                                iphydr,                       &
                                frcxt, cpro_gradp)
endif

!    Calcul des efforts aux parois (partie 2/5), si demande
!    La pression a la face est calculee comme dans gradrc/gradmc
!    et on la transforme en pression totale
!    On se limite a la premiere iteration (pour faire simple par
!      rapport a la partie issue de condli, hors boucle)
if (iforbr.ge.0 .and. iterns.eq.1) then
  call field_get_coefa_s (ivarfl(ipr), coefa_p)
  call field_get_coefb_s (ivarfl(ipr), coefb_p)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)
    pip = cvar_pr(iel)                          &
          + diipbx*cpro_gradp(1,iel)            &
          + diipby*cpro_gradp(2,iel)            &
          + diipbz*cpro_gradp(3,iel)
    pfac = coefa_p(ifac) + coefb_p(ifac)*pip
    pfac = pfac                                                   &
         + ro0*(gx*(cdgfbo(1,ifac)-xyzp0(1))                      &
         + gy*(cdgfbo(2,ifac)-xyzp0(2))                           &
         + gz*(cdgfbo(3,ifac)-xyzp0(3)) )                         &
         - pred0
    do isou = 1, 3
      forbr(isou,ifac) = forbr(isou,ifac) + pfac*surfbo(isou,ifac)
    enddo
  enddo
endif

if (iappel.eq.1) then

  ! Initialization
  ! NB: at the second call, trav contains the temporal increment
  do iel = 1, ncel
    trav(1,iel) = 0.d0
    trav(2,iel) = 0.d0
    trav(3,iel) = 0.d0
  enddo
endif

! FIXME : "rho g" will be second order only if extrapolated
if (iphydr.eq.1) then
  do iel = 1, ncel
    trav(1,iel) = trav(1,iel)+(frcxt(1 ,iel) - cpro_gradp(1,iel)) * cell_f_vol(iel)
    trav(2,iel) = trav(2,iel)+(frcxt(2 ,iel) - cpro_gradp(2,iel)) * cell_f_vol(iel)
    trav(3,iel) = trav(3,iel)+(frcxt(3 ,iel) - cpro_gradp(3,iel)) * cell_f_vol(iel)
  enddo
else if (iphydr.eq.2) then
  do iel = 1, ncel
    rom = crom(iel)
    trav(1,iel) = trav(1,iel) +  (rom*gx - grdphd(1,iel) - cpro_gradp(1,iel)) &
                                * cell_f_vol(iel)
    trav(2,iel) = trav(2,iel) +  (rom*gy - grdphd(2,iel) - cpro_gradp(2,iel)) &
                                * cell_f_vol(iel)
    trav(3,iel) = trav(3,iel) +  (rom*gz - grdphd(3,iel) - cpro_gradp(3,iel)) &
                                * cell_f_vol(iel)
  enddo
else if (ippmod(icompf).ge.0) then
  do iel = 1, ncel
    rom = crom(iel)
    trav(1,iel) = trav(1,iel)+(rom*gx - cpro_gradp(1,iel)) * cell_f_vol(iel)
    trav(2,iel) = trav(2,iel)+(rom*gy - cpro_gradp(2,iel)) * cell_f_vol(iel)
    trav(3,iel) = trav(3,iel)+(rom*gz - cpro_gradp(3,iel)) * cell_f_vol(iel)
  enddo
  ! Boussinesq approximation
else if (idilat.eq.0) then

  !FIXME make it dependant on the scalar and use is_buoyant field
  call field_get_val_s(ibeta, cpro_beta)
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_t)

  ! Delta rho = - rho_0 beta (T-T0)
  tref = t0
  ! for atmospheric flows, variable is potential temperature
  if (ippmod(iatmos).gt.0) then
    rscp = rair/cp0
    tref = t0 * (ps / p0)**rscp
  endif

  do iel = 1, ncel
    drom = - crom(iel) * cpro_beta(iel) * (cvar_t(iel) - tref)
    trav(1,iel) = trav(1,iel) + (drom * gx - cpro_gradp(1,iel)) * cell_f_vol(iel)
    trav(2,iel) = trav(2,iel) + (drom * gy - cpro_gradp(2,iel)) * cell_f_vol(iel)
    trav(3,iel) = trav(3,iel) + (drom * gz - cpro_gradp(3,iel)) * cell_f_vol(iel)
  enddo

else
  do iel = 1, ncel
    if (ischtp.eq.2 .and. itpcol.eq.1) then
      drom = (3.d0/2.d0*croma(iel) - 1.d0/2.d0*cromaa(iel) - ro0)
    else
      drom = (crom(iel)-ro0)
    endif
    trav(1,iel) = trav(1,iel)+(drom*gx - cpro_gradp(1,iel) ) * cell_f_vol(iel)
    trav(2,iel) = trav(2,iel)+(drom*gy - cpro_gradp(2,iel) ) * cell_f_vol(iel)
    trav(3,iel) = trav(3,iel)+(drom*gz - cpro_gradp(3,iel) ) * cell_f_vol(iel)
  enddo
endif

! Free memory
if (allocated(grad)) deallocate(grad)


!   Pour IAPPEL = 1 (ie appel standard sans les estimateurs)
!     TRAV rassemble les termes sources  qui seront recalcules
!       a toutes les iterations sur navsto
!     Si on n'itere pas sur navsto et qu'on n'extrapole pas les
!       termes sources, TRAV contient tous les termes sources
!       jusqu'au basculement dans SMBR
!     A ce niveau, TRAV contient -grad P et rho g
!       P est suppose pris a n+1/2
!       rho est eventuellement interpole a n+1/2

!-------------------------------------------------------------------------------
! ---> Initialize trava array and source terms at the first call (iterns=1)

!     trava contains all source terms needed from the first sub iteration
!       (iterns=1) for the other iterations.
!     When there is only one iteration, we build source terms directly in trav
!       array.
!     Explicit source terms will be used at the next time step in case of
!       extrapolation (if there is only one or many iteration on navtsv)

! At the first iteration on navstv
if (iterns.eq.1) then

  ! Si on   extrapole     les T.S. : -theta*valeur precedente
  if (isno2t.gt.0) then
    ! S'il n'y a qu'une    iter : TRAV  incremente
    if (nterup.eq.1) then
      do iel = 1, ncel
        do ii = 1, ndim
          trav (ii,iel) = trav (ii,iel) - thets*c_st_vel(ii,iel)
        enddo
      enddo
      ! S'il   y a plusieurs iter : TRAVA initialise
    else
      do iel = 1, ncel
        do ii = 1, ndim
          trava(ii,iel) = - thets*c_st_vel(ii,iel)
        enddo
      enddo
    endif
    ! Et on initialise le terme source pour le remplir ensuite
    do iel = 1, ncel
      do ii = 1, ndim
        c_st_vel(ii,iel) = 0.d0
      enddo
    enddo

  ! Si on n'extrapole pas les T.S.
  else
    ! S'il   y a plusieurs iter : TRAVA initialise
    !  sinon TRAVA n'existe pas
    if(nterup.gt.1) then
      do ii = 1, ndim
        do iel = 1, ncel
          trava(ii,iel)  = 0.d0
        enddo
      enddo
    endif
  endif

endif

!-------------------------------------------------------------------------------
! Initialization of the implicit terms

if (iappel.eq.1) then

  do iel = 1, ncel
    do isou = 1, 3
      fimp(isou,isou,iel) = vcopt_u%istat*pcrom(iel)/dt(iel)*cell_f_vol(iel)
      do jsou = 1, 3
        if(jsou.ne.isou) fimp(jsou,isou,iel) = 0.d0
      enddo
    enddo
  enddo

!     Le remplissage de FIMP est toujours indispensable,
!       meme si on peut se contenter de n'importe quoi pour IAPPEL=2.
else
  do iel = 1, ncel
    do isou = 1, 3
      do jsou = 1, 3
        fimp(jsou,isou,iel) = 0.d0
      enddo
    enddo
  enddo
endif

!-------------------------------------------------------------------------------
! ---> 2/3 RHO * GRADIENT DE K SI k-epsilon ou k-omega
!      NB : ON NE PREND PAS LE GRADIENT DE (RHO K), MAIS
!           CA COMPLIQUERAIT LA GESTION DES CL ...
!     On peut se demander si l'extrapolation en temps sert a
!       quelquechose

!     Ce terme explicite est calcule une seule fois,
!       a la premiere iter sur navsto : il est stocke dans un champ si on
!       doit l'extrapoler en temps ; il va dans TRAVA si on n'extrapole
!       pas mais qu'on itere sur navsto. Il va dans TRAV si on
!       n'extrapole pas et qu'on n'itere pas sur navsto.
if(     (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) &
   .and. igrhok.eq.1 .and. iterns.eq.1) then

  ! Allocate a work array for the gradient calculation
  allocate(grad(3,ncelet))

  iprev  = 1
  inc    = 1

  call field_gradient_scalar(ivarfl(ik), iprev, inc, grad)

  d2s3 = 2.d0/3.d0

  ! Si on extrapole les termes source en temps
  if (isno2t.gt.0) then
    ! Calcul de rho^n grad k^n      si rho non extrapole
    !           rho^n grad k^n      si rho     extrapole

    do iel = 1, ncel
      romvom = -croma(iel)*cell_f_vol(iel)*d2s3
      do isou = 1, 3
        c_st_vel(isou,iel) = c_st_vel(isou,iel)+grad(isou,iel)*romvom
      enddo
    enddo
  ! Si on n'extrapole pas les termes sources en temps : TRAV ou TRAVA
  else
    if(nterup.eq.1) then
      do iel = 1, ncel
        romvom = -crom(iel)*cell_f_vol(iel)*d2s3
        do isou = 1, 3
          trav(isou,iel) = trav(isou,iel) + grad(isou,iel) * romvom
        enddo
      enddo
    else
      do iel = 1, ncel
        romvom = -crom(iel)*cell_f_vol(iel)*d2s3
        do isou = 1, 3
          trava(isou,iel) = trava(isou,iel) + grad(isou,iel) * romvom
        enddo
      enddo
    endif
  endif

  ! Calcul des efforts aux parois (partie 3/5), si demande
  if (iforbr.ge.0) then
    call field_get_coefa_s (ivarfl(ik), coefa_k)
    call field_get_coefb_s (ivarfl(ik), coefb_k)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      xkb = cvara_k(iel) + diipbx*grad(1,iel)                      &
           + diipby*grad(2,iel) + diipbz*grad(3,iel)
      xkb = coefa_k(ifac)+coefb_k(ifac)*xkb
      xkb = d2s3*crom(iel)*xkb
      do isou = 1, 3
        forbr(isou,ifac) = forbr(isou,ifac) + xkb*surfbo(isou,ifac)
      enddo
    enddo
  endif

  ! Free memory
  deallocate(grad)

endif

!-------------------------------------------------------------------------------
! ---> Transpose of velocity gradient in the diffusion term

!     These terms are taken into account in cs_balance_vector.
!     We only compute here the secondary viscosity.

if (ivisse.eq.1) then

  call visecv(secvif, secvib)

endif

!-------------------------------------------------------------------------------
! ---> Head losses
!      (if iphydr=1 this term has already been taken into account)

! ---> Explicit part
if ((ncepdc.gt.0).and.(iphydr.ne.1)) then

  ! Les termes diagonaux sont places dans TRAV ou TRAVA,
  !   La prise en compte de velk a partir de la seconde iteration
  !   est faite directement dans cs_equation_iterative_solve_vector.
  if (iterns.eq.1) then

    allocate(hl_exp(3, ncepdc))

    call tspdcv(ncepdc, icepdc, vela, ckupdc, hl_exp)

    ! If we have inner iterations, we use trava, otherwise trav
    if(nterup.gt.1) then
      do ielpdc = 1, ncepdc
        iel    = icepdc(ielpdc)
        trava(1,iel) = trava(1,iel) + hl_exp(1,ielpdc)
        trava(2,iel) = trava(2,iel) + hl_exp(2,ielpdc)
        trava(3,iel) = trava(3,iel) + hl_exp(3,ielpdc)
      enddo
    else
      do ielpdc = 1, ncepdc
        iel    = icepdc(ielpdc)
        trav(1,iel) = trav(1,iel) + hl_exp(1,ielpdc)
        trav(2,iel) = trav(2,iel) + hl_exp(2,ielpdc)
        trav(3,iel) = trav(3,iel) + hl_exp(3,ielpdc)
      enddo
    endif

    deallocate(hl_exp)
  endif

endif

! ---> Implicit part

!  At the second call, fimp is not needed anymore
if (iappel.eq.1) then
  if (ncepdc.gt.0) then
    ! The theta-scheme for the head loss is the same as the other terms
    thetap = vcopt_u%thetav
    do ielpdc = 1, ncepdc
      iel = icepdc(ielpdc)
      romvom = crom(iel)*cell_f_vol(iel)*thetap

      ! Diagonal part
      do isou = 1, 3
        fimp(isou,isou,iel) = fimp(isou,isou,iel) + romvom*ckupdc(isou,ielpdc)
      enddo
      ! Extra-diagonal part
      cpdc12 = ckupdc(4,ielpdc)
      cpdc23 = ckupdc(5,ielpdc)
      cpdc13 = ckupdc(6,ielpdc)

      fimp(1,2,iel) = fimp(1,2,iel) + romvom*cpdc12
      fimp(2,1,iel) = fimp(2,1,iel) + romvom*cpdc12
      fimp(1,3,iel) = fimp(1,3,iel) + romvom*cpdc13
      fimp(3,1,iel) = fimp(3,1,iel) + romvom*cpdc13
      fimp(2,3,iel) = fimp(2,3,iel) + romvom*cpdc23
      fimp(3,2,iel) = fimp(3,2,iel) + romvom*cpdc23
    enddo
  endif
endif

!-------------------------------------------------------------------------------
! ---> Surface tension force for VoF

if (ivofmt.gt.0.and.sigmaS.gt.0.d0) then
  allocate(stf(3,ncel))
  call vof_surface_tension(stf)
endif

!-------------------------------------------------------------------------------
! ---> Coriolis force
!     (if iphydr=1 then this term is already taken into account)

! --->  Explicit part

if ((icorio.eq.1.or.iturbo.eq.1) .and. iphydr.ne.1) then

  ! A la premiere iter sur navsto, on ajoute la partie issue des
  ! termes explicites
  if (iterns.eq.1) then

    ! Si on n'itere pas sur navsto : TRAV
    if (nterup.eq.1) then

      ! Reference frame rotation
      do iel = 1, ncel
        romvom = -2.d0*crom(iel)*cell_f_vol(iel)
        call add_coriolis_v(0, romvom, vela(:,iel), trav(:,iel))
      enddo
      ! Turbomachinery frozen rotors rotation
      if (iturbo.eq.1) then
        do iel = 1, ncel
          if (irotce(iel).gt.0) then
            romvom = -crom(iel)*cell_f_vol(iel)
            call add_coriolis_v(irotce(iel), romvom, vela(:,iel), trav(:,iel))
          endif
        enddo
      endif

    ! Si on itere sur navsto : TRAVA
    else

      ! Reference frame rotation
      do iel = 1, ncel
        romvom = -2.d0*crom(iel)*cell_f_vol(iel)
        call add_coriolis_v(0, romvom, vela(:,iel), trava(:,iel))
      enddo
      ! Turbomachinery frozen rotors rotation
      if (iturbo.eq.1) then
        do iel = 1, ncel
          if (irotce(iel).gt.0) then
            romvom = -crom(iel)*cell_f_vol(iel)
            call add_coriolis_v(irotce(iel), romvom, vela(:,iel), trava(:,iel))
          endif
        enddo
      endif

    endif
  endif
endif

! --->  Implicit part

!  At the second call, fimp is not needed anymore
if (iappel.eq.1) then
  if (icorio.eq.1 .or. iturbo.eq.1) then
    ! The theta-scheme for the Coriolis term is the same as the other terms
    thetap = vcopt_u%thetav

    ! Reference frame rotation
    do iel = 1, ncel
      romvom = -2.d0*crom(iel)*cell_f_vol(iel)*thetap
      call add_coriolis_t(0, romvom, fimp(:,:,iel))
    enddo
    ! Turbomachinery frozen rotors rotation
    if (iturbo.eq.1) then
      do iel = 1, ncel
        if (irotce(iel).gt.0) then
          romvom = -crom(iel)*cell_f_vol(iel)*thetap
          call add_coriolis_t(irotce(iel), romvom, fimp(:,:,iel))
        endif
      enddo
    endif
  endif
endif

!-------------------------------------------------------------------------------
! ---> - Divergence of tensor Rij
! ---> - Non linear part of Rij for non-liear Eddy Viscosity Models

if((itytur.eq.3.or.iturb.eq.23.or.iturb.eq.24).and.iterns.eq.1) then

  allocate(rij(6,ncelet))
  allocate(coefat(6,nfabor))
  allocate(coefbt(6,6,nfabor))

  ! Reynolds Stress Models
  if (itytur.eq.3) then

    do iel = 1, ncelet
      rij(1,iel) = cvara_rij(1,iel)
      rij(2,iel) = cvara_rij(2,iel)
      rij(3,iel) = cvara_rij(3,iel)
      rij(4,iel) = cvara_rij(4,iel)
      rij(5,iel) = cvara_rij(5,iel)
      rij(6,iel) = cvara_rij(6,iel)
    enddo
    ! --- Boundary conditions on the components of the tensor Rij
    call field_get_coefad_v(ivarfl(irij),coefap)
    coefat = coefap

    do ifac = 1, nfabor
      do ii = 1, 6
        do jj = 1, 6
          coefbt(jj,ii,ifac) = 0.d0
        enddo
      enddo
    enddo

    call field_get_coefbd_v(ivarfl(irij),coefbp)
    coefbt = coefbp

  ! Baglietto et al. quadratic k-epislon model
  else if(iturb.eq.23) then

    ! --- Compute the non linear part of Rij
    call cs_turbulence_ke_q(rij)
    ! --- Boundary conditions : Homogeneous Neumann
    do ifac = 1, nfabor
      do ii = 1, 6
        coefat(ii,ifac) = 0.d0
        do jj = 1, 6
          coefbt(jj,ii,ifac) = 1.d0
        enddo
      enddo
    enddo
  ! Baglietto et al. cubic k-epsilon model
  else if(iturb.eq.24) then

    ! --- Compute the non linear part of Rij
    call cs_turbulence_ke_c(rij)
    ! --- Boundary conditions : Homogeneous Neumann
    do ifac = 1, nfabor
      do ii = 1, 6
        coefat(ii,ifac) = 0.d0
        do jj = 1, 6
          coefbt(jj,ii,ifac) = 1.d0
        enddo
      enddo
    enddo
  end if

  ! Flux computation options
  f_id = -1
  init = 1
  inc  = 1
  iflmb0 = 0
  if(itytur.eq.3) then
    call field_get_key_struct_var_cal_opt(ivarfl(irij), vcopt)
  else
    call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
  end if
  imrgrp = vcopt%imrgra
  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  itypfl = 1

  allocate(tflmas(3,nfac))
  allocate(tflmab(3,nfabor))

  call divrij &
 ( f_id   , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgrp , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   rij    ,                                                       &
   coefat , coefbt ,                                              &
   tflmas , tflmab )

  deallocate(rij)
  deallocate(coefat, coefbt)

  !     Calcul des efforts aux bords (partie 5/5), si necessaire

  if (iforbr.ge.0) then
    do ifac = 1, nfabor
      do isou = 1, 3
        forbr(isou,ifac) = forbr(isou,ifac) + tflmab(isou,ifac)
      enddo
    enddo
  endif


  call field_get_id_try("reynolds_stress_divergence", f_id)
  if (f_id.ge.0) then
    call field_get_val_v(f_id, cpro_divr)
  else
    allocate(divt(3,ncelet))
    cpro_divr => divt
  endif

  init = 1
  call divmat(init,tflmas,tflmab, cpro_divr)

  deallocate(tflmas, tflmab)

  ! (if iphydr=1 then this term is already taken into account)
  if (iphydr.ne.1.or.igprij.ne.1) then

    ! If extrapolation of source terms
    if (isno2t.gt.0) then
      do iel = 1, ncel
        do isou = 1, 3
          c_st_vel(isou,iel) = c_st_vel(isou,iel) - cpro_divr(isou,iel)
        enddo
      enddo

    ! No extrapolation of source terms
    else

      ! No inner iteration
      if (nterup.eq.1) then
        do iel = 1, ncel
          do isou = 1, 3
            trav(isou,iel) = trav(isou,iel) - cpro_divr(isou,iel)
          enddo
        enddo
      ! Inner iterations
      else
        do iel = 1, ncel
          do isou = 1, 3
            trava(isou,iel) = trava(isou,iel) - cpro_divr(isou,iel)
          enddo
        enddo
      endif
    endif
  endif

endif

!-------------------------------------------------------------------------------
! ---> Face diffusivity for the velocity

if (vcopt_u%idiff.ge. 1) then

  call field_get_val_s(iviscl, viscl)
  call field_get_val_s(ivisct, visct)

  if (itytur.eq.3) then
    do iel = 1, ncel
      w1(iel) = viscl(iel)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = viscl(iel) + vcopt_u%idifft*visct(iel)
    enddo
  endif

  ! Scalar diffusivity (Default)
  if (iand(idftnp, ISOTROPIC_DIFFUSION).ne.0) then

    call viscfa &
   ( imvisp ,                                                       &
     w1     ,                                                       &
     viscf  , viscb  )

    ! When using Rij-epsilon model with the option irijnu=1, the face
    ! viscosity for the Matrix (viscfi and viscbi) is increased
    if (itytur.eq.3.and.irijnu.eq.1) then

      do iel = 1, ncel
        w1(iel) = viscl(iel) + vcopt_u%idifft*visct(iel)
      enddo

      call viscfa &
   ( imvisp ,                                                       &
     w1     ,                                                       &
     viscfi , viscbi )

    endif

  ! Tensorial diffusion of the velocity (in case of tensorial porosity)
  else if (iand(idftnp, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then

    do iel = 1, ncel
      do isou = 1, 3
        viscce(isou, iel) = w1(iel)
      enddo
      do isou = 4, 6
        viscce(isou, iel) = 0.d0
      enddo
    enddo

    call vistnv &
     ( imvisp ,                                                       &
       viscce ,                                                       &
       viscf  , viscb  )

    ! When using Rij-epsilon model with the option irijnu=1, the face
    ! viscosity for the Matrix (viscfi and viscbi) is increased
    if (itytur.eq.3.and.irijnu.eq.1) then

      do iel = 1, ncel
        w1(iel) = viscl(iel) + vcopt_u%idifft*visct(iel)
      enddo

      do iel = 1, ncel
        do isou = 1, 3
          viscce(isou, iel) = w1(iel)
        enddo
        do isou = 4, 6
          viscce(isou, iel) = 0.d0
        enddo
      enddo

      call vistnv &
       ( imvisp ,                                                       &
         viscce ,                                                       &
         viscfi , viscbi )

    endif
  endif


! --- If no diffusion, viscosity is set to 0.
else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  if(itytur.eq.3.and.irijnu.eq.1) then
    do ifac = 1, nfac
      viscfi(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscbi(ifac) = 0.d0
    enddo
  endif

endif

if (irijnu.eq.2) then
  ! Add Rusanov
  do ifac = 1, nfac
    viscf(ifac) = viscf(ifac) + ipro_rusanov(ifac)
  enddo

  call field_get_val_s_by_name('b_rusanov_diff', bpro_rusanov)
  do ifac = 1, nfabor
    n(1) = surfbo(1,ifac) / surfbn(ifac) ! WARNING Normalised here
    n(2) = surfbo(2,ifac) / surfbn(ifac)
    n(3) = surfbo(3,ifac) / surfbn(ifac)

    do isou = 1, 3
      do jsou = 1, 3
        cofbfv(jsou, isou, ifac) = cofbfv(jsou, isou, ifac) &
          + 2.d0 * bpro_rusanov(ifac) * n(isou)*n(jsou)
      enddo
    enddo
  enddo
endif

!-------------------------------------------------------------------------------
! ---> Take external forces partially equilibrated with the pressure gradient
!      into account (only for the first call, the second one is dedicated
!      to error estimators)

if (iappel.eq.1.and.iphydr.eq.1) then

  ! External forces at previous time step:
  !     frcxt was initialised to 0
  !     NB: frcxt was used in typecl, and will be updated
  !         at the end of navstv

  ! External force variation between time step n and n+1
  ! (used in the correction step)
  !-----------------------------

  ! Boussinesq approximation
  if (idilat.eq.0) then

    !FIXME make it dependant on the scalar and use is_buoyant field
    call field_get_val_s(ibeta, cpro_beta)
    call field_get_val_s(ivarfl(isca(iscalt)), cvar_t)

    ! Delta rho = - rho_0 beta (T-T0)
    tref = t0
    ! for atmospheric flows, variable is potential temperature
    if (ippmod(iatmos).gt.0) then
      rscp = rair/cp0
      tref = t0 * (ps / p0)**rscp
    endif

    do iel = 1, ncel
      drom = - crom(iel) * cpro_beta(iel) * (cvar_t(iel) - tref) * cell_is_active(iel)
      dfrcxt(1, iel) = drom*gx - frcxt(1, iel) * cell_is_active(iel)
      dfrcxt(2, iel) = drom*gy - frcxt(2, iel) * cell_is_active(iel)
      dfrcxt(3, iel) = drom*gz - frcxt(3, iel) * cell_is_active(iel)
    enddo

  else
    do iel = 1, ncel
      if (ischtp.eq.2 .and. itpcol.eq.1) then
        drom = (3.d0/2.d0*croma(iel) - 1.d0/2.d0*cromaa(iel) - ro0)            &
             * cell_is_active(iel)
      else
        drom = (crom(iel)-ro0) * cell_is_active(iel)
      endif

      dfrcxt(1, iel) = drom*gx - frcxt(1, iel) * cell_is_active(iel)
      dfrcxt(2, iel) = drom*gy - frcxt(2, iel) * cell_is_active(iel)
      dfrcxt(3, iel) = drom*gz - frcxt(3, iel) * cell_is_active(iel)
    enddo
  endif

  ! Add head losses
  if (ncepdc.gt.0) then
    do ielpdc = 1, ncepdc
      iel=icepdc(ielpdc)
      vit1   = vela(1,iel) * cell_is_active(iel)
      vit2   = vela(2,iel) * cell_is_active(iel)
      vit3   = vela(3,iel) * cell_is_active(iel)
      cpdc11 = ckupdc(1,ielpdc)
      cpdc22 = ckupdc(2,ielpdc)
      cpdc33 = ckupdc(3,ielpdc)
      cpdc12 = ckupdc(4,ielpdc)
      cpdc23 = ckupdc(5,ielpdc)
      cpdc13 = ckupdc(6,ielpdc)
      dfrcxt(1 ,iel) = dfrcxt(1 ,iel) &
                    - crom(iel)*(cpdc11*vit1+cpdc12*vit2+cpdc13*vit3)
      dfrcxt(2 ,iel) = dfrcxt(2 ,iel) &
                    - crom(iel)*(cpdc12*vit1+cpdc22*vit2+cpdc23*vit3)
      dfrcxt(3 ,iel) = dfrcxt(3 ,iel) &
                    - crom(iel)*(cpdc13*vit1+cpdc23*vit2+cpdc33*vit3)
    enddo
  endif

  ! Add Coriolis force
  if (icorio.eq.1 .or. iturbo.eq.1) then

    ! Reference frame rotation
    do iel = 1, ncel
      rom = -2.d0*crom(iel) * cell_is_active(iel)
      call add_coriolis_v(0, rom, vela(:,iel), dfrcxt(:,iel))
    enddo
    ! Turbomachinery frozen rotors rotation
    if (iturbo.eq.1) then
      do iel = 1, ncel
        if (irotce(iel).gt.0) then
          rom = -crom(iel) * cell_is_active(iel)
          call add_coriolis_v(irotce(iel), rom, vela(:,iel), dfrcxt(:,iel))
        endif
      enddo
    endif
  endif

  ! Add -div( rho R) as external force
  if (itytur.eq.3.and.igprij.eq.1) then
    do iel = 1, ncel
      dvol = 0.d0
      ! If it is not a solid cell
      if (cell_is_active(iel).eq.1) dvol = 1.d0 / cell_f_vol(iel)
      do isou = 1, 3
        dfrcxt(isou, iel) = dfrcxt(isou, iel) - cpro_divr(isou, iel)*dvol
      enddo
    enddo
  endif

  ! Surface tension force for VoF
  if (ivofmt.gt.0.and.sigmaS.gt.0.d0) then
    do iel = 1, ncel
      dvol = 0.d0
      ! If it is not a solid cell
      if (cell_is_active(iel).eq.1) dvol = 1.d0 / cell_f_vol(iel)
      do isou = 1, 3
        dfrcxt(isou, iel) = dfrcxt(isou, iel) + stf(isou, iel)*dvol
      enddo
    enddo
  endif

  ! ---> Use user source terms

  if (igpust.eq.1) then
    do iel = 1, ncel
      dvol = 0.d0
      ! If it is not a solid cell
      if (cell_is_active(iel).eq.1) dvol = 1.d0 / cell_f_vol(iel)

      do isou = 1, 3
        dfrcxt(isou, iel) = dfrcxt(isou, iel) + tsexp(isou, iel)*dvol
      enddo
    enddo

  endif

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(dfrcxt)
  endif

endif

!===============================================================================
! 3. Solving of the 3x3xNcel coupled system
!===============================================================================

! ---> AU PREMIER APPEL,
!      MISE A ZERO DE L'ESTIMATEUR POUR LA VITESSE PREDITE
!      S'IL DOIT ETRE CALCULE

if (iappel.eq.1) then
  if (iespre.ge.0) then
    call field_get_val_s(iespre, c_estim)
    do iel = 1, ncel
      c_estim(iel) =  0.d0
    enddo
  endif
endif

! ---> AU DEUXIEME APPEL,
!      MISE A ZERO DE L'ESTIMATEUR TOTAL POUR NAVIER-STOKES
!      (SI ON FAIT UN DEUXIEME APPEL, ALORS IL DOIT ETRE CALCULE)

if (iappel.eq.2) then
  call field_get_val_s(iestot, c_estim)
  do iel = 1, ncel
    c_estim(iel) =  0.d0
  enddo
endif

!-------------------------------------------------------------------------------
! ---> Use user source terms

! ---> Explicit contribution due to implicit terms

if (iterns.eq.1) then
  if (nterup.gt.1) then
    do iel = 1, ncel
      do isou = 1, 3
        do jsou = 1, 3
          trava(isou,iel) = trava(isou,iel)                                  &
                          + tsimp(jsou,isou,iel)*vela(jsou,iel)
        enddo
      enddo
    enddo
  else
    do iel = 1, ncel
      do isou = 1, 3
        do jsou = 1, 3
          trav(isou,iel) = trav(isou,iel)                                    &
                         + tsimp(jsou,isou,iel)*vela(jsou,iel)
        enddo
      enddo
    enddo
  endif
endif

! Explicit user source terms are added
if ((iphydr.ne.1.or.igpust.ne.1)) then
  ! If source terms are time-extrapolated, they are stored in fields
  if (isno2t.gt.0) then
    if (iterns.eq.1) then
      do iel = 1, ncel
        do isou = 1, 3
          c_st_vel(isou,iel) = c_st_vel(isou,iel) + tsexp(isou,iel)
        enddo
      enddo
    endif

  else
    ! Always in the current work array because this may be updated
    ! during inner iterations
     do iel = 1, ncel
       do isou = 1, 3
         trav(isou,iel) = trav(isou,iel) + tsexp(isou,iel)
       enddo
    enddo
  endif
endif

! Surface tension is added
if (iphydr.ne.1.and.sigmaS.gt.0.d0) then
  ! If source terms are time-extrapolated, they are stored in fields
  if (isno2t.gt.0) then
    if (iterns.eq.1) then
      do iel = 1, ncel
        do isou = 1, 3
          c_st_vel(isou,iel) = c_st_vel(isou,iel) + stf(isou,iel)
        enddo
      enddo
    endif

  else
    ! Always in the current work array because this may be updated
    ! during inner iterations
     do iel = 1, ncel
       do isou = 1, 3
         trav(isou,iel) = trav(isou,iel) + stf(isou,iel)
       enddo
    enddo
  endif
endif

! ---> Implicit terms
if (iappel.eq.1) then
  ! If source terms are time-extrapolated
  if (isno2t.gt.0) then
    thetap = vcopt_u%thetav
    do iel = 1, ncel
      do isou = 1, 3
        do jsou = 1, 3
          fimp(jsou,isou,iel) = fimp(jsou,isou,iel)                      &
                              - tsimp(jsou,isou,iel)*thetap
        enddo
      enddo
    enddo
  else
    do iel = 1, ncel
      do isou = 1, 3
        do jsou = 1, 3
          fimp(jsou,isou,iel) = fimp(jsou,isou,iel)                      &
                              + max(-tsimp(jsou,isou,iel),zero)
        enddo
      enddo
    enddo
  endif
endif

!-------------------------------------------------------------------------------
! --->  Mass source terms

if (ncetsm.gt.0) then

!     On calcule les termes Gamma (uinj - u)
!       -Gamma u a la premiere iteration est mis dans
!          TRAV ou TRAVA selon qu'on itere ou non sur navsto
!       Gamma uinj a la premiere iteration est placee dans W1
!       ROVSDT a chaque iteration recoit Gamma
  allocate(gavinj(3,ncelet))
  if (nterup.eq.1) then
    call catsmv(ncetsm, iterns, icetsm, itypsm(1,iu),               &
                cell_f_vol, vela, smacel(:,iu), smacel(:,ipr),      &
                trav, fimp, gavinj)
  else
    call catsmv(ncetsm, iterns, icetsm, itypsm(1,iu),               &
                cell_f_vol, vela, smacel(:,iu), smacel(:,ipr),      &
                trava, fimp, gavinj)
  endif

  ! At the first inner iteration, the explicit part "Gamma u^{in}" is added
  if (iterns.eq.1) then
    ! If source terms are extrapolated, stored in fields
    if(isno2t.gt.0) then
      do iel = 1, ncel
        do isou = 1, 3
          c_st_vel(isou,iel) = c_st_vel(isou,iel) + gavinj(isou,iel)
        enddo
      enddo

    else
      ! If no inner iteration: in trav
      if (nterup.eq.1) then
        do iel = 1,ncel
          do isou = 1, 3
            trav(isou,iel)  = trav(isou,iel) + gavinj(isou,iel)
          enddo
        enddo
      ! Otherwise, in trava
      else
        do iel = 1,ncel
          do isou = 1, 3
            trava(isou,iel) = trava(isou,iel) + gavinj(isou,iel)
          enddo
        enddo
      endif
    endif
  endif

  deallocate(gavinj)

endif

! ---> Right Hand Side initialization

! If source terms are extrapolated in time
if (isno2t.gt.0) then
  thetp1 = 1.d0 + thets
  ! If no inner iteration: trav
  if (nterup.eq.1) then
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel) + thetp1*c_st_vel(isou,iel)
      enddo
    enddo

  else
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel) + trava(isou,iel)       &
                       + thetp1*c_st_vel(isou,iel)
      enddo
    enddo
  endif

! No time extrapolation
else
  ! No inner iteration
  if (nterup.eq.1) then
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel)
      enddo
    enddo
  ! Inner iterations
  else
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel) + trava(isou,iel)
      enddo
    enddo
  endif
endif

! ---> LAGRANGIEN : COUPLAGE RETOUR

!     L'ordre 2 sur les termes issus du lagrangien necessiterait de
!       decomposer TSLAGR(IEL,ISOU) en partie implicite et
!       explicite, comme c'est fait dans cs_user_source_terms.
!     Pour le moment, on n'y touche pas.

if (iilagr.eq.2 .and. ltsdyn.eq.1)  then

  call field_get_val_v_by_name('velocity_st_lagr', lagr_st_vel)

  do iel = 1, ncel
    do isou = 1, 3
      smbr(isou,iel) = smbr(isou,iel) + lagr_st_vel(isou,iel)
    enddo
  enddo
  ! At the second call, fimp is unused
  if(iappel.eq.1) then
    do iel = 1, ncel
      do isou = 1, 3
        fimp(isou,isou,iel) = fimp(isou,isou,iel) + max(-tslagr(iel,itsli),zero)
      enddo
    enddo
  endif

endif

! ---> Electric Arc (Laplace Force)
!      (No 2nd order in time yet)
if (ippmod(ielarc).ge.1) then
  call field_get_val_v_by_name('laplace_force', lapla)
  do iel = 1, ncel
    smbr(1,iel) = smbr(1,iel) + cell_f_vol(iel) * lapla(1,iel)
    smbr(2,iel) = smbr(2,iel) + cell_f_vol(iel) * lapla(2,iel)
    smbr(3,iel) = smbr(3,iel) + cell_f_vol(iel) * lapla(3,iel)
  enddo
endif

! Solver parameters

if (ippmod(icompf).ge.0) then
  ! impose boundary convective flux at some faces (face indicator icvfli)
  icvflb = 1
else
  ! all boundary convective flux with upwind
  icvflb = 0
endif

if (staggered.eq.1) then
  do iel = 1, ncel
    smbr(1,iel) = 0
    smbr(2,iel) = 0
    smbr(3,iel) = 0
  enddo
endif

if (iappel.eq.1) then

  ! Store fimp as the velocity matrix is stored in it in codtiv call
  do iel = 1, ncel
    do isou = 1, 3
      do jsou = 1, 3
        fimpcp(jsou,isou,iel) = fimp(jsou,isou,iel)
      enddo
    enddo
  enddo

  iescap = 0
  if (iespre.ge.0) iescap = 1

  vcopt_loc = vcopt_u

  vcopt_loc%istat  = -1
  vcopt_loc%idifft = -1
  vcopt_loc%iwgrec = 0
  vcopt_loc%blend_st = 0 ! Warning, may be overwritten if a field

  p_k_value => vcopt_loc
  c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

  ! Warning: in case of convergence estimators, eswork give the estimator
  ! of the predicted velocity
  call cs_equation_iterative_solve_vector                     &
   ( idtvar , iterns ,                                        &
     ivarfl(iu)      , c_null_char     ,                      &
     ivisse , iescap , c_k_value       ,                      &
     vela   , velk   ,                                        &
     coefav , coefbv , cofafv , cofbfv ,                      &
     imasfl , bmasfl , viscfi ,                               &
     viscbi , viscf  , viscb  , secvif ,                      &
     secvib , rvoid  , rvoid  , rvoid  ,                      &
     icvflb , icvfli ,                                        &
     fimp   , smbr   , vel    , eswork )


  ! Compute kinetic energy balance for compressible algorithme
  ! See H. Amino thesis
  call cs_thermal_model_kinetic_st_prepare &
    ( imasfl, bmasfl, vela, vel)

  ! Store inverse of the velocity matrix for the correction step
  !  if needed (otherwise vitenp is used in cs_pressure_correction)
  if (rcfact.eq.1) then

    do iel = 1, ncel

      tensor(1) = fimp(1,1,iel)/crom(iel)
      tensor(2) = fimp(2,2,iel)/crom(iel)
      tensor(3) = fimp(3,3,iel)/crom(iel)
      tensor(4) = fimp(1,2,iel)/crom(iel)
      tensor(5) = fimp(2,3,iel)/crom(iel)
      tensor(6) = fimp(1,3,iel)/crom(iel)

      call symmetric_matrix_inverse(tensor, da_uu(:, iel))

      do ii = 1, 6
        da_uu(ii,iel) = cell_f_vol(iel)*da_uu(ii,iel)
      enddo

    enddo

    call syntis(da_uu)

  endif

  ! Velocity-pression coupling: compute the vector T, stored in tpucou,
  ! cs_equation_iterative_solve_vector is called, only one sweep is done,
  ! and tpucou is initialized by 0, so that the advection/diffusion added
  ! by cs_balance_vector is 0.
  !  nswrsp = -1 indicated that only one sweep is required and inc=0
  !  for boundary contitions on the weight matrix.
  if (ipucou.eq.1) then

    ! Allocate temporary arrays for the velocity-pressure resolution
    allocate(vect(3,ncelet))

    nswrsp = -1
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = cell_f_vol(iel)
      enddo
    enddo
    do iel = 1, ncelet
      do isou = 1, 3
        vect(isou,iel) = 0.d0
      enddo
    enddo
    iescap = 0

    ! We do not take into account transpose of grad
    ivisep = 0

    vcopt_loc%nswrsm = nswrsp

    p_k_value => vcopt_loc
    c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

    call cs_equation_iterative_solve_vector                 &
 ( idtvar , iterns ,                                        &
   ivarfl(iu)      , c_null_char     ,                      &
   ivisep , iescap , c_k_value       ,                      &
   vect   , vect   ,                                        &
   coefav , coefbv , cofafv , cofbfv ,                      &
   imasfl , bmasfl , viscfi ,                               &
   viscbi , viscf  , viscb  , secvif ,                      &
   secvib , rvoid  , rvoid  , rvoid  ,                      &
   icvflb , ivoid  ,                                        &
   fimpcp , smbr   , vect   , rvoid  )

    do iel = 1, ncelet
      rom = crom(iel)
      do isou = 1, 3
        tpucou(isou,iel) = rom*vect(isou,iel)
      enddo
      do isou = 4, 6
        tpucou(isou,iel) = 0.d0
      enddo
    enddo

    ! Free memory
    deallocate(vect)

  endif

  ! ---> The estimator on the predicted velocity is summed up over the components
  if (iespre.ge.0) then
    call field_get_val_s(iespre, c_estim)
    do iel = 1, ncel
      do isou = 1, 3
        c_estim(iel) =  c_estim(iel) + eswork(isou,iel)
      enddo
    enddo
  endif


! ---> End of the construction of the total estimator:
!       RHS residual of (U^{n+1}, P^{n+1}) + rho*volume*(U^{n+1} - U^n)/dt
else if (iappel.eq.2) then

  inc = 1
  ! Pas de relaxation en stationnaire
  idtva0 = 0
  imasac = 0

  vcopt_loc = vcopt_u

  vcopt_loc%istat  = -1
  vcopt_loc%idifft = -1
  vcopt_loc%iswdyn = -1
  vcopt_loc%nswrsm = -1
  vcopt_loc%iwgrec = 0
  vcopt_loc%blend_st = 0 ! Warning, may be overwritten if a field
  vcopt_loc%epsilo = -1
  vcopt_loc%epsrsm = -1

  p_k_value => vcopt_loc
  c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

  call cs_balance_vector &
 ( idtva0 , ivarfl(iu)      , imasac , inc    , ivisse ,                   &
   c_k_value                , vel    , vel    , coefav , coefbv , cofafv , &
   cofbfv , imasfl , bmasfl , viscf  , viscb  , secvif ,                   &
   secvib , rvoid  , rvoid  , rvoid  ,                                     &
   icvflb , icvfli ,c_null_ptr  ,c_null_ptr , smbr   )

  call field_get_val_s(iestot, c_estim)
  do iel = 1, ncel
    do isou = 1, 3
      c_estim(iel) = c_estim(iel) + (smbr(isou,iel)/volume(iel))**2
    enddo
  enddo
endif

! ---> Finilaze estimators + Printings

call field_get_id_try("predicted_velocity", f_id)
if (f_id.ge.0) then
  call field_get_val_v(f_id, cpro_pred_vel)
  do iel = 1, ncel
    do isou = 1, 3
      cpro_pred_vel(isou, iel) = vel(isou, iel)
    enddo
  enddo
endif

if (iappel.eq.1) then

  ! ---> Estimator on the predicted velocity:
  !      square root (norm) or square root of the sum times the volume (L2 norm)
  if (iespre.ge.0) then
    call field_get_val_s(iespre, c_estim)
    do iel = 1, ncel
      c_estim(iel) = sqrt(c_estim(iel)*volume(iel))
    enddo
  endif

  ! ---> Norm printings
  if (vcopt_u%iwarni.ge.2) then
    rnorm = -1.d0
    do iel = 1, ncel
      vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
      rnorm = max(rnorm,vitnor)
    enddo

    if (irangp.ge.0) call parmax (rnorm)

    write(nfecra,1100) rnorm

    do iel = 1, ncel
      vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
      rnorm = min(rnorm,vitnor)
    enddo

    if (irangp.ge.0) call parmin (rnorm)

    write(nfecra,1200) rnorm

  endif

! ---> Estimator on the whole Navier-Stokes:
!      square root (norm) or square root of the sum times the volume (L2 norm)
else if (iappel.eq.2) then

  call field_get_val_s(iestot, c_estim)
  if (iestot.ge.0) then
    do iel = 1, ncel
      c_estim(iel) = sqrt(c_estim(iel)*volume(iel))
    enddo
  endif

endif

! Free memory
!------------
deallocate(smbr)
deallocate(fimp)
deallocate(fimpcp)
deallocate(w1)
if (allocated(loctsexp)) deallocate(loctsexp)
if (allocated(loctsimp)) deallocate(loctsimp)
if (allocated(viscce)) deallocate(viscce)
if (allocated(divt)) deallocate(divt)
if (allocated(cproa_rho_tc)) deallocate(cproa_rho_tc)
if (allocated(stf)) deallocate(stf)

!--------
! Formats
!--------

 1100 format(/,                                                   &
 1X,'Maximum velocity after prediction ',E12.4)

 1200 format(/,                                                   &
 1X,'Minimum velocity after prediction ',E12.4)

!----
! End
!----

return

end subroutine

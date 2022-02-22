!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file reseps.f90
!>
!> \brief This subroutine performs the solving of epsilon in
!>        \f$ R_{ij} - \varepsilon \f$ RANS turbulence model.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for each variable
!> \param[in]     dt            time step (per cell)
!> \param[in]     gradv         work array for the term grad
!>                               of velocity only for iturb=31
!> \param[in]     produc        work array for production (without
!>                               rho volume) only for iturb=30
!> \param[in]     gradro        work array for \f$ \grad{rom} \f$
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate
!> \param[in]     viscf         visc*surface/dist at internal faces
!> \param[in]     viscb         visc*surface/dist at edge faces
!> \param[in]     smbr          working array
!> \param[in]     rovsdt        working array
!_______________________________________________________________________________

subroutine reseps &
 ( nvar   , ncesmp ,                                              &
   icetsm , itypsm ,                                              &
   dt     ,                                                       &
   gradv  , produc , gradro ,                                     &
   smacel ,                                                       &
   viscf  , viscb  ,                                              &
   smbr   , rovsdt )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use optcal
use cstnum
use cstphy
use parall
use period
use lagran
use mesh
use field
use field_operator
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar
integer          ncesmp

integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision produc(6,ncelet), gradv(3, 3, ncelet)
double precision gradro(3,ncelet)
double precision smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbr(ncelet), rovsdt(ncelet)

! Local variables

integer          init
integer          ivar
integer          iel
integer          iflmas, iflmab
integer          iwarnp
integer          imvisp
integer          iescap
integer          st_prv_id
integer          imucpp
integer          icvflb
integer          ivoid(1)
integer          key_t_ext_id
integer          iroext
double precision alpha3
double precision trprod , trrij
double precision tseps , kseps , ceps2
double precision tuexpe, thets , thetv , thetap, thetp1
double precision prdeps, xttdrb, xttke , xttkmg
double precision normp
double precision sigmae

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, cprod
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer ::  crom, cromo
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:,:), pointer :: visten
double precision, dimension(:), pointer :: cvar_ep, cvara_ep, cvar_al
double precision, dimension(:,:), pointer :: cvara_rij, lagr_st_rij
double precision, dimension(:), pointer :: viscl, visct, c_st_prv

character(len=80) :: label

type(var_cal_opt) :: vcopt
type(var_cal_opt), target :: vcopt_loc
type(var_cal_opt), pointer :: p_k_value
type(c_ptr) :: c_k_value

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_key_int(icrom, key_t_ext_id, iroext)

ivar = iep

! Allocate work arrays
allocate(w1(ncelet))
allocate(cprod(ncelet))
allocate(dpvar(ncelet))
allocate(viscce(6,ncelet))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1000) label
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

call field_get_val_s(ivarfl(iep), cvar_ep)
call field_get_val_prev_s(ivarfl(iep), cvara_ep)
call field_get_key_double(ivarfl(iep), ksigmas, sigmae)
if (iturb.eq.32) call field_get_val_s(ivarfl(ial), cvar_al)
call field_get_val_prev_v(ivarfl(irij), cvara_rij)
call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Constant Ce2, which worths Ce2 for iturb=30 and CSSGE2 for itrub=31
if (iturb.eq.30) then
  ceps2 = ce2
elseif (iturb.eq.31) then
  ceps2 = cssge2
else
  ceps2 = cebme2
endif

! S as Source, V as Variable
thets  = thetst
thetv  = vcopt%thetav

call field_get_key_int(ivarfl(ivar), kstprv, st_prv_id)
if (st_prv_id.ge.0) then
  call field_get_val_s(st_prv_id, c_st_prv)
else
  c_st_prv=> null()
endif

if (st_prv_id.ge.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 2. User source terms
!===============================================================================

call user_source_terms(ivarfl(ivar), smbr, rovsdt)

!     If we extrapolate the source terms
if (st_prv_id.ge.0) then
  do iel = 1, ncel
    !       Save for exchange
    tuexpe = c_st_prv(iel)
    !       For the continuation and the next time step
    c_st_prv(iel) = smbr(iel)
    !       Second member of previous time step
    !       We suppose -rovsdt > 0: we implicit
    !          the user source term (the rest)
    smbr(iel) = rovsdt(iel)*cvara_ep(iel) - thets*tuexpe
    !       Diagonal
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel)   = rovsdt(iel)*cvara_ep(iel) + smbr(iel)
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

!===============================================================================
! 3. Lagrangian source terms
!===============================================================================

!     Second order is not taken into account
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  call field_get_val_v_by_name('rij_st_lagr', lagr_st_rij)

  do iel = 1, ncel
    ! Source terms with eps
    tseps = -0.5d0 * (  lagr_st_rij(1,iel)                        &
                      + lagr_st_rij(2,iel)                        &
                      + lagr_st_rij(3,iel))
    ! quotient k/eps
    kseps = 0.5d0 * (  cvara_rij(1,iel)                           &
                     + cvara_rij(2,iel)                           &
                     + cvara_rij(3,iel))                          &
                   / cvara_ep(iel)
    smbr(iel)   = smbr(iel) + ce4 *tseps *cvara_ep(iel) /kseps
    rovsdt(iel) = rovsdt(iel) + max( (-ce4*tseps/kseps) , zero)
  enddo

endif

!===============================================================================
! 4. Mass source term
!===============================================================================

if (ncesmp.gt.0) then

  ! We increment smbr with -Gamma.var_prev. and rovsdt with Gamma
  call catsma(ncesmp, 1, icetsm, itypsm(:,ivar),                    &
              cell_f_vol, cvara_ep, smacel(:,ivar), smacel(:,ipr),  &
              smbr, rovsdt, w1)

  ! If we extrapolate the source terms, we put Gamma Pinj in c_st_prv
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      c_st_prv(iel) = c_st_prv(iel) + w1(iel)
    enddo
  ! Otherwise we put it directly in smbr
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w1(iel)
    enddo
  endif

endif

!===============================================================================
! 5. Unsteady term
!===============================================================================

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                          &
              + vcopt%istat*(crom(iel)/dt(iel))*cell_f_vol(iel)
enddo

!===============================================================================
! 6. Production (rho * Ce1 * epsilon / k * P)
!    Dissipation (rho*Ce2.epsilon/k*epsilon)
!===============================================================================

if (st_prv_id.ge.0) then
  thetap = thetv
else
  thetap = 1.d0
endif

! ---> Calculation the production trace, depending we are in standard
!     Rij or in SSG (use of produc or grdvit)
if (iturb.eq.30) then
  do iel = 1, ncel
    cprod(iel) = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
  enddo
else
  do iel = 1, ncel
    cprod(iel) = -( cvara_rij(1,iel)*gradv(1, 1, iel) +               &
                    cvara_rij(4,iel)*gradv(2, 1, iel) +               &
                    cvara_rij(6,iel)*gradv(3, 1, iel) +               &
                    cvara_rij(4,iel)*gradv(1, 2, iel) +               &
                    cvara_rij(2,iel)*gradv(2, 2, iel) +               &
                    cvara_rij(5,iel)*gradv(3, 2, iel) +               &
                    cvara_rij(6,iel)*gradv(1, 3, iel) +               &
                    cvara_rij(5,iel)*gradv(2, 3, iel) +               &
                    cvara_rij(3,iel)*gradv(3, 3, iel) )
  enddo
endif


! EBRSM
if (iturb.eq.32) then

  do iel = 1, ncel
    ! Half-traces
    trprod = cprod(iel)
    trrij  = 0.5d0 * (cvara_rij(1,iel) + cvara_rij(2,iel) + cvara_rij(3,iel))
    ! Calculation of the Durbin time scale
    xttke  = trrij/cvara_ep(iel)
    xttkmg = xct*sqrt(viscl(iel)/crom(iel)/cvara_ep(iel))
    xttdrb = max(xttke,xttkmg)

    prdeps = trprod/cvara_ep(iel)
    alpha3 = cvar_al(iel)**3

    ! Production (explicit)
    ! Compute of C_eps_1'
    w1(iel) = cromo(iel)*cell_f_vol(iel)*                                 &
              ce1*(1.d0+xa1*(1.d0-alpha3)*prdeps)*trprod/xttdrb


    ! Dissipation (implicit)
    smbr(iel) = smbr(iel) - crom(iel)*cell_f_vol(iel)*                    &
                             ceps2*cvara_ep(iel)/xttdrb

    rovsdt(iel) = rovsdt(iel)                                         &
                + ceps2*crom(iel)*cell_f_vol(iel)*thetap/xttdrb
  enddo

! SSG and LRR
else

  do iel = 1, ncel
    ! Half-traces
    trprod = cprod(iel)
    trrij  = 0.5d0 * (cvara_rij(1,iel) + cvara_rij(2,iel) + cvara_rij(3,iel))
    xttke  = trrij/cvara_ep(iel)
    ! Production (explicit, a part might be implicit)
    rovsdt(iel) = rovsdt(iel)                                        &
                + max(- cromo(iel)*cell_f_vol(iel)*ce1*trprod/trrij, 0.d0)
    w1(iel) = cromo(iel)*cell_f_vol(iel)*ce1/xttke*trprod

    ! Dissipation (implicit)
    smbr(iel) = smbr(iel)                                            &
              - crom(iel)*cell_f_vol(iel)*ceps2*cvara_ep(iel)**2/trrij
    rovsdt(iel) = rovsdt(iel)                                        &
                + ceps2*crom(iel)*cell_f_vol(iel)/xttke*thetap
  enddo

endif

! Extrapolation of source terms (2nd order in time)
if (st_prv_id.ge.0) then
  do iel = 1, ncel
    c_st_prv(iel) = c_st_prv(iel) + w1(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w1(iel)
  enddo
endif

!===============================================================================
! 7. Buoyancy term
!===============================================================================

!FIXME use beta ... WARNING
if (igrari.eq.1) then

  ! Allocate a work array
  allocate(w7(ncelet))

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  ! Note that w7 is positive.
  if (irijco.eq.1) then
    call rijtheps(gradro, w7)
  else
    call rijthe(ivar, gradro, w7)
  endif

  ! Extrapolation of source terms (2nd order in time)
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      c_st_prv(iel) = c_st_prv(iel) + w7(iel) * cell_f_vol(iel)
    enddo
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w7(iel) * cell_f_vol(iel)
    enddo
  endif

  ! Free memory
  deallocate(w7)

endif

!===============================================================================
! 8. Diffusion term (Daly Harlow: generalized gradient hypothesis method)
!===============================================================================

! Symmetric tensor diffusivity (GGDH)
if (iand(vcopt%idften, ANISOTROPIC_DIFFUSION).ne.0) then

  call field_get_val_v(ivsten, visten)

  do iel = 1, ncel
    viscce(1,iel) = visten(1,iel)/sigmae + viscl(iel)
    viscce(2,iel) = visten(2,iel)/sigmae + viscl(iel)
    viscce(3,iel) = visten(3,iel)/sigmae + viscl(iel)
    viscce(4,iel) = visten(4,iel)/sigmae
    viscce(5,iel) = visten(5,iel)/sigmae
    viscce(6,iel) = visten(6,iel)/sigmae
  enddo

  iwarnp = vcopt%iwarni

  call vitens &
 ( viscce , iwarnp ,             &
   weighf , weighb ,             &
   viscf  , viscb  )

! Scalar diffusivity
else

  do iel = 1, ncel
    w1(iel) = viscl(iel) + vcopt%idifft*visct(iel)/sigmae
  enddo

  imvisp = vcopt%imvisf

  call viscfa                    &
 ( imvisp ,                      &
   w1     ,                      &
   viscf  , viscb  )

endif

!===============================================================================
! 9. Solving
!===============================================================================

if (st_prv_id.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*c_st_prv(iel)
  enddo
endif

iescap = 0
imucpp = 0
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0
init   = 1

vcopt_loc = vcopt

vcopt_loc%istat  = -1
vcopt_loc%icoupl = -1
vcopt_loc%idifft = -1
vcopt_loc%iwgrec = 0 ! Warning, may be overwritten if a field
vcopt_loc%thetav = thetv
vcopt_loc%blend_st = 0 ! Warning, may be overwritten if a field

p_k_value => vcopt_loc
c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

call cs_equation_iterative_solve_scalar          &
 ( idtvar , init   ,                             &
   ivarfl(ivar)     , c_null_char ,              &
   iescap , imucpp , normp  , c_k_value       ,  &
   cvara_ep        , cvara_ep        ,           &
   coefap , coefbp , cofafp , cofbfp ,           &
   imasfl , bmasfl ,                             &
   viscf  , viscb  , viscf  , viscb  ,           &
   viscce , weighf , weighb ,                    &
   icvflb , ivoid  ,                             &
   rovsdt , smbr   , cvar_ep, dpvar  ,           &
   rvoid  , rvoid  )

! Free memory
deallocate(w1)
deallocate(cprod)
deallocate(viscce)
deallocate(weighf, weighb)

!--------
! Formats
!--------

 1000 format(/,'           Solving variable ', a8          ,/)

!----
! End
!----

return

end subroutine

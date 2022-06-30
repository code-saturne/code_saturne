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

!> \file resssg2.f90
!>
!> \brief This subroutine prepares the solving of the coupled Reynolds stress
!> components in \f$ R_{ij} - \varepsilon \f$ RANS (SSG) turbulence model.
!>
!> \remark
!> - cvar_var(1,*) for \f$ R_{11} \f$
!> - cvar_var(2,*) for \f$ R_{22} \f$
!> - cvar_var(3,*) for \f$ R_{33} \f$
!> - cvar_var(4,*) for \f$ R_{12} \f$
!> - cvar_var(5,*) for \f$ R_{23} \f$
!> - cvar_var(6,*) for \f$ R_{13} \f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ivar          variable number
!> \param[in]     gradv         work array for the velocity grad term
!>                                 only for iturb=31
!> \param[in]     produc        work array for production
!> \param[in]     gradro        work array for grad rom
!>                              (without rho volume) only for iturb=30
!> \param[in]     viscf         visc*surface/dist at internal faces
!> \param[in]     viscb         visc*surface/dist at edge faces
!> \param[out]    viscce        Daly Harlow diffusion term
!> \param[in]     smbr          working array
!> \param[in]     rovsdt        working array
!> \param[out]    weighf        working array
!> \param[out]    weighb        working array
!_______________________________________________________________________________

subroutine resssg2 &
 ( ivar   ,                                                       &
   gradv  , produc , gradro ,                                     &
   viscf  , viscb  , viscce ,                                     &
   smbr   , rovsdt ,                                              &
   weighf , weighb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use lagran
use mesh
use field
use field_operator
use cs_f_interfaces
use rotation
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ivar

double precision gradv(3, 3, ncelet)
double precision produc(6, ncelet)
double precision gradro(3, ncelet)
double precision viscf(nfac), viscb(nfabor), viscce(6, ncelet)
double precision smbr(6, ncelet)
double precision rovsdt(6, 6, ncelet)
double precision weighf(2, nfac), weighb(nfabor)

! Local variables

integer          iel, isou, jsou
integer          ii    , jj    , kk    , iii   , jjj
integer          iwarnp
integer          imvisp
integer          st_prv_id
integer          iprev , inc, iccocg, ll
integer          t2v(3,3)
integer          iv2t(6), jv2t(6)
integer          key_t_ext_id, f_id, rot_id
integer          iroext

double precision trprod, trrij
double precision deltij(6), dij(3,3)
double precision thets , thetv
double precision aiksjk, aikrjk, aii ,aklskl, aikakj
double precision xaniso(3,3), xstrai(3,3), xrotac(3,3), xprod(3,3), matrot(3,3)
double precision sym_strain(6)
double precision xrij(3,3), xnal(3), xnoral, xnnd
double precision d1s2, d1s3, d2s3
double precision alpha3
double precision pij, phiij1, phiij2, epsij
double precision phiijw, epsijw, epsijw_imp
double precision ccorio
double precision rctse
double precision eigen_max
double precision eigen_vals(3)
double precision turb_schmidt
double precision gradchk, gradro_impl, const
double precision kseps
double precision matrn(6), oo_matrn(6)
double precision ceps_impl, cphi3impl, cphiw_impl
double precision impl_drsm(6,6)
double precision implmat2add(3,3)
double precision impl_lin_cst, impl_id_cst
double precision gkks3
double precision grav(3)
double precision, dimension(3,3) :: cvara_r

character(len=80) :: label
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2
double precision, allocatable, dimension(:,:), target :: buoyancy
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:,:), pointer :: visten
double precision, dimension(:), pointer :: cvara_ep, cvar_al
double precision, dimension(:,:), pointer :: cvar_var, cvara_var
double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:,:), pointer :: c_st_prv
double precision, dimension(:,:), pointer :: cpro_buoyancy

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_key_int(icrom, key_t_ext_id, iroext)

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet))

! Generating the tensor to vector (t2v) and vector to tensor (v2t) mask arrays

! a) t2v
t2v(1,1) = 1; t2v(1,2) = 4; t2v(1,3) = 6;
t2v(2,1) = 4; t2v(2,2) = 2; t2v(2,3) = 5;
t2v(3,1) = 6; t2v(3,2) = 5; t2v(3,3) = 3;

! b) i index of v2t
iv2t(1) = 1; iv2t(2) = 2; iv2t(3) = 3;
iv2t(4) = 1; iv2t(5) = 2; iv2t(6) = 1;

! c) j index of v2t
jv2t(1) = 1; jv2t(2) = 2; jv2t(3) = 3;
jv2t(4) = 2; jv2t(5) = 3; jv2t(6) = 3;

! d) kronecker symbol
dij(:,:) = 0.0d0;
dij(1,1) = 1.0d0; dij(2,2) = 1.0d0; dij(3,3) = 1.0d0;

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1000) label
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

call field_get_val_prev_s(ivarfl(iep), cvara_ep)
if (iturb.ne.31) call field_get_val_s(ivarfl(ial), cvar_al)

call field_get_val_v(ivarfl(ivar), cvar_var)
call field_get_val_prev_v(ivarfl(ivar), cvara_var)

d1s2 = 1.d0/2.d0
d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

do isou = 1, 3
  deltij(isou) = 1.0d0
enddo
do isou = 4, 6
  deltij(isou) = 0.0d0
enddo

!     S as Source, V as Variable
thets  = thetst
thetv  = vcopt%thetav

call field_get_key_int(ivarfl(ivar), kstprv, st_prv_id)
if (st_prv_id .ge. 0) then
  call field_get_val_v(st_prv_id, c_st_prv)
else
  c_st_prv=> null()
endif

if (st_prv_id.ge.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif

! Coefficient of the "Coriolis-type" term
if (icorio.eq.1) then
  ! Relative velocity formulation
  ccorio = 2.d0
elseif (iturbo.eq.1) then
  ! Mixed relative/absolute velocity formulation
  ccorio = 1.d0
else
  ccorio = 0.d0
endif

!===============================================================================
! Production, Pressure-Strain correlation, dissipation, Coriolis
!===============================================================================

! Source term
!  -rho*epsilon*( Cs1*aij + Cs2*(aikajk -1/3*aijaij*deltaij))
!  -Cr1*P*aij + Cr2*rho*k*sij - Cr3*rho*k*sij*sqrt(aijaij)
!  +Cr4*rho*k(aik*sjk+ajk*sik-2/3*akl*skl*deltaij)
!  +Cr5*rho*k*(aik*rjk + ajk*rik)
!  -2/3*epsilon*deltaij

! EBRSM
if (iturb.eq.32) then

  allocate(grad(3,ncelet))

  ! Compute the gradient of Alpha
  iprev  = 1
  inc    = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(ial), iprev, 0, inc, iccocg, grad)

endif

do iel = 1, ncel

  ! EBRSM
  if (iturb.eq.32) then
    ! Compute the magnitude of the Alpha gradient
    xnoral = ( grad(1,iel)*grad(1,iel)          &
           +   grad(2,iel)*grad(2,iel)          &
           +   grad(3,iel)*grad(3,iel) )
    xnoral = sqrt(xnoral)
   ! Compute the unitary vector of Alpha
    if (xnoral.le.epzero) then
      xnal(1) = 0.d0
      xnal(2) = 0.d0
      xnal(3) = 0.d0
    else
      xnal(1) = grad(1,iel)/xnoral
      xnal(2) = grad(2,iel)/xnoral
      xnal(3) = grad(3,iel)/xnoral
    endif
  endif

  ! Pij
  xprod(1,1) = produc(1, iel)
  xprod(1,2) = produc(4, iel)
  xprod(1,3) = produc(6, iel)
  xprod(2,2) = produc(2, iel)
  xprod(2,3) = produc(5, iel)
  xprod(3,3) = produc(3, iel)

  ! Rotating frame of reference => "Coriolis production" term

  rot_id = icorio
  if (iturbo.eq.1) rot_id = irotce(iel)

  if (rot_id .ge. 1) then

    call coriolis_t(rot_id, 1.d0, matrot)
    cvara_r(1,1) = cvara_var(1,iel)
    cvara_r(2,2) = cvara_var(2,iel)
    cvara_r(3,3) = cvara_var(3,iel)
    cvara_r(1,2) = cvara_var(4,iel)
    cvara_r(2,3) = cvara_var(5,iel)
    cvara_r(1,3) = cvara_var(6,iel)
    cvara_r(2,1) = cvara_var(4,iel)
    cvara_r(3,2) = cvara_var(5,iel)
    cvara_r(3,1) = cvara_var(6,iel)

    do ii = 1, 3
      do jj = ii, 3
        do kk = 1, 3
          xprod(ii,jj) = xprod(ii,jj)                             &
                          - ccorio*( matrot(ii,kk)*cvara_r(jj,kk) &
                          + matrot(jj,kk)*cvara_r(ii,kk) )
        enddo
      enddo
    enddo
  endif

  xprod(2,1) = xprod(1,2)
  xprod(3,1) = xprod(1,3)
  xprod(3,2) = xprod(2,3)

  trprod = d1s2 * (xprod(1,1) + xprod(2,2) + xprod(3,3) )
  trrij  = d1s2 * (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel))

  !-----> aII = aijaij
  aii    = 0.d0
  aklskl = 0.d0
  aiksjk = 0.d0
  aikrjk = 0.d0
  aikakj = 0.d0

  ! aij
  xaniso(1,1) = cvara_var(1,iel)/trrij - d2s3
  xaniso(2,2) = cvara_var(2,iel)/trrij - d2s3
  xaniso(3,3) = cvara_var(3,iel)/trrij - d2s3
  xaniso(1,2) = cvara_var(4,iel)/trrij
  xaniso(1,3) = cvara_var(6,iel)/trrij
  xaniso(2,3) = cvara_var(5,iel)/trrij
  xaniso(2,1) = xaniso(1,2)
  xaniso(3,1) = xaniso(1,3)
  xaniso(3,2) = xaniso(2,3)

  ! Sij
  xstrai(1,1) = gradv(1, 1, iel)
  xstrai(1,2) = d1s2*(gradv(2, 1, iel)+gradv(1, 2, iel))
  xstrai(1,3) = d1s2*(gradv(3, 1, iel)+gradv(1, 3, iel))
  xstrai(2,1) = xstrai(1,2)
  xstrai(2,2) = gradv(2, 2, iel)
  xstrai(2,3) = d1s2*(gradv(3, 2, iel)+gradv(2, 3, iel))
  xstrai(3,1) = xstrai(1,3)
  xstrai(3,2) = xstrai(2,3)
  xstrai(3,3) = gradv(3, 3, iel)

  ! omegaij
  xrotac(1,1) = 0.d0
  xrotac(1,2) = d1s2*(gradv(2, 1, iel)-gradv(1, 2, iel))
  xrotac(1,3) = d1s2*(gradv(3, 1, iel)-gradv(1, 3, iel))
  xrotac(2,1) = -xrotac(1,2)
  xrotac(2,2) = 0.d0
  xrotac(2,3) = d1s2*(gradv(3, 2, iel)-gradv(2, 3, iel))
  xrotac(3,1) = -xrotac(1,3)
  xrotac(3,2) = -xrotac(2,3)
  xrotac(3,3) = 0.d0

  do ii=1,3
    do jj = 1,3
      ! aii = aij.aij
      aii    = aii+xaniso(ii,jj)*xaniso(ii,jj)
      ! aklskl = aij.Sij
      aklskl = aklskl + xaniso(ii,jj)*xstrai(ii,jj)
    enddo
  enddo

  if (irijco .ne. 0) then

    ! Initalize implicit matrices at 0
    do isou = 1, 6
      do jsou = 1, 6
        impl_drsm(isou, jsou) = 0.0d0
      end do
    end do
    do isou = 1, 3
      do jsou = 1, 3
        implmat2add(isou, jsou) = 0.0d0
      end do
    end do

    impl_lin_cst = 0.0d0
    impl_id_cst  = 0.0d0

    ! Computation of implicit components
    ! -----------------------------------

    sym_strain(1) = xstrai(1,1)
    sym_strain(2) = xstrai(2,2)
    sym_strain(3) = xstrai(3,3)
    sym_strain(4) = xstrai(1,2)
    sym_strain(5) = xstrai(2,3)
    sym_strain(6) = xstrai(1,3)

    ! Global variables needed for SSG and EBRSM
    ! -------------
    ! Computing the inverse matrix of R^n
    ! Scaling by tr(R) in order to dodge inversion errors
    do isou = 1, 6
      matrn(isou) = cvara_var(isou,iel)/trrij
      oo_matrn(isou) = 0.0d0
    end do

    ! Inversing the matrix
    call symmetric_matrix_inverse(matrn, oo_matrn)
    do isou = 1, 6
      oo_matrn(isou) = oo_matrn(isou)/trrij
    end do

    ! Computing the maximal eigenvalue (in terms of norm!) of S
    call calc_symtens_eigvals(sym_strain, eigen_vals)
    eigen_max = maxval(abs(eigen_vals))

    ! Constant for the dissipation
    ceps_impl = d1s3 * cvara_ep(iel)

    if (iturb .eq. 31 ) then ! SSG - Epsilon

      ! Identity constant for phi3
      cphi3impl = abs(cssgr2 - cssgr3*sqrt(aii))

      ! Identity constant
      impl_id_cst = - d2s3*cssgr1*min(trprod,0.0d0)         & ! Phi1
                    - d1s3*cssgs2*cvara_ep(iel)*aii         & ! Phi2
                    + cphi3impl * trrij * eigen_max         & ! Phi3
                    + 2.0d0*d2s3*cssgr4*trrij*eigen_max     & ! Phi4
                    + d2s3*trrij*cssgr4*max(aklskl,0.0d0)     ! Phi4

      ! Linear constant
      impl_lin_cst = eigen_max *     ( &
                     1.0d0             & ! Production
                   + cssgr4            & ! Phi 4 linear part
                   + cssgr5          )   ! Phi 5 linear part

      do jsou = 1, 3
        do isou = 1 ,3
          iii = t2v(isou,jsou)
          implmat2add(isou,jsou) = xrotac(isou,jsou)              &
                                 + impl_lin_cst*deltij(iii)       &
                                 + impl_id_cst*d1s2*oo_matrn(iii) &
                                 + ceps_impl*oo_matrn(iii)
        end do
      end do

      impl_drsm(:,:) = 0.0d0
      call reduce_symprod33_to_66(implmat2add, impl_drsm)

    else ! EBRSM

      alpha3 = cvar_al(iel)**3

      ! Phi3 constant
      cphi3impl = abs(cebmr2 - cebmr3*sqrt(aii))

      ! PhiWall + epsilon_wall constants for EBRSM
      cphiw_impl = 6.0d0*(1.0d0-alpha3)*cvara_ep(iel)/trrij

      ! -------------
      ! The implicit components of Phi (pressure-velocity fluctuations)
      ! are split into the linear part (A*R) and Id part (A*Id).
      ! -------------

      ! Identity constant
      impl_id_cst = alpha3 * (                             &
                    - d2s3*cebmr1*min(trprod,0.0d0)        & ! Phi1
                    + cphi3impl * trrij * eigen_max        & ! Phi3
                    + 2.0d0*d2s3*cebmr4*trrij*eigen_max    & ! Phi4
                    + d2s3*trrij*cebmr4*max(aklskl,0.d0) )   ! Phi4

      ! Linear constant
      impl_lin_cst = eigen_max * (                  &
                     1.0d0                          & ! Production
                     + cebmr4 * alpha3              & ! Phi4 Linear part
                     + cebmr5 * alpha3            ) & ! Phi5 Linear part
                     + cphiw_impl                     ! Epsilon + Phi wall

      do jsou = 1, 3
        do isou = 1, 3
          iii = t2v(isou,jsou)
          implmat2add(isou,jsou) = xrotac(isou,jsou)                  &
                                 + impl_lin_cst*deltij(iii)           &
                                 + impl_id_cst*d1s2*oo_matrn(iii)     &
                                 + alpha3*ceps_impl*oo_matrn(iii)
        end do
      end do

      impl_drsm(:,:) = 0.0d0
      call reduce_symprod33_to_66(implmat2add, impl_drsm)

    end if   ! EBRSM

  endif   ! (irijco .ne. 0)

  ! Rotating frame of reference => "absolute" vorticity
  if (icorio.eq.1) then
    do ii = 1, 3
      do jj = 1, 3
        xrotac(ii,jj) = xrotac(ii,jj) + matrot(ii,jj)
      enddo
    enddo
  endif

  do isou = 1, 6
    iii = iv2t(isou)
    jjj = jv2t(isou)
    aiksjk = 0
    aikrjk = 0
    aikakj = 0
    do kk = 1,3
      ! aiksjk = aik.Sjk+ajk.Sik
      aiksjk =   aiksjk + xaniso(iii,kk)*xstrai(jjj,kk)              &
               + xaniso(jjj,kk)*xstrai(iii,kk)
      ! aikrjk = aik.Omega_jk + ajk.omega_ik
      aikrjk =   aikrjk + xaniso(iii,kk)*xrotac(jjj,kk)              &
               + xaniso(jjj,kk)*xrotac(iii,kk)
      ! aikakj = aik*akj
      aikakj = aikakj + xaniso(iii,kk)*xaniso(kk,jjj)
    enddo

    ! If we extrapolate the source terms (rarely),
    ! we put all in the previous ST.
    ! We do not implicit the term with Cs1*aij neither the term with Cr1*P*aij.
    ! Otherwise, we put all in smbr and we can implicit Cs1*aij
    ! and Cr1*P*aij. Here we store the second member and the implicit term
    ! in W1 and W2, to avoid the test(ST_PRV_ID.GE.0)
    ! in the ncel loop
    ! In the term with W1, which is dedicated to be extrapolated, we use
    ! cromo.
    ! The implicitation of the two terms can also be done in the case of
    ! extrapolation, by isolating those two terms and by putting it in
    ! the RHS but not in the prev. ST and by using ipcrom .... to be modified
    ! if needed.

    if (iturb.eq.31) then

      ! Explicit terms
      pij = xprod(iii,jjj)
      phiij1 = -cvara_ep(iel)* &
         (cssgs1*xaniso(iii,jjj)+cssgs2*(aikakj-d1s3*deltij(isou)*aii))
      phiij2 = - cssgr1*trprod*xaniso(iii,jjj)                             &
             +   trrij*xstrai(iii,jjj)*(cssgr2-cssgr3*sqrt(aii))           &
             +   cssgr4*trrij*(aiksjk-d2s3*deltij(isou)*aklskl)            &
             +   cssgr5*trrij* aikrjk
      epsij = -d2s3*cvara_ep(iel)*deltij(isou)

      w1(iel) = cromo(iel)*cell_f_vol(iel)*(pij+phiij1+phiij2+epsij)

      ! Implicit terms
      w2(iel) =   cell_f_vol(iel)/trrij*crom(iel)                          &
                * (cssgs1*cvara_ep(iel) + cssgr1*max(trprod,0.d0))

    ! EBRSM
    else

      xrij(1,1) = cvara_var(1,iel)
      xrij(2,2) = cvara_var(2,iel)
      xrij(3,3) = cvara_var(3,iel)
      xrij(1,2) = cvara_var(4,iel)
      xrij(2,3) = cvara_var(5,iel)
      xrij(1,3) = cvara_var(6,iel)
      xrij(2,1) = xrij(1,2)
      xrij(3,1) = xrij(1,3)
      xrij(3,2) = xrij(2,3)

      ! Compute the explicit term

      ! Calculation of the terms near the walls and et almost homogeneous
      ! of phi and epsilon

      ! Calculation of the term near the wall \f$ \Phi_{ij}^w \f$ --> W3
      !
      ! Phiw = -5.0 * (eps/k) * [ R*Xn + Xn^T*R - 0.5*tr(Xn*R)*(Xn + Id) ]
      !
      phiijw = 0.d0
      xnnd = d1s2*( xnal(iii)*xnal(jjj) + deltij(isou) )
      do kk = 1, 3
        phiijw = phiijw + xrij(iii,kk)*xnal(jjj)*xnal(kk)
        phiijw = phiijw + xrij(jjj,kk)*xnal(iii)*xnal(kk)
        do ll = 1, 3
          phiijw = phiijw - xrij(kk,ll)*xnal(kk)*xnal(ll)*xnnd
        enddo
      enddo
      phiijw = -5.d0*cvara_ep(iel)/trrij * phiijw

      ! Calculation of the almost homogeneous term \f$ \phi_{ij}^h \f$ --> W4
      phiij1 = -cvara_ep(iel)*cebms1*xaniso(iii,jjj)
      phiij2 = -cebmr1*trprod*xaniso(iii,jjj)                        &
                 +trrij*xstrai(iii,jjj)*(cebmr2-cebmr3*sqrt(aii))    &
                 +cebmr4*trrij   *(aiksjk-d2s3*deltij(isou)*aklskl)  &
                 +cebmr5*trrij   * aikrjk

      ! Calculation of \f $\e_{ij}^w \f$ --> W5 (Rotta model)
      ! Rij/k*epsilon
      epsijw =  xrij(iii,jjj)/trrij   *cvara_ep(iel)

      ! Calcul de \e_{ij}^h --> W6
      epsij =  d2s3*cvara_ep(iel)*deltij(isou)

      ! Calcul du terme source explicite de l'equation des Rij
      !  \f[ P_{ij} + (1-\alpha^3)\Phi_{ij}^w + \alpha^3\Phi_{ij}^h
      !            - (1-\alpha^3)\e_{ij}^w   - \alpha^3\e_{ij}^h  ]\f$ --> W1
      alpha3 = cvar_al(iel)**3

      w1(iel) = cell_f_vol(iel)*crom(iel)*(                            &
                 xprod(iii,jjj)                                        &
              + (1.d0-alpha3)*phiijw + alpha3*(phiij1+phiij2)          &
              - (1.d0-alpha3)*epsijw - alpha3*epsij)

      !  Implicit term

      if (irijco.eq.0) then
        ! Implicitation of epsijw
        ! (the factor 5 appears when we calculate \f$ Phi_{ij}^w - epsijw\f$)
        ! epsijw_imp = 5.d0 * (1.d0-alpha3)*cvara_ep(iel)/trrij           &
        !                   + (1.d0-alpha3)*cvara_ep(iel)/trrij
        epsijw_imp = 6.d0 * (1.d0-alpha3)*cvara_ep(iel)/trrij
      else
        ! FIXME
        epsijw_imp = 0.d0
      endif

      ! The term below corresponds to the implicit part of SSG
      ! in the context of elliptical weighting, it is multiplied by
      ! \f$ \alpha^3 \f$
      w2(iel) =   cell_f_vol(iel)*crom(iel)                             &
                * (  cebms1*cvara_ep(iel)/trrij*alpha3                  &
                   + cebmr1*max(trprod/trrij,0.d0)*alpha3               &
                   + epsijw_imp)

    endif ! EBRSM

    if (st_prv_id.ge.0) then
      c_st_prv(isou,iel) = c_st_prv(isou,iel) + w1(iel)
    else
      smbr(isou,iel) = smbr(isou,iel) + w1(iel)
      rovsdt(isou,isou,iel) = rovsdt(isou,isou,iel) + w2(iel)

      if (irijco.ne.0) then
        ! Careful ! Inversion of the order of the coefficients since
        ! rovsdt matrix is then used by a c function for the linear solving
        do jsou = 1, 6
          rovsdt(jsou,isou,iel) = rovsdt(jsou,isou,iel) + cell_f_vol(iel) &
                                  *crom(iel) * impl_drsm(isou,jsou)
        end do
      endif
    endif
  enddo
enddo

if (iturb.eq.32) then
  deallocate(grad)
endif

!===============================================================================
! Buoyancy source term
!===============================================================================

if (igrari.eq.1) then

  call field_get_id_try("rij_buoyancy", f_id)
  if (f_id.ge.0) then
    call field_get_val_v(f_id, cpro_buoyancy)
  else
    ! Allocate a work array
    allocate(buoyancy(6,ncelet))
    cpro_buoyancy => buoyancy
  endif

  call rijthe2(gradro, cpro_buoyancy)

  ! If we extrapolate the source terms: previous ST
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      do isou = 1, 6
        c_st_prv(isou,iel) = c_st_prv(isou,iel) + cpro_buoyancy(isou,iel) * cell_f_vol(iel)
      enddo
    enddo
  ! Otherwise smbr
  else
    do iel = 1, ncel
      do isou = 1, 6
        smbr(isou,iel) = smbr(isou,iel) + cpro_buoyancy(isou,iel) * cell_f_vol(iel)
      enddo
    enddo
  endif

  ! Free memory
  if (allocated(buoyancy)) deallocate(buoyancy)

  if (irijco.ne.0) then

    grav(1) = gx
    grav(2) = gy
    grav(3) = gz

    ! Implicit buoyancy term
    if (iscalt.gt.0) then
      call field_get_key_double(ivarfl(isca(iscalt)), ksigmas, turb_schmidt)
      const = -1.5d0 * cmu / turb_schmidt
    else
      const = -1.5d0 * cmu
    end if

    do iel = 1, ncel

      do jsou = 1, 6
        do isou = 1, 6
          impl_drsm(isou,jsou) = 0.0d0
        end do
      end do
      implmat2add(:,:) = 0.0d0

      kseps = (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel)) &
        / (2.0d0*cvara_ep(iel))

      xrij(1,1) = cvara_var(1,iel)
      xrij(2,2) = cvara_var(2,iel)
      xrij(3,3) = cvara_var(3,iel)
      xrij(1,2) = cvara_var(4,iel)
      xrij(2,3) = cvara_var(5,iel)
      xrij(1,3) = cvara_var(6,iel)
      xrij(2,1) = xrij(1,2)
      xrij(3,1) = xrij(1,3)
      xrij(3,2) = xrij(2,3)

      trrij = 0.5d0 * (xrij(1,1) + xrij(2,2) + xrij(3,3))

      gkks3 = 0.d0
      do jsou = 1, 3
        do isou = 1, 3
          gkks3 = gkks3 + grav(isou) * gradro(jsou,iel) * xrij(isou, jsou)
        enddo
      enddo
      gkks3 = const * kseps * gkks3 * crij3 * 2.d0 / 3.d0

      if (gkks3.le.0.d0) then
        ! Term "C3 tr(G) Id"
        ! Computing the inverse matrix of R^n
        ! Scaling by tr(R) in order to dodge inversion errors
        do isou = 1, 6
          matrn(isou) = cvara_var(isou,iel)/trrij
          oo_matrn(isou) = 0.0d0
        end do

        ! Inversing the matrix
        call symmetric_matrix_inverse(matrn, oo_matrn)
        do isou = 1, 6
          oo_matrn(isou) = oo_matrn(isou)/trrij
        end do

        do jsou = 1, 3
          do isou = 1, 3
            iii = t2v(isou,jsou)
            implmat2add(isou,jsou) = - 0.5d0 * gkks3 * oo_matrn(iii)
          end do
        end do

      endif

      gradchk = grav(1)*gradro(1,iel) + grav(2)*gradro(2,iel) + grav(3)*gradro(3,iel)
      if (gradchk .gt. 0.0d0) then

        ! Implicit term written as:
        !   Po . R^n+1 + R^n+1 . Po^t
        ! with Po proportional to "g (x) Grad rho"
        gradro_impl = const * (1.0d0-crij3) * kseps
        do jsou = 1, 3
          do isou = 1, 3
            implmat2add(isou,jsou) =   implmat2add(isou,jsou)   &
                                     - gradro_impl * grav(isou) * gradro(jsou,iel)
          end do
        end do
      end if

      ! Compute the 6x6 matrix A which verifies
      ! A . R = M . R + R . M^t
      call reduce_symprod33_to_66(implmat2add, impl_drsm)

      do isou = 1, 6
        do jsou = 1, 6
          ! Careful ! Inversion of the order of the coefficients since
          ! rovsdt matrix is then used by a c function for the linear solving
          rovsdt(jsou,isou,iel) =   rovsdt(jsou,isou,iel)       &
                                  + cell_f_vol(iel) * impl_drsm(isou,jsou)
        end do
      end do

    end do

  endif ! (irijco .ne. 0)

endif

!===============================================================================
! Diffusion term (Daly Harlow: generalized gradient hypothesis method)
!===============================================================================

! Symmetric tensor diffusivity (GGDH)
if (iand(vcopt%idften, ANISOTROPIC_RIGHT_DIFFUSION).ne.0) then

  call field_get_val_v(ivsten, visten)

  do iel = 1, ncel
    viscce(1,iel) = vcopt%idifft*visten(1,iel) + viscl(iel)
    viscce(2,iel) = vcopt%idifft*visten(2,iel) + viscl(iel)
    viscce(3,iel) = vcopt%idifft*visten(3,iel) + viscl(iel)
    viscce(4,iel) = vcopt%idifft*visten(4,iel)
    viscce(5,iel) = vcopt%idifft*visten(5,iel)
    viscce(6,iel) = vcopt%idifft*visten(6,iel)
  enddo

  iwarnp = vcopt%iwarni

  call vitens(viscce, iwarnp, weighf, weighb, viscf, viscb)

! Scalar diffusivity
else

  do iel = 1, ncel
    rctse = csrij * visct(iel) / cmu
    w1(iel) = viscl(iel) + vcopt%idifft*rctse
  enddo

  imvisp = vcopt%imvisf

  call viscfa(imvisp, w1, viscf, viscb)

endif

! Free memory

deallocate(w1, w2)

!--------
! Formats
!--------

 1000 format(/,'           Solving variable ', a8          ,/)

!----
! End
!----

return

end subroutine

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

!> \file resrij.f90
!>
!> \brief This subroutine prepares the solving of the Reynolds stress components
!> in \f$ R_{ij} - \varepsilon \f$ RANS (LRR) turbulence model.
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

subroutine resrij &
 ( ivar   ,                                                       &
   produc , gradro ,                                              &
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
use cs_f_interfaces
use rotation
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ivar

double precision produc(6, ncelet)
double precision gradro(3, ncelet)
double precision viscf(nfac), viscb(nfabor), viscce(6, ncelet)
double precision smbr(6, ncelet)
double precision rovsdt(6, 6, ncelet)
double precision weighf(2, nfac), weighb(nfabor)

! Local variables

integer          iel
integer          isou
integer          ii    , jj    , kk
integer          iflmas, iflmab
integer          iwarnp
integer          imvisp
integer          st_prv_id
integer          t2v(3,3)
integer          iv2t(6), jv2t(6)
integer          f_id
integer          key_t_ext_id, rot_id
integer          iroext

double precision trprod, trrij
double precision deltij(6), dij(3,3), cvara_r(3,3)
double precision thets , thetv
double precision matrot(3,3)
double precision d1s2, d1s3, d2s3
double precision ccorio
double precision rctse

character(len=80) :: label
double precision, allocatable, dimension(:,:), target :: buoyancy
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:,:) :: w2
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:,:), pointer :: visten
double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:,:), pointer :: cvar_var, cvara_var
double precision, dimension(:), pointer :: viscl
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
allocate(w1(ncelet))
allocate(w2(6,ncelet))

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

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_v(ivarfl(ivar), cvar_var)
call field_get_val_prev_v(ivarfl(ivar), cvara_var)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

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
! Production, Pressure-Strain correlation, dissipation
!===============================================================================

! Source term

!      (1-CRIJ2) Pij (for all components of Rij)

!      DELTAIJ*(2/3.CRIJ2.P+2/3.CRIJ1.EPSILON)
!                    (diagonal terms for R11, R22 et R33)

!      -DELTAIJ*2/3*EPSILON

!     If we extrapolate the source terms
!     We modify the implicit part:
!     In PHI1, we will only take rho CRIJ1 epsilon/k and not
!                                rho CRIJ1 epsilon/k (1-2/3 DELTAIJ)
!     It allow to keep  k^n instead of (R11^(n+1)+R22^n+R33^n)

!     If we extrapolate the source terms
if (st_prv_id.ge.0) then

  do iel = 1, ncel

    ! Half-traces of Prod and R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = 0.5d0 * (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel))

    do isou = 1, 6
      ! Calculation of Prod+Phi1+Phi2-Eps
      !  = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
      ! In c_st_prv:
      !  = rhoPij-C1rho eps/k(   -2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
      !  = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij           }
      c_st_prv(isou,iel) = c_st_prv(isou,iel) + cromo(iel) * cell_f_vol(iel)  &
                  *(   deltij(isou)*d2s3*                                     &
                       (  crij2*trprod                                        &
                        +(crij1-1.d0)* cvara_ep(iel)  )                       &
                  +(1.0d0-crij2)*produc(isou,iel)               )
      ! In smbr
      !  =       -C1rho eps/k(Rij         )
      !  = rho{                                     -C1eps/kRij}
      smbr(isou,iel) = smbr(isou,iel) + crom(iel) * cell_f_vol(iel)           &
        *( -crij1*cvara_ep(iel)/trrij * cvara_var(isou,iel))

      ! Calculation of the implicit part coming from Phil
      !  = C1rho eps/k(1        )
      rovsdt(isou,isou,iel) =   rovsdt(isou,isou,iel)                         &
                              + crom(iel) * cell_f_vol(iel)                   &
                                *crij1*cvara_ep(iel)/trrij*thetv
    enddo

    ! If we want to implicit a part of -C1rho eps/k(   -2/3k dij)
    ! FIXME: check if we want to use this or if it should be removed
    !        previously "isoluc = 2", never called

    if (.false.) then

      do isou = 1, 6
        ! We remove cromo
        !  =       -C1rho eps/k(   -1/3Rij dij)
        c_st_prv(isou,iel) =   c_st_prv(isou,iel)                             &
                             - cromo(iel) * cell_f_vol(iel)                   &
                               * (deltij(isou)*d1s3*crij1*cvara_ep(iel)/trrij &
                                  * cvara_var(isou,iel))
        ! We add to smbr (with crom)
        !  =       -C1rho eps/k(   -1/3Rij dij)
        smbr(isou,iel) =   smbr(isou,iel)                                     &
                         + crom(iel) * cell_f_vol(iel)                        &
                         * (deltij(isou)*d1s3*crij1*cvara_ep(iel)/trrij       &
                            * cvara_var(isou,iel))
        !  We add to rovsdt (woth crom)
        ! =        C1rho eps/k(   -1/3    dij)
        rovsdt(isou,isou,iel) =   rovsdt(isou,isou,iel)                       &
                                + crom(iel) * cell_f_vol(iel)                 &
                                *(deltij(isou)*d1s3*crij1*cvara_ep(iel)/trrij)
      enddo

    endif

  enddo

! If we do not extrapolate the source terms
else

  do iel = 1, ncel

    !     Half-traces of Prod and R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = 0.5d0 * (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel))

    do isou = 1, 6
      ! Calculation of Prod+Phi1+Phi2-Eps
      !  = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
      !  = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij-C1eps/kRij}
      smbr(isou,iel) = smbr(isou,iel) + crom(iel) * cell_f_vol(iel)           &
        *(   deltij(isou)*d2s3*                                               &
             (  crij2*trprod                                                  &
              +(crij1-1.d0)* cvara_ep(iel)  )                                 &
           +(1.0d0-crij2)*produc(isou,iel)                                    &
           -crij1*cvara_ep(iel)/trrij * cvara_var(isou,iel)  )

      ! Calculation of the implicit part coming from Phi1
      !  = C1rho eps/k(1-1/3 dij)
      rovsdt(isou,isou,iel) =   rovsdt(isou,isou,iel)                         &
                              + crom(iel) * cell_f_vol(iel)                   &
                                *(1.d0-d1s3*deltij(isou))*crij1               &
                                *cvara_ep(iel)/trrij
    enddo

  enddo

endif

!===============================================================================
! Coriolis terms in the Phi1 and production
!===============================================================================

if (icorio.eq.1 .or. iturbo.eq.1) then

  do iel = 1, ncel

    rot_id = icorio
    if (iturbo.eq.1) rot_id = irotce(iel)

    if (rot_id .lt. 1) cycle

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

    ! Compute Gij: (i,j) component of the Coriolis production
    do isou = 1, 6
      ii = iv2t(isou)
      jj = jv2t(isou)

      w2(isou,iel) = 0.d0
      do kk = 1, 3

        w2(isou,iel) = w2(isou,iel) - ccorio*(  matrot(ii,kk)*cvara_r(jj,kk) &
                                    + matrot(jj,kk)*cvara_r(ii,kk) )
      enddo
    enddo
  enddo

  ! Coriolis contribution in the Phi1 term: (1-C2/2)Gij
  if (icorio.eq.1) then
    do iel = 1, ncel
      do isou = 1, 6
        w2(isou,iel) = crom(iel)*cell_f_vol(iel)*(1.d0 - 0.5d0*crij2)*w2(isou,iel)
      enddo
    enddo
  endif

  ! If source terms are extrapolated
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      do isou = 1, 6
        c_st_prv(isou,iel) = c_st_prv(isou,iel) + w2(isou,iel)
      enddo
    enddo
  ! Otherwise, directly in smbr
  else
    do iel = 1, ncel
      do isou = 1, 6
        smbr(isou,iel) = smbr(isou,iel) + w2(isou,iel)
      enddo
    enddo
  endif

endif

!===============================================================================
! Wall echo terms
!===============================================================================

if (irijec.eq.1) then !todo

  do iel = 1, ncel
    do isou = 1, 6
      w2(isou,iel) = 0.d0
    enddo
  enddo

  call rijech2(produc, w2)

  ! If we extrapolate the source terms: c_st_prv
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      do isou = 1, 6
        c_st_prv(isou,iel) = c_st_prv(isou,iel) + w2(isou,iel)
      enddo
    enddo
  ! Otherwise smbr
  else
    do iel = 1, ncel
      do isou = 1, 6
        smbr(isou,iel) = smbr(isou,iel) + w2(isou,iel)
      enddo
    enddo
  endif

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
    trrij = 0.5d0 * (cvara_var(1,iel) + cvara_var(2,iel) + cvara_var(3,iel))
    rctse = crom(iel) * csrij * trrij**2 / cvara_ep(iel)
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

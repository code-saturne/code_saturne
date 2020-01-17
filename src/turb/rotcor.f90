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

!===============================================================================
! Function:
! ---------

!> \file rotcor.f90
!>
!> \brief Computing rotation/curvature correction for eddy-viscosity models.
!>        The subroutine is called for the linear eddy viscosity RANS models,
!>        when the option irccor = 1 is verified.
!>
!> Two type of rotation/curvature correction are computed, depending on
!> the specific eddy-viscosity model:
!>
!> - itycor = 1: - Cazalbou correction (variable Ce2 coefficient in the
!>                 destruction term of dissipation equation)
!>               - default correction for \f$ k - \epsilon \f$ type models,
!>                 including elliptic relaxation/blending models
!>                 (iturb = 20, 21, 50 or 51)
!>
!> - itycor = 2: - Spalart-Shur correction (production terms are multiplied
!>                 by a rotation function)
!>               - default correction for \f$ k - \omega \f$ SST or
!>                 Spalart-Allmaras (iturb = 60 or 70)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     dt            time step (per cell)
!> \param[out]    rotfct        rotation function of Spalart-Shur correction
!>                               at cell center
!> \param[out]    ce2rc         modified ce2 coeficient of Cazalbou correction
!>                               at cell center
!______________________________________________________________________________!

subroutine rotcor &
 ( dt     , rotfct , ce2rc  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use cstnum
use parall
use period
use mesh
use field
use field_operator
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision dt(ncelet)
double precision rotfct(ncel), ce2rc(ncel)

! Local variables

integer          iel, ivar, ifac, isou, f_id
integer          iccocg, inc, imrgrp, nswrgp, imligp, iwarnp
integer          istrai(3,3), ivorab(3,3)
integer          ii, jj, kk, iprev

double precision epsrgp, climgp, extrap
double precision matrot(3,3), sigvor(3,3)
double precision dsijdt, trrota, wiksjk, rstar, echtm2
double precision stilde, wtilde, rotild
double precision xe, xk

double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: strain, vortab
double precision, allocatable, dimension(:,:) :: grdsij
double precision, allocatable, dimension(:) :: coeas, coebs
double precision, allocatable, dimension(:) :: brtild, eta1, eta2

double precision, dimension(:,:), pointer :: vela

double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:,:), pointer :: cpro_straio

integer          ipass
data             ipass /0/
save             ipass

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

call field_get_val_v(iprpfl(istraio), cpro_straio)

ipass = ipass + 1
if(ipass.eq.1) then
  do isou = 1, 6
    do iel = 1, ncelet
      cpro_straio(isou,iel) = 0.d0
    enddo
  enddo
endif

! Map field arrays
call field_get_val_prev_v(ivarfl(iu), vela)

if (itycor.eq.1) then
  call field_get_val_prev_s(ivarfl(ik), cvara_k)
  call field_get_val_prev_s(ivarfl(iep), cvara_ep)
else if (itycor.eq.2) then
  if (iturb.eq.60) call field_get_val_prev_s(ivarfl(iomg), cvara_omg)
endif

if (icorio.eq.1) then
  ! In case of rotating frame, all cells belong to the same "rotor"
  call coriolis_t(1, 1.d0, matrot)
else
  do ii = 1, 3
    do jj = 1, 3
      matrot(ii,jj) = 0.d0
    enddo
  enddo
endif

!===============================================================================
! 1. Preliminary calculations
!===============================================================================

!-------------------------------------------------------------------------------
! 1.1 Compute the strain rate and absolute vorticity tensor
!-------------------------------------------------------------------------------

! Allocate temporary arrays
allocate(strain(ncelet,6),vortab(ncelet,3))

allocate(gradv(3, 3, ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, 0, inc, gradv)

! Compute the strain rate tensor (symmetric)
!          S_ij = 0.5(dU_i/dx_j+dU_j/dx_i)
! and the absolute vorticity tensor (anti-symmetric)
!          W_ij = 0.5(dU_i/dx_j-dU_j/dx_i) + e_imj*Omega_m

! Only the non zero components in the upper triangle are stored

do iel = 1, ncel

  ! S11
  strain(iel,1) = gradv(1, 1, iel)
  ! S22
  strain(iel,2) = gradv(2, 2, iel)
  ! S33
  strain(iel,3) = gradv(3, 3, iel)
  ! S12
  strain(iel,4) = 0.5d0*(gradv(2, 1, iel) + gradv(1, 2, iel))
  ! S13
  strain(iel,5) = 0.5d0*(gradv(3, 1, iel) + gradv(1, 3, iel))
  ! S23
  strain(iel,6) = 0.5d0*(gradv(3, 2, iel) + gradv(2, 3, iel))
  ! W12
  vortab(iel,1) = 0.5d0*(gradv(2, 1, iel) - gradv(1, 2, iel)) + matrot(1,2)
  ! W13
  vortab(iel,2) = 0.5d0*(gradv(3, 1, iel) - gradv(1, 3, iel)) + matrot(1,3)
  ! W23
  vortab(iel,3) = 0.5d0*(gradv(3, 2, iel) - gradv(2, 3, iel)) + matrot(2,3)

enddo

! Free memory (strain and vortab arrays are deallocated later)
deallocate(gradv)

!-------------------------------------------------------------------------------
! 1.2 Computation of :
!
!   brtild = 2.W_ik.S_jk(DS_ij/Dt + (e_imn.S_jn + e_jmn.S_in)*Omega_m)
!     eta1 = S_ij.S_ij
!     eta2 = W_ij.W_ij
!
! ------------------------------------------------------------------------------

! Allocate temporary arrays
allocate(grdsij(3,ncelet))
allocate(coeas(nfabor),coebs(nfabor))
allocate(brtild(ncel),eta1(ncel),eta2(ncel))

! Index connectivity

! istrai(i,j) : position of the (i,j) component of the tensor
!               in the strain and straio arrays
! ivorab(i,j) : position of the (i,j) component of the tensor
!               in the vortab array
! sigvor(i,j) : sign of the (i,j) component of the absolute vorticity tensor
!               = 1  if i > j
!               = -1 if i < j
!               = 0  if i = j

istrai(1,1) = 1
istrai(2,2) = 2
istrai(3,3) = 3
istrai(1,2) = 4
istrai(1,3) = 5
istrai(2,3) = 6
istrai(2,1) = istrai(1,2)
istrai(3,1) = istrai(1,3)
istrai(3,2) = istrai(2,3)

ivorab(1,1) = 1
ivorab(2,2) = 1
ivorab(3,3) = 1
ivorab(1,2) = 1
ivorab(1,3) = 2
ivorab(2,3) = 3
ivorab(2,1) = ivorab(1,2)
ivorab(3,1) = ivorab(1,3)
ivorab(3,2) = ivorab(2,3)

do ii = 1, 3
  do jj = 1, 3
    if (ii.lt.jj) then
      sigvor(ii,jj) = 1.d0
    elseif (ii.eq.jj) then
      sigvor(ii,jj) = 0.d0
    else
      sigvor(ii,jj) = -1.d0
    endif
  enddo
enddo

! Boundary conditions for S_ij -> homogeneous Neumann
do ifac = 1, nfabor
  coeas(ifac) = 0.d0
  coebs(ifac) = 1.d0
enddo

do iel = 1, ncel
  brtild(iel) = 0.d0
  eta1(iel) = 0.d0
  eta2(iel) = 0.d0
enddo

do ii = 1, 3

  do jj = 1, 3

    iccocg = 1
    inc = 1

    if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
      ivar = ik
    elseif (iturb.eq.70) then
      ivar = inusa
    endif

    call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

    imrgrp = vcopt%imrgra
    nswrgp = vcopt%nswrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    epsrgp = vcopt%epsrgr
    climgp = vcopt%climgr
    extrap = vcopt%extrag

    f_id = -1

    call gradient_s                                                &
  ( f_id   , imrgrp , inc    , iccocg , nswrgp , imligp ,          &
    iwarnp , epsrgp , climgp , extrap ,                            &
    strain(1,istrai(ii,jj))  , coeas  , coebs ,                    &
    grdsij )

    do iel = 1, ncel

      ! material derivative of S_ij
      if (idtvar.lt.0) then
        dsijdt = 0.d0
      else
        dsijdt = ((strain(iel,istrai(ii,jj))        &
                 - cpro_straio(istrai(ii,jj),iel))  &
                /dt(iel))
      endif
      dsijdt = dsijdt + vela(1,iel)*grdsij(1,iel)        &
                      + vela(2,iel)*grdsij(2,iel)        &
                      + vela(3,iel)*grdsij(3,iel)

      ! (e_imn.S_jn+e_jmn.S_in)*Omega_m term
      trrota = 0.d0
      do kk = 1, 3
        trrota = trrota                                       &
              + matrot(ii,kk)*strain(iel,istrai(jj,kk))       &
              + matrot(jj,kk)*strain(iel,istrai(ii,kk))
      enddo

      ! W_ik.S_jk term
      wiksjk = 0.d0
      do kk = 1, 3
        wiksjk = wiksjk + sigvor(ii,kk)*vortab(iel,ivorab(ii,kk))   &
              *strain(iel,istrai(jj,kk))
      enddo

      ! brtild, eta1, eta2 (see the definitions above)
      brtild(iel) = brtild(iel)                 &
           + 2.d0*wiksjk*(dsijdt + trrota)
      eta1(iel) = eta1(iel)                   &
           + strain(iel,istrai(ii,jj))**2
      eta2(iel) = eta2(iel)                   &
           + (sigvor(ii,jj)*vortab(iel,ivorab(ii,jj)))**2

    enddo

  enddo

enddo

!===============================================================================
! 2. Effective computation of the rotation correction
!===============================================================================

if (itycor.eq.1) then

!-------------------------------------------------------------------------------
! 2.1 Cazalbou correction
!-------------------------------------------------------------------------------

  do iel = 1, ncel

    ! Computation of stilde = sqrt(2.S_ij.S_ij) et wtilde = sqrt(W_ij.W_ij/2)
    stilde = max(sqrt(eta1(iel)*2.d0),1.d-15)
    wtilde = max(sqrt(eta2(iel)/2.d0),1.d-15)

    xk = max(cvara_k(iel),1.d-15)
    xe = max(cvara_ep(iel),1.d-15)
    rotild = xe/wtilde/xk
    brtild(iel) = -brtild(iel)*xk/xe/stilde**3

    ! Variable C_eps_2 coefficient of Cazalbou
    ce2rc(iel) = ccaze2                                 &
         + (ccaze2 - 1.d0)/(1 + ccaza*rotild**1.5d0)             &
         + ccaze2*ccazsc*stilde*xk/xe                         &
         *(tanh(ccazb*brtild(iel) + ccazc) - ccazd)

    ce2rc(iel) = max(ce2rc(iel),0.d0)

  enddo

elseif (itycor.eq.2) then

!-------------------------------------------------------------------------------
! 2.2 Spalart-Shur correction (including modifications of
!   Smirnov & Menter, ASME, 2009)
!-------------------------------------------------------------------------------

  do iel = 1, ncel

    ! Computation of stilde = 2.S_ij.S_ij and wtilde = 2.W_ij.W_ij
    stilde = max(eta1(iel)*2.d0,1.d-15)
    wtilde = max(eta2(iel)*2.d0,1.d-15)

    echtm2 = stilde

    if (iturb.eq.60)  echtm2 = max(echtm2,cmu*cvara_omg(iel)**2)

    brtild(iel) = brtild(iel)/sqrt(wtilde*echtm2**3)

    rstar = sqrt(stilde)/sqrt(wtilde)

    ! Rotation function of Spalart & Shur
    rotfct(iel) = (1.d0 + cssr1)*2.d0*rstar/(1.d0 + rstar)     &
         *(1.d0 - cssr3*atan(cssr2*brtild(iel))) - cssr1

    rotfct(iel) = min(max(rotfct(iel),0.d0),1.25d0)

  enddo

endif

!===============================================================================
! 3. Finalizations
!===============================================================================

! Save de strain rate tensor for the next time step
if (idtvar.ge.0) then
  do isou = 1, 6
    do iel = 1, ncelet
      cpro_straio(isou,iel) = strain(iel,isou)
    enddo
  enddo
endif

! Free memory
deallocate(strain,vortab)
deallocate(grdsij)
deallocate(coeas,coebs)
deallocate(brtild,eta1,eta2)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine

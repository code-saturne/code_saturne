!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file turrij.f90
!>
!> \brief Solving the \f$ R_{ij} - \epsilon \f$ for incompressible flows or
!>  slightly compressible flows for one time step.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        index of the ncepdp cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        mass source type for the variables (cf. cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     tslagr        coupling term of the lagangian module
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        values of the variables associated to the
!>                               mass source
!>                               (for ivar=ipr, smacel is the mass flux)
!_______________________________________________________________________________

subroutine turrij &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   tslagr ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use lagdim, only: ntersl
use numvar
use entsor
use cstphy
use optcal
use lagran
use pointe, only: rvoid1
use ppincl
use mesh
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp)

integer, dimension(ncesmp,nvar), target :: itypsm

double precision dt(ncelet)
double precision ckupdc(ncepdp,6)

double precision, dimension(ncesmp,nvar), target ::  smacel
double precision, dimension(ncelet,ntersl), target :: tslagr
double precision, dimension(:), pointer :: bromo, cromo

! Local variables

integer          ifac  , iel   , ivar  , isou  , ii
integer          inc   , iccocg
integer          iwarnp, iclip
integer          nswrgp, imligp
integer          f_id0
integer          iitsla
integer          iprev
double precision epsrgp, climgp, extrap
double precision rhothe

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: produc
double precision, allocatable, dimension(:,:) :: gradro
double precision, allocatable, dimension(:,:) :: smbrts
double precision, allocatable, dimension(:,:,:) ::rovsdtts

double precision, pointer, dimension(:) :: tslage, tslagi
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvara_r12, cvara_r13, cvara_r23
double precision, dimension(:), pointer :: cvara_scalt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

tslage => rvoid1
tslagi => rvoid1

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(ncelet), rovsdt(ncelet))
allocate(gradv(3, 3, ncelet))
if (irijco.eq.1) then
  allocate(smbrts(6,ncelet))
  allocate(rovsdtts(6,6,ncelet))
endif

! Allocate other arrays, depending on user options
if (iturb.eq.30) then
  allocate(produc(6,ncelet))
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

if (iturb.eq.30) then
  call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
  call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
  call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
  call field_get_val_prev_s(ivarfl(ir12), cvara_r12)
  call field_get_val_prev_s(ivarfl(ir13), cvara_r13)
  call field_get_val_prev_s(ivarfl(ir23), cvara_r23)
endif

if(iwarni(iep).ge.1) then
  if (iturb.eq.30) then
    write(nfecra,1000)
  elseif (iturb.eq.31) then
    write(nfecra,1001)
  else
    write(nfecra,1002)
  endif
endif

!===============================================================================
! 2.1 Compute the velocity gradient
!===============================================================================

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

!===============================================================================
! 2.2 Compute the production term for Rij LRR (iturb =30)
!===============================================================================

if (iturb.eq.30) then

  do ii = 1 , 6
    do iel = 1, ncel
      produc(ii,iel) = 0.0d0
    enddo
  enddo

  do iel = 1 , ncel

    ! grad u

    produc(1,iel) = produc(1,iel)                                  &
                  - 2.0d0*(cvara_r11(iel)*gradv(1, 1, iel) +           &
                           cvara_r12(iel)*gradv(2, 1, iel) +           &
                           cvara_r13(iel)*gradv(3, 1, iel) )

    produc(4,iel) = produc(4,iel)                                  &
                  - (cvara_r12(iel)*gradv(1, 1, iel) +                 &
                     cvara_r22(iel)*gradv(2, 1, iel) +                 &
                     cvara_r23(iel)*gradv(3, 1, iel) )

    produc(6,iel) = produc(6,iel)                                  &
                  - (cvara_r13(iel)*gradv(1, 1, iel) +                 &
                     cvara_r23(iel)*gradv(2, 1, iel) +                 &
                     cvara_r33(iel)*gradv(3, 1, iel) )

    ! grad v

    produc(2,iel) = produc(2,iel)                                  &
                  - 2.0d0*(cvara_r12(iel)*gradv(1, 2, iel) +           &
                           cvara_r22(iel)*gradv(2, 2, iel) +           &
                           cvara_r23(iel)*gradv(3, 2, iel) )

    produc(4,iel) = produc(4,iel)                                  &
                  - (cvara_r11(iel)*gradv(1, 2, iel) +                 &
                     cvara_r12(iel)*gradv(2, 2, iel) +                 &
                     cvara_r13(iel)*gradv(3, 2, iel) )

    produc(5,iel) = produc(5,iel)                                  &
                  - (cvara_r13(iel)*gradv(1, 2, iel) +                 &
                     cvara_r23(iel)*gradv(2, 2, iel) +                 &
                     cvara_r33(iel)*gradv(3, 2, iel) )

    ! grad w

    produc(3,iel) = produc(3,iel)                                  &
                  - 2.0d0*(cvara_r13(iel)*gradv(1, 3, iel) +           &
                           cvara_r23(iel)*gradv(2, 3, iel) +           &
                           cvara_r33(iel)*gradv(3, 3, iel) )

    produc(6,iel) = produc(6,iel)                                  &
                  - (cvara_r11(iel)*gradv(1, 3, iel) +                 &
                     cvara_r12(iel)*gradv(2, 3, iel) +                 &
                     cvara_r13(iel)*gradv(3, 3, iel) )

    produc(5,iel) = produc(5,iel)                                  &
                  - (cvara_r12(iel)*gradv(1, 3, iel) +                 &
                     cvara_r22(iel)*gradv(2, 3, iel) +                 &
                     cvara_r23(iel)*gradv(3, 3, iel) )

  enddo

endif

!===============================================================================
! 3. Compute the density gradient for buoyant terms
!===============================================================================

! Note that the buoyant term is normally expressed in temr of
! (u'T') or (u'rho') here modelled with a GGDH:
!   (u'rho') = C * k/eps * R_ij Grad_j(rho)

! Buoyant term for the Atmospheric module
! (function of the potential temperature)
if (igrari.eq.1 .and. ippmod(iatmos).ge.1) then
  ! Allocate a temporary array for the gradient calculation
  ! Warning, grad(theta) here
  allocate(gradro(3,ncelet))

  call field_get_val_s(icrom, cromo)

  call field_get_val_prev_s(ivarfl(isca(iscalt)), cvara_scalt)

  iprev = 1
  inc = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(isca(iscalt)), 1, imrgra, inc, &
                             iccocg,                               &
                             gradro)

  ! gradro stores: rho grad(theta)/theta
  do iel = 1, ncel
    rhothe = cromo(iel)/cvara_scalt(iel)
    gradro(1, iel) = rhothe*gradro(1, iel)
    gradro(2, iel) = rhothe*gradro(2, iel)
    gradro(3, iel) = rhothe*gradro(3, iel)
  enddo

else if (igrari.eq.1) then
  ! Allocate a temporary array for the gradient calculation
  allocate(gradro(3,ncelet))

! Boundary conditions: Dirichlet romb
!   We use viscb to store the relative coefficient of rom
!   We impose in Dirichlet (coefa) the value romb

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

! The choice below has the advantage to be simple

  nswrgp = nswrgr(ir11)
  imligp = imligr(ir11)
  iwarnp = iwarni(ir11)
  epsrgp = epsrgr(ir11)
  climgp = climgr(ir11)
  extrap = extrag(ir11)

  f_id0 = -1

  ! If we extrapolate the source terms and rho, we use cpdt rho^n
  if(isto2t.gt.0.and.iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
    call field_get_val_prev_s(ibrom, bromo)
  else
    call field_get_val_s(icrom, cromo)
    call field_get_val_s(ibrom, bromo)
  endif

  call gradient_s                                                 &
 ( f_id0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   cromo  , bromo  , viscb           ,                            &
   gradro )

endif

!===============================================================================
! 4.  Loop on the variables Rij (6 variables)
!     The order is R11 R22 R33 R12 R13 R23 (The place of those variables
!      is IR11.    ..
!     We solve the equation in a routine similar to covofi.f90
!===============================================================================

if (irijco.eq.1) then
  ivar = irij

  ! Rij-epsilon standard (LRR)
  if (iturb.eq.30) then !TODO
    call resrij &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   produc , gradro ,                                              &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   tslage , tslagi ,                                              &
   smbr   , rovsdt )

  ! Rij-epsilon SSG or EBRSM
  elseif (iturb.eq.31.or.iturb.eq.32) then

      call resssg2 &
      !==========
    ( nvar   , nscal  , ncepdp , ncesmp ,                            &
      ivar   ,                                              &
      icepdc , icetsm , itypsm ,                                     &
      dt     ,                                                       &
      gradv  , gradro ,                                              &
      ckupdc , smacel ,                                              &
      viscf  , viscb  ,                                              &
      tslage , tslagi ,                                              &
      smbrts   , rovsdtts )
  endif
else
  do isou = 1, 6
    if    (isou.eq.1) then
      ivar   = ir11
    elseif(isou.eq.2) then
      ivar   = ir22
    elseif(isou.eq.3) then
      ivar   = ir33
    elseif(isou.eq.4) then
      ivar   = ir12
    elseif(isou.eq.5) then
      ivar   = ir23
    elseif(isou.eq.6) then
      ivar   = ir13
    endif

    if (iilagr.eq.2) then
      iitsla = itsr11 + (isou-1)
      tslage => tslagr(1:ncelet,iitsla)
      tslagi => tslagr(1:ncelet,itsli)
    endif

    ! Rij-epsilon standard (LRR)
    if (iturb.eq.30) then
      call resrij &
      !==========
   ( nvar   , nscal  , ncepdp , ncesmp ,                            &
     ivar   , isou   ,                                              &
     icepdc , icetsm , itypsm ,                                     &
     dt     ,                                                       &
     produc , gradro ,                                              &
     ckupdc , smacel ,                                              &
     viscf  , viscb  ,                                              &
     tslage , tslagi ,                                              &
     smbr   , rovsdt )

    ! Rij-epsilon SSG or EBRSM
    elseif (iturb.eq.31.or.iturb.eq.32) then
        call resssg &
        !==========
      ( nvar   , nscal  , ncepdp , ncesmp ,                            &
        ivar   , isou   ,                                              &
        icepdc , icetsm , itypsm ,                                     &
        dt     ,                                                       &
        gradv  , gradro ,                                              &
        ckupdc , smacel ,                                              &
        viscf  , viscb  ,                                              &
        tslage , tslagi ,                                              &
        smbr   , rovsdt )
    endif

  enddo
endif
!===============================================================================
! 5. Solve Epsilon
!===============================================================================

call reseps &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   gradv  , produc , gradro ,                                     &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   smbr   , rovsdt )

!===============================================================================
! 6. Clipping
!===============================================================================

if (iturb.eq.32) then
  iclip = 1
else
  iclip = 2
endif

if (irijco.eq.1) then
  call clprij2 &
  !==========
   ( ncelet , ncel   , nvar   ,                                     &
     iclip  )

else

  call clprij &
  !==========
   ( ncelet , ncel   , nvar   ,                                     &
     iclip  )
endif

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, rovsdt)
if (allocated(gradro)) deallocate(gradro)
if (allocated(produc)) deallocate(produc)
deallocate(gradv)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                      &
'   ** Resolution du Rij-EPSILON LRR'                   ,/,&
'      -----------------------------'             ,/)
 1001 format(/,                                      &
'   ** Resolution du Rij-EPSILON SSG'                   ,/,&
'      -----------------------------'             ,/)
 1002 format(/,                                      &
'   ** Resolution du Rij-EPSILON EBRSM'                 ,/,&
'      -------------------------------'           ,/)

#else

 1000 format(/,                                      &
'   ** Solving Rij-EPSILON LRR'                   ,/,&
'      -----------------------'                   ,/)
 1001 format(/,                                      &
'   ** Solving Rij-EPSILON SSG'                   ,/,&
'      -----------------------'                   ,/)
 1002 format(/,                                      &
'   ** Solving Rij-EPSILON EBRSM'                 ,/,&
'      -------------------------'                 ,/)

#endif

!----
! End
!----

return

end subroutine

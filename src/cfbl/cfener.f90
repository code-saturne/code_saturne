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

!> \file cfener.f90
!> \brief Perform the solving of the convection/diffusion equation (with
!> eventual source terms) for total energy over a time step. It is the third
!> step of the compressible algorithm at each time iteration.
!>
!> Please refer to the <a href="../../theory.pdf#cfener"><b>cfener</b></a> section
!> of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     iscal         scalar number
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!_______________________________________________________________________________

subroutine cfener &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe, only:rvoid1
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field
use field_operator
use cs_c_bindings
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          use_previous

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)

! Local variables

character(len=80) :: chaine
integer          ivar  , ivarsp, iesp
integer          ifac  , iel
integer          init  , iii
integer          ifcvsl, iflmas, iflmab
integer          icvflb
integer          imrgrp, nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision sclnor, thetap, epsrsp, relaxp
double precision turb_schmidt

integer          inc    , iccocg , imucpp , idftnp , iswdyp
integer          f_id0  , ii, jj
integer          iel1  , iel2
integer          iterns

double precision flux, yip, yjp, gradnb, tip
double precision dijpfx, dijpfy, dijpfz, pnd  , pip   , pjp
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz
double precision normp

double precision mk, cpk, cvk
double precision rvoid(1)

double precision, allocatable, dimension(:) :: wb, iprtfl, bprtfl, viscf, viscb
double precision, allocatable, dimension(:) :: viscfk, viscbk
double precision, allocatable, dimension(:) :: dpvar, smbrs, rovsdt
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: grad_dd
double precision, allocatable, dimension(:) :: gapinj
double precision, allocatable, dimension(:) :: w1, w7, w9
double precision, allocatable, dimension(:) :: kspe, btemp

double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, cromo
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: coefat, coefbt, coefayk, coefbyk
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:,:), pointer :: vel

double precision, dimension(:), pointer :: cvar_pr, cvar_energ, cvara_energ
double precision, dimension(:), pointer :: cvar_tempk, cvar_yk
double precision, dimension(:), pointer :: visct, cpro_cp, cpro_cv, cpro_viscls

type(gas_mix_species_prop) :: s_k
type(var_cal_opt) :: vcopt_u, vcopt_p, vcopt_e
double precision, dimension(:), pointer :: cvar_fracv, cvar_fracm, cvar_frace

!===============================================================================

! Computation number and post-treatment number of the scalar total energy
ivar   = isca(iscal)

! Map field arrays
call field_get_val_prev_s(ivarfl(ivar), cvara_energ)
call field_get_val_s(ivarfl(ivar), cvar_energ)
call field_get_val_s(ivarfl(isca(itempk)), cvar_tempk)
call field_get_val_v(ivarfl(iu), vel)

if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
else
  cpro_cp => rvoid1
endif

if (icv.ge.0) then
  call field_get_val_s(icv, cpro_cv)
else
  cpro_cv => rvoid1
endif

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array
allocate(wb(nfabor))
allocate(smbrs(ncelet), rovsdt(ncelet))

! Allocate work arrays
allocate(grad(3,ncelet))
allocate(w7(ncelet), w9(ncelet))
allocate(dpvar(ncelet))

! Physical property numbers
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_prev_s(icrom, cromo)

call field_get_val_s(ivarfl(ipr), cvar_pr)

call field_get_val_s(ivisct, visct)

if (ippmod(icompf).eq.2) then
  call field_get_val_s(ivarfl(isca(ifracv)), cvar_fracv)
  call field_get_val_s(ivarfl(isca(ifracm)), cvar_fracm)
  call field_get_val_s(ivarfl(isca(ifrace)), cvar_frace)
else
  cvar_fracv => rvoid1
  cvar_fracm => rvoid1
  cvar_frace => rvoid1
endif

call field_get_key_int(ivarfl(ivar), kimasf, iflmas)
call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

! Prints
call field_get_label(ivarfl(ivar), chaine)

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)
call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt_e)

if(vcopt_e%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

! Barotropic version
if (ippmod(icompf).eq.1) goto 33

!===============================================================================
! 2. Source terms
!===============================================================================

! Theta-scheme:
! for now, theta=1 is assumed and the theta-scheme is not implemented

! --> Initialization

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! HEAT VOLUMIC SOURCE TERM: RHO * PHI *VOLUME
! =================================

call ustssc                                                                    &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                                          &
  iscal  ,                                                                     &
  icepdc , icetsm , itypsm ,                                                   &
  dt     ,                                                                     &
  ckupdc , smacel , smbrs  , rovsdt )

! C version
call user_source_terms(ivarfl(isca(iscal)), smbrs, rovsdt)

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + rovsdt(iel)*cvar_energ(iel)
  rovsdt(iel) = max(-rovsdt(iel),zero)
enddo


! MASS SOURCE TERMS
! =================

! GAMMA(IEL) = SMACEL(IEL,IPR)

! Implicit term : GAMMA*VOLUME
!                                                        n
! Explicit term : GAMMA*VOLUME*e   - GAMMA*VOLUME*e
!                                     inj
if (ncesmp.gt.0) then
  iterns = 1
  allocate(gapinj(ncelet))
  call catsma ( ncelet , ncel , ncesmp , iterns ,                              &
                isno2t,                                                        &
                icetsm , itypsm(1,ivar),                                       &
                cell_f_vol    , cvara_energ   , smacel(1,ivar),                &
                smacel(1,ipr) , smbrs , rovsdt , gapinj)
  deallocate(gapinj)
endif


!                          RHO*VOLUME
! UNSTEADY IMPLICIT TERM : ----------
! ======================       DT

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                                    &
                + vcopt_e%istat*(cromo(iel)/dt(iel))*cell_f_vol(iel)
enddo

!                                       __        v
!     VISCOUS DISSIPATION TERM        : >  ((SIGMA *U).n)  *S
!     ========================          --               ij  ij

if (vcopt_u%idiff.ge. 1) then
  call cfdivs(smbrs, vel)
endif

!                              __   P        n+1
! PRESSURE TRANSPORT TERM  : - >  (---)  *(Q    .n)  *S
! =======================      --  RHO ij   pr     ij  ij

allocate(iprtfl(nfac))
allocate(bprtfl(nfabor))

! No reconstruction yet

! Internal faces
do ifac = 1, nfac
  iel1 = ifacel(1,ifac)
  iel2 = ifacel(2,ifac)
  iprtfl(ifac) =                                                             &
       - cvar_pr(iel1)/crom(iel1) * 0.5d0*(imasfl(ifac) +abs(imasfl(ifac)))  &
       - cvar_pr(iel2)/crom(iel2) * 0.5d0*(imasfl(ifac) -abs(imasfl(ifac)))
enddo


! Boundary faces: for the faces where a flux (Rusanov or analytical) has been
! computed, the standard contribution is replaced by this flux in bilsc2.

call field_get_coefa_s(ivarfl(ipr), coefa_p)
call field_get_coefb_s(ivarfl(ipr), coefb_p)

do ifac = 1, nfabor
  if (icvfli(ifac).eq.0) then
    iel = ifabor(ifac)
    bprtfl(ifac) = - bmasfl(ifac)                                            &
                    * (coefa_p(ifac) + coefb_p(ifac)*cvar_pr(iel))           &
                    / brom(ifac)
  else
    bprtfl(ifac) = 0.d0
  endif
enddo

!     Divergence
init = 0
call divmas(init, iprtfl, bprtfl, smbrs)

deallocate(iprtfl)
deallocate(bprtfl)

! GRAVITATION FORCE TERM : RHO*g.U *VOLUME
! ======================

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + crom(iel)*cell_f_vol(iel)                          &
                           *( gx*vel(1,iel)                                    &
                            + gy*vel(2,iel)                                    &
                            + gz*vel(3,iel) )
enddo

!                                  Kij*Sij           LAMBDA   Cp      MUT
!     FACE DIFFUSION "VELOCITY" : --------- avec K = ------ + -- .------------
!     =========================    IJ.nij              Cv     Cv  TURB_SCHMIDT

! Only SGDH available

allocate(w1(ncelet))
allocate(viscf(nfac))
allocate(viscb(nfabor))

if (vcopt_e%idiff.ge. 1) then

  call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

!     MUT/TURB_SCHMIDT
  do iel = 1, ncel
    w1(iel) = visct(iel)/turb_schmidt
  enddo
!     CP*MUT/TURB_SCHMIDT
  if(icp.ge.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)*cpro_cp(iel)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)*cp0
    enddo
  endif
!     (CP/CV)*MUT/TURB_SCHMIDT
  if(icv.ge.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)/cpro_cv(iel)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)/cv0
    enddo
  endif
!     (CP/CV)*MUT/TURB_SCHMIDT+LAMBDA/CV
  if(ifcvsl.lt.0)then
    do iel = 1, ncel
      w1(iel) = w1(iel) + visls0(iscal)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel) + cpro_viscls(iel)
    enddo
  endif

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )


!     COMPLEMENTARY DIFFUSIVE TERM : - div( K grad ( epsilon - Cv.T ) )
!     ============================                   1  2
!                                    - div( K grad ( -.u  ) )
!                                                    2

! Complementary term at cell centers

  ! Compute e - CvT

  ! At cell centers
  call cs_cf_thermo_eps_sup(crom, w9, ncel)
  !========================

  ! At boundary faces centers
  call cs_cf_thermo_eps_sup(brom, wb, nfabor)
  !========================

! Divergence computation with reconstruction


! Computation of the gradient of (0.5*u*u+EPSILONsup)

  do iel = 1, ncel
    w7(iel) = 0.5d0*( vel(1,iel)**2                                            &
                     +vel(2,iel)**2                                            &
                     +vel(3,iel)**2 ) + w9(iel)
  enddo

! Note : by default, since the parameters are unknowns, the velocity parameters
! are taken

  iii = iu
  inc = 1
  iccocg = 1
  imrgrp = vcopt_u%imrgra
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_u%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr
  extrap = vcopt_u%extrag

! Allocate temporary arrays
  allocate(coefap(nfabor))
  allocate(coefbp(nfabor))

  do ifac = 1, nfabor
    coefap(ifac) = zero
    coefbp(ifac) = 1.d0
  enddo

!  f_id0 = -1 (indicates, for the rotation periodicity,
!              that the variable is not Rij)
  f_id0 = -1
  call gradient_s                                                   &
  !==========
   ( f_id0  , imrgrp , inc    , iccocg , nswrgp , imligp ,          &
     iwarnp , epsrgp , climgp , extrap ,                            &
     w7     , coefap , coefbp ,                                     &
     grad   )

! Free memory
  deallocate(coefap, coefbp)

! Internal faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    pnd   = pond(ifac)

! Computation of II' and JJ'

    diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
    diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
    diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
    djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) +  pnd  * dijpfx
    djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) +  pnd  * dijpfy
    djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) +  pnd  * dijpfz

    pip = w7(ii) + grad(1,ii)*diipfx+grad(2,ii)*diipfy+grad(3,ii)*diipfz
    pjp = w7(jj) + grad(1,jj)*djjpfx+grad(2,jj)*djjpfy+grad(3,jj)*djjpfz

    flux = viscf(ifac)*(pip-pjp)

    smbrs(ii) = smbrs(ii) + flux
    smbrs(jj) = smbrs(jj) - flux

  enddo

  if (ippmod(igmix).gt.0) then

    ! Diffusion flux for the species at internal faces

    allocate(kspe(ncelet),viscfk(nfac),viscbk(nfabor))

    ! Diffusion coefficient  T*lambda*Cvk/Cv
    do iel =1, ncel
      kspe(iel) = w1(iel)* cvar_tempk(iel)
    enddo

    call viscfa(imvisf, kspe, viscfk, viscbk)

    deallocate(kspe)

    allocate(grad_dd(nfac))

    do ifac = 1, nfac
      grad_dd(ifac) = 0.d0
    enddo

    do iesp = 1, nscasp

      ivarsp = isca(iscasp(iesp))
      call field_get_val_s(ivarfl(ivarsp), cvar_yk)
      call field_get_key_struct_gas_mix_species_prop(ivarfl(ivarsp), s_k)

      mk =  s_k%mol_mas
      cpk = s_k%cp
      cvk = cpk - cs_physical_constants_r/mk

      use_previous = 0
      call field_gradient_scalar(ivarfl(ivarsp), use_previous, 0, inc, &
                                 iccocg, grad)

      do ifac = 1, nfac

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        dijpfx = dijpf(1,ifac)
        dijpfy = dijpf(2,ifac)
        dijpfz = dijpf(3,ifac)

        pnd   = pond(ifac)

        ! Computation of II' and JJ'
        diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
        diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
        diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
        djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) +  pnd  * dijpfx
        djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) +  pnd  * dijpfy
        djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) +  pnd  * dijpfz

        yip = cvar_yk(ii) + grad(1,ii)*diipfx &
                          + grad(2,ii)*diipfy &
                          + grad(3,ii)*diipfz
        yjp = cvar_yk(jj) + grad(1,jj)*djjpfx &
                          + grad(2,jj)*djjpfy &
                          + grad(3,jj)*djjpfz

        ! Gradient of deduced species
        grad_dd(ifac) = grad_dd(ifac)-(yjp-yip)

        flux = viscfk(ifac)*cvk*(yip-yjp)

        smbrs(ii) = smbrs(ii) + flux
        smbrs(jj) = smbrs(jj) - flux

      enddo

    enddo

    ! Diffusion flux for the deduced species

    call field_get_key_struct_gas_mix_species_prop(iddgas, s_k)
    mk =  s_k%mol_mas
    cpk = s_k%cp
    cvk = cpk - cs_physical_constants_r/mk

    do ifac = 1, nfac

      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)

      flux = viscf(ifac)*grad_dd(ifac)*cvk

      smbrs(ii) = smbrs(ii) + flux
      smbrs(jj) = smbrs(jj) - flux

    enddo

    deallocate(grad_dd)

  endif ! ippmod(igmix)

  ! Assembling based on boundary faces
  ! for the faces where a flux or a temperature is imposed,
  ! all is taken into account by the energy diffusion term.
  ! Hence the contribution of the terms in u2 and e-CvT shouldn't be taken into
  ! account when ifbet(ifac).ne.0

  call field_get_coefa_v(ivarfl(iu), coefau)
  call field_get_coefb_v(ivarfl(iu), coefbu)

  do ifac = 1, nfabor

    if (ifbet(ifac).eq.0) then

      iel = ifabor(ifac)

      flux = viscb(ifac)*(w1(iel)/distb(ifac))*                  &
            ( w9(iel) - wb(ifac)                                 &
            + 0.5d0*( vel(1,iel)**2 -                            &
              ( coefau(1, ifac) + coefbu(1, 1, ifac)*vel(1,iel)  &
                                + coefbu(1, 2, ifac)*vel(2,iel)  &
                                + coefbu(1, 3, ifac)*vel(3,iel)  &
              )**2                                               &
                   + vel(2,iel)**2 -                             &
              ( coefau(2, ifac) + coefbu(2, 1, ifac)*vel(1,iel)  &
                                + coefbu(2, 2, ifac)*vel(2,iel)  &
                                + coefbu(2, 3, ifac)*vel(3,iel)  &
              )**2                                               &
                   + vel(3,iel)**2 -                             &
              ( coefau(3, ifac) + coefbu(3, 1, ifac)*vel(1,iel)  &
                                + coefbu(3, 2, ifac)*vel(2,iel)  &
                                + coefbu(3, 3, ifac)*vel(3,iel)  &
              )**2                                               &
            ))

      smbrs(iel) = smbrs(iel) + flux

    endif

  enddo

  if (ippmod(igmix).gt.0) then

    call field_get_coefa_s(ivarfl(isca(itempk)), coefat)
    call field_get_coefb_s(ivarfl(isca(itempk)), coefbt)

    allocate(grad_dd(nfabor), btemp(nfabor))

    use_previous = 0
    call field_gradient_scalar(ivarfl(isca(itempk)), use_previous, 0, inc,&
                               iccocg, grad)

    do ifac = 1, nfabor
      grad_dd(ifac) = 0.d0

      tip = cvar_tempk(iel) + grad(1,iel)*diipb(1,ifac) &
                            + grad(2,iel)*diipb(2,ifac) &
                            + grad(3,iel)*diipb(3,ifac)
      btemp(ifac) = coefat(ifac)+coefbt(ifac)*tip
    enddo

    do iesp = 1, nscasp
      ivarsp = isca(iscasp(iesp))
      call field_get_coefa_s(ivarfl(ivarsp), coefayk)
      call field_get_coefb_s(ivarfl(ivarsp), coefbyk)
      call field_get_val_s(ivarfl(ivarsp), cvar_yk)
      call field_get_key_struct_gas_mix_species_prop(ivarfl(ivarsp), s_k)

      mk =  s_k%mol_mas
      cpk = s_k%cp
      cvk = cpk - cs_physical_constants_r/mk

      use_previous = 0
      call field_gradient_scalar(ivarfl(ivarsp), use_previous, 0, inc, &
                                 iccocg, grad)

      do ifac = 1, nfabor
        if (ifbet(ifac).eq.0) then
          iel = ifabor(ifac)

          yip = cvar_yk(iel) + grad(1,iel)*diipb(1,ifac) &
                             + grad(2,iel)*diipb(2,ifac) &
                             + grad(3,iel)*diipb(3,ifac)

          gradnb = coefayk(ifac)+(coefbyk(ifac)-1)*yip

          grad_dd(ifac) =  grad_dd(ifac) - gradnb

          flux = viscbk(ifac)*w1(iel)*btemp(ifac)*cvk/distb(ifac)*(-gradnb)

          smbrs(iel) = smbrs(iel) + flux
        endif
      enddo ! end ifac loop
    enddo ! end iesp loop

    ! Boundary diffusion flux for the deduced species
    call field_get_key_struct_gas_mix_species_prop(iddgas, s_k)

    mk =  s_k%mol_mas
    cpk = s_k%cp
    cvk = cpk - cs_physical_constants_r/mk

    do ifac = 1, nfabor
      if (ifbet(ifac).eq.0) then
        iel = ifabor(ifac)

        flux = viscbk(ifac)*w1(iel)*btemp(ifac)*cvk/distb(ifac)*grad_dd(ifac)

        smbrs(iel) = smbrs(iel) + flux
      endif
    enddo

    deallocate(grad_dd, btemp)
    deallocate(viscfk, viscbk)

  endif ! ippmod(igmix)

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

!===============================================================================
! 4. Solving
!===============================================================================

iconvp = vcopt_e%iconv
idiffp = vcopt_e%idiff
ndircp = vcopt_e%ndircl
imrgrp = vcopt_e%imrgra
nswrsp = vcopt_e%nswrsm
nswrgp = vcopt_e%nswrgr
imligp = vcopt_e%imligr
ircflp = vcopt_e%ircflu
ischcp = vcopt_e%ischcv
isstpp = vcopt_e%isstpc
iwarnp = vcopt_e%iwarni
blencp = vcopt_e%blencv
epsilp = vcopt_e%epsilo
epsrsp = vcopt_e%epsrsm
epsrgp = vcopt_e%epsrgr
climgp = vcopt_e%climgr
extrap = vcopt_e%extrag
relaxp = vcopt_e%relaxv
thetap = vcopt_e%thetav
iescap = 0

!  idtvar = 1  => unsteady
imucpp = 0  ! not a thermal scalar
idftnp = ISOTROPIC_DIFFUSION
iswdyp = 0  ! no dynamic relaxation

! impose boundary convective at some faces (face indicator icvfli)
icvflb = 1
normp = -1.d0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits                                                      &
!==========
( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
  imrgrp , nswrsp , nswrgp , imligp , ircflp ,                   &
  ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
  iwarnp , normp  ,                                              &
  blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
  relaxp , thetap ,                                              &
  cvara_energ     , cvara_energ     ,                            &
  coefap , coefbp , cofafp , cofbfp ,                            &
  imasfl, bmasfl,                                                &
  viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
  rvoid  , rvoid  ,                                              &
  icvflb , icvfli ,                                              &
  rovsdt , smbrs  , cvar_energ      , dpvar  ,                   &
  rvoid  , rvoid  )

deallocate(viscf)
deallocate(viscb)

!===============================================================================
! 5. Printings and clippings
!===============================================================================

call clpsca(iscal)

! --- Traitement utilisateur pour gestion plus fine des bornes
!       et actions correctives eventuelles.
call cs_cf_check_internal_energy(cvar_energ, ncel, vel)

! Explicit balance (see codits : the increment is removed)

if (vcopt_e%iwarni.ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                    &
            - vcopt_e%istat*(crom(iel)/dt(iel))*cell_f_vol(iel)                  &
                *(cvar_energ(iel)-cvara_energ(iel))                            &
                * max(0,min(vcopt_e%nswrsm-2,1))
  enddo
  sclnor = sqrt(cs_gdot(ncel,smbrs,smbrs))
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. Final updating of the pressure (and temperature)
!===============================================================================
!                             n+1      n+1  n+1
! The state equation is used P   =P(RHO   ,H   )

! Computation of P and T at cell centers
call cs_cf_thermo_pt_from_de(cpro_cp, cpro_cv, crom, cvar_energ, cvar_pr, &
                             cvar_tempk, vel, &
                             cvar_fracv, cvar_fracm, cvar_frace, ncel)

33 continue
! Barotropic version
if (ippmod(icompf).eq.1) then
  do iel = 1, ncel
    cvar_energ(iel) = eint0
  enddo
endif
!                             n+1      n+1  n+1
! The state equation is used P   =P(rho   ,e   )

!===============================================================================
! 7. Communication of pressure, energy and temperature
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(cvar_pr)
  !==========
  call synsca(cvar_energ)
  !==========
  call synsca(cvar_tempk)
  !==========
endif

! Free memory
if (allocated(wb)) deallocate(wb)
if (allocated(smbrs)) deallocate(smbrs, rovsdt)
if (allocated(grad)) deallocate(grad)
if (allocated(w1)) deallocate(w1)
if (allocated(w7)) deallocate(w7, w9)

!--------
! Formats
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION FOR THE VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)


!----
! End
!----

return

end subroutine cfener

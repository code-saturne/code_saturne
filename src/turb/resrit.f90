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

!> \file resrit.f90
!>
!> \brief This subroutine perform the solving of the transport equation
!> of the turbulent heat fluxes.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in]     iscal         number of the scalar used
!> \param[in]     xcpp          \f$ C_p \f$
!> \param[in,out] xut, xuta     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     dt            time step (per cell)
!> \param[in]     gradv         mean velocity gradient
!> \param[in]     gradt         mean scalar gradient
!> \param[in]     grad_al       alpha scalar gradient
!______________________________________________________________________________!

subroutine resrit &
 ( nscal  ,                                                       &
   iscal  , xcpp   , xut    , xuta   ,                            &
   dt     ,                                                       &
   gradv  , gradt  , grad_al)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstnum
use cstphy
use parall
use period
use field
use mesh
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nscal , iscal

double precision dt(ncelet)
double precision xcpp(ncelet), xut(3,ncelet), xuta(3,ncelet)
double precision gradv(3,3,ncelet)
double precision gradt(3,ncelet)
double precision grad_al(3,ncelet)

! Local variables

integer          iel
integer          ii, ivar
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          st_prv_id
integer          ivisep, ifcvsl
integer          isou, jsou
integer          itt
integer          idftnp, iswdyp, icvflb
integer          f_id
integer          keyvar, iut

integer          ivoid(1)

double precision blencp, epsilp, epsrgp, climgp, relaxp
double precision epsrsp
double precision trrij
double precision thets , thetv , thetp1
double precision xttke , prdtl
double precision grav(3)
double precision xrij(3,3),phiith(3), phiitw(3)
double precision xnal(3), xnoral
double precision alpha, xttdrbt, xttdrbw
double precision pk, gk, xxc1, xxc2, xxc3, imp_term
double precision rctse

double precision rvoid(1)

character(len=80) :: fname, name

double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:) :: viscb
double precision, allocatable, dimension(:) :: viscf
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: smbrut
double precision, allocatable, dimension(:,:,:) :: fimp
double precision, allocatable, dimension(:) :: w1

double precision, dimension(:,:), pointer :: coefav, cofafv, visten
double precision, dimension(:,:,:), pointer :: coefbv, cofbfv
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, cpro_beta
double precision, dimension(:), pointer :: cvar_al
double precision, dimension(:), pointer :: cvar_ep
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: cvar_tt, cvara_tt
double precision, dimension(:), pointer :: viscl, visct, viscls, c_st_prv

type(var_cal_opt) :: vcopt, vcopt_ut

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet))
allocate(viscce(6,ncelet))
allocate(smbrut(3,ncelet))
allocate(fimp(3,3,ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

if ((itytur.eq.2).or.(itytur.eq.5).or.(iturb.eq.60)) then
  write(nfecra,*)'Utiliser un modele Rij avec ces modeles de thermiques'!FIXME
  call csexit(1)
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

call field_get_val_s(ivarfl(iep), cvar_ep)

if (irijco.eq.1) then
  call field_get_val_v(ivarfl(irij), cvar_rij)
else
  call field_get_val_s(ivarfl(ir11), cvar_r11)
  call field_get_val_s(ivarfl(ir22), cvar_r22)
  call field_get_val_s(ivarfl(ir33), cvar_r33)
  call field_get_val_s(ivarfl(ir12), cvar_r12)
  call field_get_val_s(ivarfl(ir13), cvar_r13)
  call field_get_val_s(ivarfl(ir23), cvar_r23)
endif

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

if (ibeta.ge.0) then
  call field_get_val_s(ibeta, cpro_beta)
endif

call field_get_val_v(ivsten, visten)

ivar = isca(iscal)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

call field_get_key_id("variable_id", keyvar)
! Get name of the reference scalar
call field_get_name(ivarfl(ivar), fname)
! Index of the corresponding turbulent flux
call field_get_id(trim(fname)//'_turbulent_flux', f_id)
! Set pointer values of turbulent fluxes
call field_get_key_int(f_id, keyvar, iut)

call field_get_key_struct_var_cal_opt(ivarfl(iut), vcopt_ut)

if (vcopt%iwarni.ge.1) then
  call field_get_name(ivarfl(ivar), name)
  write(nfecra,1000) trim(name)//'_turbulent_flux'
endif

! S pour Source, V pour Variable
thets  = thetst
thetv  = vcopt%thetav

call field_get_key_int(ivarfl(ivar), kstprv, st_prv_id)
if (st_prv_id.ge.0) then
  call field_get_val_s(st_prv_id, c_st_prv)
else
  c_st_prv=> null()
endif

call field_get_key_int (ivarfl(ivar), kivisl, ifcvsl)
if (ifcvsl .ge. 0) then
  call field_get_val_s(ifcvsl, viscls)
endif

do iel = 1, ncelet
  do isou = 1, 3
    smbrut(isou,iel) = 0.d0
    do jsou = 1, 3
      fimp(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

! Find the corresponding variance of the scalar iscal
itt = -1
if ((abs(gx)+abs(gy)+abs(gz)).gt.0) then
  grav(1) = gx
  grav(2) = gy
  grav(3) = gz
  do ii = 1, nscal
    if (iscavr(ii).eq.iscal) itt = ii
  enddo
endif

!===============================================================================
! 2. Mass source terms FIXME
!===============================================================================

if (st_prv_id.ge.0) then
  do iel = 1, ncel
    do isou = 1,3
      smbrut(isou,iel) = fimp(isou,isou,iel)*xuta(isou,iel)
      fimp(isou,isou,iel) = - thetv*fimp(isou,isou,iel)
    enddo
  enddo
! If we do not extrapolate the source terms
else
  do iel = 1, ncel
    do isou = 1, 3
      ! User source term
      smbrut(isou,iel) = smbrut(isou,iel) + fimp(isou,isou,iel)*xuta(isou,iel)
      ! Diagonal
      fimp(isou,isou,iel) = max(-fimp(isou,isou,iel),zero)
    enddo
  enddo
endif

!===============================================================================
! 3. Unsteady term
!===============================================================================

do iel = 1, ncel
  do isou = 1, 3
    fimp(isou,isou,iel) = fimp(isou,isou,iel)                                  &
                        + vcopt%istat*(crom(iel)/dt(iel))*volume(iel)
  enddo
enddo

!===============================================================================
! 4. Right Hand Side of the thermal fluxes:
!     rho*(Pit + Git + Phi*_it - eps_it)
!===============================================================================

if (itt.gt.0) then
  call field_get_val_s(ivarfl(isca(itt)), cvar_tt)
  call field_get_val_prev_s(ivarfl(isca(itt)), cvara_tt)
endif
if (iturt(iscal).eq.31) then
  ! Name of the scalar
  call field_get_name(ivarfl(isca(iscal)), fname)

  ! Index of the corresponding turbulent flux
  call field_get_id(trim(fname)//'_alpha', f_id)

  call field_get_val_s(f_id, cvar_al)
endif

do iel = 1, ncel

  if (irijco.eq.1) then
    xrij(1,1) = cvar_rij(1,iel)
    xrij(2,2) = cvar_rij(2,iel)
    xrij(3,3) = cvar_rij(3,iel)
    xrij(1,2) = cvar_rij(4,iel)
    xrij(2,3) = cvar_rij(5,iel)
    xrij(1,3) = cvar_rij(6,iel)
  else
    xrij(1,1) = cvar_r11(iel)
    xrij(2,2) = cvar_r22(iel)
    xrij(3,3) = cvar_r33(iel)
    xrij(1,2) = cvar_r12(iel)
    xrij(1,3) = cvar_r13(iel)
    xrij(2,3) = cvar_r23(iel)
  endif
  xrij(2,1) = xrij(1,2)
  xrij(3,1) = xrij(1,3)
  xrij(3,2) = xrij(2,3)

  if (ifcvsl.ge.0) then
    prdtl = viscl(iel)*xcpp(iel)/viscls(iel)
  else
    prdtl = viscl(iel)*xcpp(iel)/visls0(iscal)
  endif

  trrij  = 0.5d0*(xrij(1,1) + xrij(2,2) + xrij(3,3))
  ! --- Compute Durbin time scheme
  xttke  = trrij/cvar_ep(iel)

  if (iturt(iscal).eq.31) then
    alpha = cvar_al(iel)
    !FIXME Warning / rhebdfm**0.5 compared to F Dehoux
    xttdrbt = xttke * sqrt( (1.d0-alpha)*prdtl/ rhebdfm + alpha )
    xttdrbw = xttdrbt * sqrt(rhebdfm/prdtl) ! And so multiplied by (R/Prandt)^0.5

    ! Compute the unite normal vector
    xnoral = &
      ( grad_al(1, iel)*grad_al(1, iel) &
      + grad_al(2, iel)*grad_al(2, iel) &
      + grad_al(3, iel)*grad_al(3, iel))
    xnoral = sqrt(xnoral)

    if (xnoral.le.epzero/cell_f_vol(iel)**(1.d0/3.d0)) then
      xnal(1) = 0.d0
      xnal(2) = 0.d0
      xnal(3) = 0.d0
    else
      xnal(1) = grad_al(1, iel)/xnoral
      xnal(2) = grad_al(2, iel)/xnoral
      xnal(3) = grad_al(3, iel)/xnoral
    endif

    ! Production and buoyancy for TKE
    pk = -( xrij(1,1)*gradv(1,1,iel) +&
            xrij(1,2)*gradv(1,2,iel) +&
            xrij(1,3)*gradv(1,3,iel) +&
            xrij(2,1)*gradv(2,1,iel) +&
            xrij(2,2)*gradv(2,2,iel) +&
            xrij(2,3)*gradv(2,3,iel) +&
            xrij(3,1)*gradv(3,1,iel) +&
            xrij(3,2)*gradv(3,2,iel) +&
            xrij(3,3)*gradv(3,3,iel) )

    if (ibeta.ge.0) then
      gk = cpro_beta(iel)*(xuta(1,iel)*gx + xuta(2,iel)*gy + xuta(3,iel)*gz)!FIXME make buoyant term coherent elsewhere
    else
      gk = 0.d0
    endif
    xxc1 = 1.d0+2.d0*(1.d0 - cvar_al(iel))*(pk+gk)/cvar_ep(iel)
    xxc2 = 0.5d0*(1.d0+1.d0/prdtl)*(1.d0-0.3d0*(1.d0 - cvar_al(iel)) &
      *(pk+gk)/cvar_ep(iel))
    xxc3 = xxc2

  else
    xttdrbt = xttke
    xttdrbw = xttke
    alpha = 1.d0
    xxc1 = 0.d0
    xxc2 = 0.d0
    xxc3 = 0.d0
    xnal(1) = 0.d0
    xnal(2) = 0.d0
    xnal(3) = 0.d0
  endif

  do isou = 1, 3
    phiith(isou) = -c1trit / xttdrbt * xuta(isou,iel)               &
                 + c2trit*(xuta(1,iel)*gradv(1,isou,iel)            &
                          +xuta(2,iel)*gradv(2,isou,iel)            &
                          +xuta(3,iel)*gradv(3,isou,iel))           &
                 + c4trit*(-xrij(isou,1)*gradt(1,iel)               &
                           -xrij(isou,2)*gradt(2,iel)               &
                           -xrij(isou,3)*gradt(3,iel))
    if (itt.gt.0.and.ibeta.ge.0) then
      phiith(isou) = phiith(isou)                                              &
             + c3trit*(cpro_beta(iel)*grav(isou)*cvar_tt(iel))
    endif

    phiitw(isou) = -1.d0/xttdrbw * xxc1 *             & !FIXME full implicit
                 ( xuta(1,iel)*xnal(1)*xnal(isou)     &
                 + xuta(2,iel)*xnal(2)*xnal(isou)     &
                 + xuta(3,iel)*xnal(3)*xnal(isou))

    ! Pressure/thermal fluctuation correlation term
    !----------------------------------------------
    smbrut(isou,iel) = smbrut(isou,iel) +                                      &
                volume(iel)*crom(iel)*(  alpha *         phiith(isou)          &
                                      + (1.d0 - alpha) * phiitw(isou))

    imp_term = max(volume(iel)*crom(iel)*(                                     &
              alpha        * (c1trit/xttdrbt-c2trit*gradv(isou,isou,iel))      &
            + (1.d0-alpha) * (xxc1*xnal(isou)*xnal(isou) / xttdrbw) &
            ), 0.d0)

    fimp(isou,isou,iel) = fimp(isou,isou,iel) + imp_term

    ! Production terms
    !-----------------
    smbrut(isou,iel) = smbrut(isou,iel)                              &
                     + volume(iel)*crom(iel)                         &
                       ! Production term due to the mean velcoity
                       *( -xuta(1,iel)*gradv(1,isou,iel)             &
                          -xuta(2,iel)*gradv(2,isou,iel)             &
                          -xuta(3,iel)*gradv(3,isou,iel)             &
                       ! Production term due to the mean temperature
                         -xrij(isou,1)*gradt(1,iel)                  &
                         -xrij(isou,2)*gradt(2,iel)                  &
                         -xrij(isou,3)*gradt(3,iel)                  &
                        )

    ! Production term due to the gravity
    if (itt.gt.0.and.ibeta.ge.0) then
      smbrut(isou,iel) = smbrut(isou,iel)                            &
                       + volume(iel)*crom(iel)*(            &
               -grav(isou)*cpro_beta(iel)*cvara_tt(iel))
    endif

    ! Dissipation (Wall term only because "h" term is zero
    smbrut(isou,iel) = smbrut(isou,iel) -                                      &
                volume(iel)*crom(iel)*(1.d0 - alpha) / xttdrbw *               &
                       ( xxc2 * xuta(isou, iel)                                &
                       + xxc3*( xuta(1,iel)*xnal(1)*xnal(isou)                 &
                              + xuta(2,iel)*xnal(2)*xnal(isou)                 &
                              + xuta(3,iel)*xnal(3)*xnal(isou)))

    ! TODO we can implicite more terms
    imp_term = max(volume(iel)*crom(iel)*(1.d0 - alpha) / xttdrbw *            &
                  ( xxc2 + xxc3 * xnal(isou)*xnal(isou)), 0.d0)
    fimp(isou,isou,iel) = fimp(isou,isou,iel) + imp_term

  enddo
enddo

!===============================================================================
! 5. Tensorial diffusion
!===============================================================================
! Symmetric tensor diffusivity (GGDH)
if (vcopt_ut%idften.eq.6) then
  do iel = 1, ncel

    if (ifcvsl.ge.0) then
      prdtl = viscl(iel)*xcpp(iel)/viscls(iel)
    else
      prdtl = viscl(iel)*xcpp(iel)/visls0(iscal)
    endif

    do isou = 1, 6
      if (isou.le.3) then
        viscce(isou,iel) = 0.5d0*(viscl(iel)*(1.d0+1.d0/prdtl))    &
                         + ctheta(iscal)*visten(isou,iel)/csrij
      else
        viscce(isou,iel) = ctheta(iscal)*visten(isou,iel)/csrij
      endif
    enddo
  enddo

  iwarnp = vcopt%iwarni

  call vitens &
  ( viscce , iwarnp ,             &
    weighf , weighb ,             &
    viscf  , viscb  )

! Scalar diffusivity
else

  do iel = 1, ncel
    if (irijco.eq.1) then
      trrij = 0.5d0 * (cvar_rij(1,iel) + cvar_rij(2,iel) + cvar_rij(3,iel))
    else
      trrij = 0.5d0 * (cvar_r11(iel) + cvar_r22(iel) + cvar_r33(iel))
    end if
    rctse = crom(iel) * csrij * trrij**2 / cvar_ep(iel)
    w1(iel) = viscl(iel) + vcopt%idifft*rctse
  enddo

  call viscfa                    &
  ( imvisf ,                     &
   w1     ,                      &
   viscf  , viscb  )

end if

!===============================================================================
! 6. Vectorial solving of the turbulent thermal fluxes
!===============================================================================

if (st_prv_id.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    do isou = 1, 3
      smbrut(isou,iel) = smbrut(isou,iel) + thetp1*c_st_prv(iel) !FIXME
    enddo
  enddo
endif

! Name of the scalar ivar
call field_get_name(ivarfl(ivar), fname)

! Index of the corresponding turbulent flux
call field_get_id(trim(fname)//'_turbulent_flux', f_id)

call field_get_coefa_v(f_id,coefav)
call field_get_coefb_v(f_id,coefbv)
call field_get_coefaf_v(f_id,cofafv)
call field_get_coefbf_v(f_id,cofbfv)

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
idftnp = 6
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
relaxp = vcopt%relaxv

! We do not take into account transpose of grad
ivisep = 0

! all boundary convective flux with upwind
icvflb = 0

call coditv &
(idtvar , f_id   , iconvp , idiffp , ndircp ,                   &
 imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
 ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
 iwarnp ,                                                       &
 blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
 relaxp , thetv  ,                                              &
 xuta   , xuta   ,                                              &
 coefav , coefbv , cofafv , cofbfv ,                            &
 imasfl , bmasfl ,                                              &
 viscf  , viscb  , viscf  , viscb  , rvoid  , rvoid  ,          &
 viscce , weighf , weighb ,                                     &
 icvflb , ivoid  ,                                              &
 fimp   ,                                                       &
 smbrut ,                                                       &
 xut    ,                                                       &
 rvoid  )

!===============================================================================
! 7. Writings
!===============================================================================

! Free memory
deallocate(viscce)
deallocate(viscf, viscb)
deallocate(smbrut)
deallocate(fimp)
deallocate(weighf, weighb)
deallocate(w1)
!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,'           Resolution pour la variable ',A23,/)

#else

 1000 format(/,'           Solving variable ',A23           ,/)
#endif

!----
! End
!----

return
end subroutine

!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file divrit.f90
!>
!> \brief This subroutine perform  add the divergence of turbulent flux
!> to the transport equation of a scalar.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nscal         total number of scalars
!> \param[in]     iscal         number of the scalar used
!> \param[in]     dt            time step (per cell)
!> \param[in]     xcpp          Cp
!> \param[out]    smbrs         Right hand side to update
!_______________________________________________________________________________

subroutine divrit &
 ( nscal  , iscal  ,                                              &
   dt     ,                                                       &
   xcpp   ,                                                       &
   vistet ,                                                       &
   smbrs )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum!
use pointe
use field
use field_operator
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nscal
integer          iscal
double precision dt(ncelet)
double precision xcpp(ncelet)
double precision vistet(6,ncelet)
double precision smbrs(ncelet)

! Local variables

integer          ifac, init, inc, iprev
integer          iccocg,iflmb0
integer          nswrgp, imligp, iwarnp
integer          itypfl
integer          ivar , iel, ii, jj
integer          itt
integer          f_id, f_id0, f_id_al
integer          ifcvsl

double precision epsrgp, climgp, extrap
double precision xk, xe, xtt
double precision grav(3),xrij(3,3), temp(3)
double precision xnal(3), xnoral
double precision alpha_theta, xR_h, xR, prdtl
double precision xpk, xgk
double precision eta_ebafm, gamma_ebafm, xi_ebafm, gamma_ebggdh
double precision xxc1, xxc2, xxc3
double precision coeff_imp

character(len=80) :: fname

double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: gradt, grad_al
double precision, allocatable, dimension(:,:) :: coefat
double precision, allocatable, dimension(:,:,:) :: coefbt
double precision, allocatable, dimension(:) :: thflxf, thflxb
double precision, allocatable, dimension(:) :: divut
double precision, allocatable, dimension(:,:) :: w1

double precision, dimension(:,:), pointer :: cofarut
double precision, dimension(:,:,:), pointer :: cofbrut
double precision, dimension(:,:), pointer :: xut
double precision, dimension(:,:), pointer :: xuta, cvara_rij
double precision, dimension(:), pointer :: brom, crom, cpro_beta
double precision, dimension(:), pointer :: viscl, viscls
double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvara_r12, cvara_r13, cvara_r23
double precision, dimension(:), pointer :: cvara_tt
double precision, dimension(:), pointer :: cvar_al

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

f_id0 = -1

! Initializations to avoid compiler warnings
xtt = 0.d0

! First component is for x,y,z  and the 2nd for u,v,w
allocate(gradv(3,3,ncelet))
allocate(gradt(3,ncelet), thflxf(nfac), thflxb(nfabor))

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

if (ibeta.ge.0) then
  call field_get_val_s(ibeta, cpro_beta)!FIXME make it dependant on the scalar
endif

call field_get_val_s(iviscl, viscl)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, viscls)
end if

! Compute scalar gradient
ivar = isca(iscal)

iprev = 1
iccocg = 1
inc = 1

call field_gradient_scalar( &
  ivarfl(ivar)    , iprev, imrgra, inc    , &
  iccocg ,                                  &
  gradt  )

! Name of the scalar ivar
call field_get_name(ivarfl(ivar), fname)

! Index of the corresponding turbulent flux
call field_get_id(trim(fname)//'_turbulent_flux', f_id)

call field_get_val_v(f_id, xut)

! EB- AFM or EB-DFM: compute the gradient of alpha of the scalar
if (iturt(iscal).eq.11 .or. iturt(iscal).eq.21 .or. iturt(iscal).eq.31) then
  ! Index of the corresponding alpha
  call field_get_id(trim(fname)//'_alpha', f_id_al)
  call field_get_val_s(f_id_al, cvar_al)

  iprev = 0
  iccocg = 1
  inc = 1

  allocate(grad_al(3,ncelet))

  call field_gradient_scalar( &
    f_id_al , iprev, imrgra, inc    , &
    iccocg ,                          &
    grad_al)

endif

! Compute velocity gradient
iprev  = 0
inc    = 1

! WARNING: gradv(xyz, uvw, iel)

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,  &
                           gradv)

! Find the variance of the thermal scalar
itt = -1
if (((abs(gx)+abs(gy)+abs(gz)).gt.epzero).and.(irovar.gt.0.or.idilat.eq.0).and.   &
    ((ityturt(iscal).eq.2).or.(ityturt(iscal).eq.3))) then
  grav(1) = gx
  grav(2) = gy
  grav(3) = gz
  do ii = 1, nscal
    if (iscavr(ii).eq.iscalt) itt = ii
  enddo
  if (itt.le.0) then
    write(nfecra,9999)
    call csexit(1)
  endif
endif

!===============================================================================
! 2. Agebraic models AFM
!===============================================================================
if (ityturt(iscal).ne.3) then

  call field_get_val_prev_s(ivarfl(iep), cvara_ep)

  if (irijco.eq.1) then
    call field_get_val_prev_v(ivarfl(irij), cvara_rij)
  else
    call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
    call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
    call field_get_val_prev_s(ivarfl(ir33), cvara_r33)
    call field_get_val_prev_s(ivarfl(ir12), cvara_r12)
    call field_get_val_prev_s(ivarfl(ir13), cvara_r13)
    call field_get_val_prev_s(ivarfl(ir23), cvara_r23)
  endif

  allocate(w1(3,ncelet))

  do ifac = 1, nfac
    thflxf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    thflxb(ifac) = 0.d0
  enddo

  if (itt.gt.0) call field_get_val_prev_s(ivarfl(isca(itt)), cvara_tt)

  do iel = 1, ncel
    !Rij
    ! Coupled version
    if(irijco.eq.1) then
      xrij(1,1) = cvara_rij(1,iel)
      xrij(2,2) = cvara_rij(2,iel)
      xrij(3,3) = cvara_rij(3,iel)
      xrij(1,2) = cvara_rij(4,iel)
      xrij(2,3) = cvara_rij(5,iel)
      xrij(1,3) = cvara_rij(6,iel)
      xrij(2,1) = xrij(1,2)
      xrij(3,1) = xrij(1,3)
      xrij(3,2) = xrij(2,3)
    else
    ! Uncoupled version
      xrij(1,1) = cvara_r11(iel)
      xrij(2,2) = cvara_r22(iel)
      xrij(3,3) = cvara_r33(iel)
      xrij(1,2) = cvara_r12(iel)
      xrij(1,3) = cvara_r13(iel)
      xrij(2,3) = cvara_r23(iel)
      xrij(2,1) = xrij(1,2)
      xrij(3,1) = xrij(1,3)
      xrij(3,2) = xrij(2,3)
    endif
    ! Epsilon
    xe = cvara_ep(iel)
    ! Kinetic turbulent energy
    xk = 0.5d0*(xrij(1,1)+xrij(2,2)+xrij(3,3))
    !  Turbulent time-scale (constant in AFM)
    xtt = xk/xe

    if (iturt(iscal).eq.11.or.iturt(iscal).eq.21) then

      alpha_theta = cvar_al(iel)

      ! Computation of production and buoyancy
       xpk = -(xrij(1,1)*gradv(1,1,iel) +&
               xrij(1,2)*gradv(1,2,iel) +&
               xrij(1,3)*gradv(1,3,iel) +&
               xrij(2,1)*gradv(2,1,iel) +&
               xrij(2,2)*gradv(2,2,iel) +&
               xrij(2,3)*gradv(2,3,iel) +&
               xrij(3,1)*gradv(3,1,iel) +&
               xrij(3,2)*gradv(3,2,iel) +&
               xrij(3,3)*gradv(3,3,iel) )

       if (ibeta.ge.0) then
         xgk = cpro_beta(iel)*(xut(1,iel)*gx + xut(2,iel)*gy + xut(3,iel)*gz)
       else
         xgk = 0.d0
       endif

       ! Computation of the thermo-mecanical scales ratio R
       xR_h = 0.5d0
       if (ifcvsl.ge.0) then
         prdtl = viscl(iel)*xcpp(iel)/viscls(iel)
       else
         prdtl = viscl(iel)*xcpp(iel)/visls0(iscal)
       endif
       xR = ( 1.d0 - alpha_theta ) * prdtl + alpha_theta * xR_h

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
    end if

    ! Constants for EB-GGDH
    if (iturt(iscal).eq.11) then
      xxc1 = 1.d0+2.d0*(1.d0 - alpha_theta)*(xpk+xgk)/cvara_ep(iel)
      xxc2 = 0.5d0*(1.d0+1.d0/prdtl)*(1.d0-0.3d0*(1.d0 - alpha_theta) &
           *(xpk+xgk)/cvara_ep(iel))
      xxc3 = xxc2

      ctheta(iscal) = (0.97d0*xR**0.5d0)/( alpha_theta * (4.15d0*0.5d0**0.5d0)  &
                   +(1.d0-alpha_theta)*(prdtl**0.5d0)*xxc2)

      gamma_ebggdh = (1.d0-alpha_theta)*(xxc1 + xxc2)
    ! Constants for EB-AFM
    elseif (iturt(iscal).eq.21) then
      xxc1 = 1.d0+2.d0*(1.d0 - alpha_theta)*(xpk+xgk)/cvara_ep(iel)
      xxc2 = 0.5d0*(1.d0+1.d0/prdtl)*(1.d0-0.3d0*(1.d0 - alpha_theta) &
           *(xpk+xgk)/cvara_ep(iel))
      xxc3 = xxc2

      ctheta(iscal) = (0.97d0*xR**0.5d0)/( alpha_theta * (4.15d0*0.5d0**0.5d0)  &
                      +(1.d0-alpha_theta)*(prdtl**0.5d0)*xxc2)
      gamma_ebafm   = (1.d0-alpha_theta)*(xxc1 + xxc2)
      eta_ebafm     = 1.d0 - alpha_theta*0.6d0
      xi_ebafm      = 1.d0 - alpha_theta*0.3d0
    end if

    ! Compute thermal flux u'T'

    do ii = 1, 3
      temp(ii) = 0.d0

      ! AFM model
      !  "-C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2)"
      if (iturt(iscal).eq.20) then
        if (itt.gt.0.and.ibeta.ge.0) then
          temp(ii) = temp(ii) - ctheta(iscal)*xtt*                           &
                     etaafm*cpro_beta(iel)*grav(ii)*cvara_tt(iel)
        endif

        do jj = 1, 3
          ! Partial implicitation of "-C_theta*k/eps*( xi* uT'.Grad u )"
          ! Only the i.ne.j  components are added.
          if (ii.ne.jj) then
            temp(ii) = temp(ii)                                              &
                     - ctheta(iscal)*xtt*xiafm*gradv(jj,ii,iel)*xut(jj,iel)
          endif
        enddo
      endif

      ! EB-AFM model
      !  "-C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2 + eps/k gamma uT' ni nj)"
      if(iturt(iscal).eq.21) then
        if (itt.gt.0.and.ibeta.ge.0) then
          temp(ii) = temp(ii) - ctheta(iscal)*xtt*                           &
                     eta_ebafm*cpro_beta(iel)*grav(ii)*cvara_tt(iel)
        endif

        do jj = 1, 3
          ! Partial implicitation of "-C_theta*k/eps*( xi* uT'.Grad u + eps/k gamma uT' ni nj)"
          ! Only the i.ne.j  components are added.
          if (ii.ne.jj) then
            temp(ii) = temp(ii)                                                 &
                     - ctheta(iscal)*xtt*xi_ebafm*gradv(jj,ii,iel)*xut(jj,iel)  &
                     - ctheta(iscal)*gamma_ebafm*xnal(ii)*xnal(jj)*xut(jj,iel)
          endif
        enddo
      end if

      ! EB-GGDH model
      !  "-C_theta*k/eps*( eps/k gamma uT' ni nj)"
      if(iturt(iscal).eq.11) then
        do jj = 1, 3
          ! Partial implicitation of "-C_theta*k/eps*( eps/k gamma uT' ni nj)"
          ! Only the i.ne.j  components are added.
          if (ii.ne.jj) then
            temp(ii) = temp(ii)                                                 &
                     - ctheta(iscal)*gamma_ebggdh*xnal(ii)*xnal(jj)*xut(jj,iel)
          endif
        enddo
      end if

    enddo

    do ii = 1, 3
      ! Add the term in "grad T" which is implicited by the GGDH part in covofi.
      !  "-C_theta*k/eps* R.grad T"
      ! The resulting XUT array is only use for post processing purpose in
      ! (EB)GGDH & (EB)AFM
      xut(ii,iel) = temp(ii) - ctheta(iscal)*xtt*( xrij(ii,1)*gradt(1,iel)  &
                                                 + xrij(ii,2)*gradt(2,iel)  &
                                                 + xrij(ii,3)*gradt(3,iel))

      ! Partial implicitation of "-C_theta*k/eps*( xi* uT'.Grad u )" for
      ! EB-GGDH & (EB)-AFM
      ! X_i = C*Y_ij*X_j -> X_i = Coeff_imp * Y_ij * X_j for i.ne.j
      ! with Coeff_imp = C/(1+C*Y_ii)
      if (iturt(iscal).eq.20) then
        ! AFM
        coeff_imp = 1.d0+ctheta(iscal)*xtt*xiafm*gradv(ii,ii,iel)

        xut(ii,iel) = xut(ii,iel)/ coeff_imp
        temp(ii)    = temp(ii)   / coeff_imp
        ! Calculation of the diffusion tensor for the implicited part
        ! of the model computed in covofi.f90
        vistet(ii,iel) = crom(iel)*ctheta(iscal)*xtt*xrij(ii,ii)/coeff_imp

      else if(iturt(iscal).eq.21) then
        ! EB-AFM
        coeff_imp = 1.d0 + ctheta(iscal)*xtt*xi_ebafm*gradv(ii,ii,iel) &
                         + ctheta(iscal)*gamma_ebafm*xnal(ii)*xnal(ii)

        xut(ii,iel) = xut(ii,iel)/ coeff_imp
        temp(ii)    = temp(ii)   / coeff_imp
        ! Calculation of the diffusion tensor for the implicited part
        ! of the model computed in covofi.f90
        vistet(ii,iel) = crom(iel)*ctheta(iscal)*xtt*xrij(ii,ii)/coeff_imp

      else if(iturt(iscal).eq.11) then
        ! EB-GGDH
        coeff_imp = 1.d0+ ctheta(iscal)*gamma_ebggdh*xnal(ii)*xnal(ii)

        xut(ii,iel) = xut(ii,iel)/ coeff_imp
        temp(ii)    = temp(ii)   / coeff_imp
        ! Calculation of the diffusion tensor for the implicited part
        ! of the model computed in covofi.f90
        vistet(ii,iel) = crom(iel)*ctheta(iscal)*xtt*xrij(ii,ii)/coeff_imp
      endif

      ! In the next step, we compute the divergence of:
      !  "-Cp*C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2)"
      !  The part "-C_theta*k/eps* R.Grad T" is computed by the GGDH part
      w1(ii,iel) = xcpp(iel)*temp(ii)
    enddo

    ! Extra diag part of the diffusion tensor for covofi
    if(iturt(iscal).eq.11.or.iturt(iscal).eq.20.or.iturt(iscal).eq.21) then
      vistet(4,iel) = crom(iel) * ctheta(iscal)* xtt * xrij(1,2)
      vistet(5,iel) = crom(iel) * ctheta(iscal)* xtt * xrij(2,3)
      vistet(6,iel) = crom(iel) * ctheta(iscal)* xtt * xrij(1,3)
    end if

  enddo ! End loop over ncel

  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

  itypfl = 1
  iflmb0 = 1
  init   = 1
  inc    = 1
  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag

  ! Local gradient boundary conditions: homogenous Neumann
  allocate(coefat(3,ndimfb))
  allocate(coefbt(3,3,ndimfb))
  do ifac = 1, nfabor
    do ii = 1, 3
    coefat(ii,ifac) = 0.d0
      do jj = 1, 3
        if (ii.eq.jj) then
          coefbt(ii,jj,ifac) = 1.d0
        else
          coefbt(ii,jj,ifac) = 0.d0
        endif
      enddo
    enddo
  enddo

  call inimav &
  ( f_id0  , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgra , nswrgp  , imligp, &
    iwarnp ,                                              &
    epsrgp , climgp ,                                     &
    crom   , brom   ,                                     &
    w1     ,                                              &
    coefat , coefbt ,                                     &
    thflxf , thflxb )

  deallocate(coefat)
  deallocate(coefbt)
  deallocate(w1)

!===============================================================================
! 3. Transport equation on turbulent thermal fluxes (DFM)
!===============================================================================
else

  call field_get_val_prev_v(f_id, xuta)

  call resrit &
( nscal  ,                                               &
  iscal  , xcpp   , xut    , xuta   ,                    &
  dt     ,                                               &
  gradv  , gradt  , grad_al)

  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

  itypfl = 1
  iflmb0 = 1
  init   = 1
  inc    = 1
  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag

  allocate(w1(3, ncelet))

  do iel = 1, ncelet
    w1(1,iel) = xcpp(iel)*xut(1,iel)
    w1(2,iel) = xcpp(iel)*xut(2,iel)
    w1(3,iel) = xcpp(iel)*xut(3,iel)
  enddo

  ! Boundary Conditions on T'u' for the divergence term of
  ! the thermal transport equation
  call field_get_coefad_v(f_id,cofarut)
  call field_get_coefbd_v(f_id,cofbrut)

  call inimav &
  ( f_id0  , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgra , nswrgp  , imligp, &
    iwarnp ,                                              &
    epsrgp , climgp ,                                     &
    crom   , brom   ,                                     &
    w1     ,                                              &
    cofarut, cofbrut,                                     &
    thflxf , thflxb )

  deallocate(w1)

endif

!===============================================================================
! 4. Add the divergence of the thermal flux to the thermal transport equation
!===============================================================================

if (iturt(iscal).eq.11.or.ityturt(iscal).eq.2.or.ityturt(iscal).eq.3) then
  allocate(divut(ncelet))

  init = 1

  call divmas(init, thflxf, thflxb, divut)

  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) - divut(iel)
  enddo

  ! Free memory
  deallocate(divut)

endif

! Free memory
deallocate(gradv)
deallocate(gradt)
if (allocated(grad_al)) deallocate(grad_al)
deallocate(thflxf)
deallocate(thflxb)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 9999 format( &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS' ,/,&
'@'                                                            ,/,&
'@  Le calcul ne sera pas execute'                             ,/,&
'@'                                                            ,/,&
'@  Le modele de flux thermique turbulent choisi        '      ,/,&
'@  necessite le calcul de la variance du scalaire thermique'  ,/,&
'@'                                                            ,/,&
'@  Verifier les donnees entrees dans l''interface'            ,/,&
'@    et dans les sous-programmes utilisateur.'                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 9999 format( &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION'                ,/,&
'@    ========'                                                ,/,&
'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLET'  ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run                  '         ,/,&
'@'                                                            ,/,&
'@  Turbulent heat flux model taken imposed that   '           ,/,&
'@  Thermal scalar variance has to be calculate.   '           ,/,&
'@'                                                            ,/,&
'@  Verify the provided data in the interface'                 ,/,&
'@    and in user subroutines.'                                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

end subroutine

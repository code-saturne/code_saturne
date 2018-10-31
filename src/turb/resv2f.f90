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
!> \file resv2f.f90
!> \brief Resolution of source convection diffusion equations
!>        for \f$\phi\f$ and diffusion for \f$ \overline{f} \f$
!>        as part of the V2F phi-model
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     icepdc        number of ncepdp cells with head losses
!> \param[in]     icetsm        number of cells with mass source
!> \param[in]     itypsm        type of masss source for the
!>                              variables (cf. cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head losses
!> \param[in]     smacel        value of variables associated to the
!>                              mass source
!>                              for ivar=ipr,smacel=flux of masse
!> \param[in]     prdv2f        storage table of term
!>                              prod of turbulence for the v2f
!______________________________________________________________________________!


subroutine resv2f &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   prdv2f )

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
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)
double precision prdv2f(ncelet)

! Local variables

integer          init  , ifac  , iel   , inc   , iprev , iccocg
integer          ivar
integer          iiun
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp, iphydp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          istprv
integer          imucpp, idftnp, iswdyp
integer          icvflb
integer          ivoid(1)
integer          key_t_ext_id
integer          iroext
integer          iviext
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision tuexpe, thets , thetv , thetap, thetp1
double precision d2s3
double precision xk, xe, xnu, xrom, ttke, ttmin, llke, llmin, tt
double precision fhomog
double precision hint
double precision l2, time_scale

double precision rvoid(1)

character(len=80) :: label
double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, allocatable, dimension(:,:) :: gradp, gradk
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w5
double precision, allocatable, dimension(:) :: dpvar
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: cvar_fb, cvara_fb
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_phi
double precision, dimension(:), pointer :: cvar_al, cvara_al
double precision, dimension(:), pointer :: cvar_var, cvara_var
double precision, dimension(:), pointer :: cpro_pcvlo, cpro_pcvto
double precision, dimension(:), pointer :: viscl, visct
double precision, dimension(:), pointer :: c_st_phi_p, c_st_a_p

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(ncelet), rovsdt(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w5(ncelet))
allocate(dpvar(ncelet))

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_prev_s(ivarfl(ik), cvara_k)
call field_get_val_prev_s(ivarfl(iep), cvara_ep)
call field_get_val_prev_s(ivarfl(iphi), cvara_phi)
if (iturb.eq.50) then
  call field_get_val_s(ivarfl(ifb), cvar_fb)
  call field_get_val_prev_s(ivarfl(ifb), cvara_fb)
elseif (iturb.eq.51) then
  call field_get_val_s(ivarfl(ial), cvar_al)
  call field_get_val_prev_s(ivarfl(ial), cvara_al)
endif

! 2nd order previous source terms

c_st_phi_p => null()
c_st_a_p => null() ! either fb or alpha

call field_get_key_int(ivarfl(iphi), kstprv, istprv)
if (istprv.ge.0) then
  call field_get_val_s(istprv, c_st_phi_p)
  if (iturb.eq.50) then
    call field_get_key_int(ivarfl(ifb), kstprv, istprv)
    if (istprv.ge.0) then
      call field_get_val_s(istprv, c_st_a_p)
    endif
  elseif (iturb.eq.51) then
    call field_get_key_int(ivarfl(ial), kstprv, istprv)
    if (istprv.ge.0) then
      call field_get_val_s(istprv, c_st_a_p)
    endif
  endif
  if (istprv.ge.0) istprv = 1
endif

d2s3 = 2.0d0/3.0d0

call field_get_key_struct_var_cal_opt(ivarfl(iphi), vcopt)

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000)
endif

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

!===============================================================================
! 2. Calculation of term grad(phi).grad(k)
!===============================================================================

! Allocate temporary arrays gradients calculation
allocate(gradp(3,ncelet), gradk(3,ncelet))

iccocg = 1
inc = 1
iprev = 1
ivar = iphi

call field_gradient_scalar(ivarfl(ivar), iprev, imrgra, inc,      &
                           iccocg,                                &
                           gradp)

iccocg = 1
inc = 1
ivar = ik

call field_gradient_scalar(ivarfl(ik), iprev, imrgra, inc,        &
                           iccocg,                                &
                           gradk)

do iel = 1, ncel
  w1(iel) = gradp(1,iel)*gradk(1,iel) &
          + gradp(2,iel)*gradk(2,iel) &
          + gradp(3,iel)*gradk(3,iel)
enddo

! Free memory
deallocate(gradp, gradk)

!===============================================================================
! 3. Resolution of the equation of f_barre / alpha
!===============================================================================

if (iturb.eq.50) then
  ivar = ifb
  cvar_var => cvar_fb
  cvara_var => cvara_fb
elseif (iturb.eq.51) then
  ivar = ial
  cvar_var => cvar_al
  cvara_var => cvara_al
endif

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1100) label
endif

!     S as Source, V as Variable
thets  = thetst
thetv  = vcopt%thetav

call field_get_val_s(icrom, cromo)
call field_get_val_s(iviscl, cpro_pcvlo)
if (istprv.ge.0) then
  call field_get_key_int(icrom, key_t_ext_id, iroext)
  if (iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
  endif

  call field_get_key_int(iviscl, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_prev_s(iviscl, cpro_pcvlo)
  endif
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 3.1 User source terms
!===============================================================================

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(ivar)    ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   smbr   , rovsdt )

!     If we extrapolate the source terms
if (istprv.ge.0) then
  do iel = 1, ncel
    !       Save for exchange
    tuexpe = c_st_a_p(iel)
    !       For the futur and the next step time
    !       We put a mark "-" because in fact we solve
    !       \f$-\div{\grad {\dfrac{\overline{f}}{\alpha}}} = ... \f$
    c_st_a_p(iel) = - smbr(iel)
    !       Second member of the previous step time
    !       we implicit the user source term (the rest)
    smbr(iel) = - rovsdt(iel)*cvara_var(iel) - thets*tuexpe
    !       Diagonal
    rovsdt(iel) = thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
  !       We put a mark "-" because in fact we solve
  !       \f$-\div\{\grad{\dfrac{\overline{f}}{\alpha}}} = ...\f$
  !       We solve by conjugated gradient, so we do not impose the mark
  !       of rovsdt
    smbr(iel)   = -rovsdt(iel)*cvara_var(iel) - smbr(iel)
  !          rovsdt(iel) =  rovsdt(iel)
  enddo
endif


!===============================================================================
! 3.2 Source term of f_barre or alpha
!   For f_barre (phi_fbar)
!     \f[ smbr =\dfrac{1}{L^2*(f_b + \dfrac{1}{T(C1-1)(phi-2/3)}
!                                               - \dfrac{C2 Pk}{k \rho}
!     -2 \dfrac{\nu}{k\grad{\phi}\cdot \grad{k} -\nu \div{\grad{\phi} )} \f]
!   For alpha (BL-V2/K)
!     \f$smbr = \dfrac{1}{L^2 (\alpha^3 - 1)} \f$
!  In fact we put a mark "-" because the solved equation is
!    \f[ -\div{\grad{ \dfrac{\overline{f}}{\alpha}}} = smbr \f]
!===============================================================================

!     We calculate the term \f$ -VOLUME \div{\grad{\phi}} \f$ with itrgrp,
!     and we store it in W2
!     Warning, the viscf and viscb calculated here are of use for itrgr but
!     for codits too

do iel = 1, ncel
  w3(iel) = 1.d0
enddo
call viscfa                                                       &
 ( imvisf ,                                                       &
   w3     ,                                                       &
   viscf  , viscb  )

! Translate coefa into cofaf and coefb into cofbf

call field_get_coefa_s(ivarfl(iphi), coefap)
call field_get_coefb_s(ivarfl(iphi), coefbp)
call field_get_coefaf_s(ivarfl(iphi), cofafp)
call field_get_coefbf_s(ivarfl(iphi), cofbfp)

do ifac = 1, nfabor

  iel = ifabor(ifac)

  hint = w3(iel)/distb(ifac)

  ! Translate coefa into cofaf and coefb into cofbf
  cofafp(ifac) = -hint*coefap(ifac)
  cofbfp(ifac) = hint*(1.d0-coefbp(ifac))

enddo

iccocg = 1
inc = 1
init = 1

call field_get_key_struct_var_cal_opt(ivarfl(iphi), vcopt)

nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
iwarnp = vcopt%iwarni
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
iphydp = 0

call itrgrp &
 ( ivarfl(iphi)    , init   , inc    , imrgra ,                      &
   iccocg , nswrgp , imligp , iphydp ,                               &
   iwarnp ,                                                          &
   epsrgp , climgp , extrap ,                                        &
   rvoid  ,                                                          &
   cvara_phi       ,                                                 &
   coefap , coefbp , cofafp , cofbfp ,                               &
   viscf  , viscb  ,                                                 &
   w3     ,                                                          &
   w2     )

!     Explicit term, stores ke temporarily in w5
!     w2 is already multipicated by the volume which already contains
!     a mark "-" (coming from itrgrp)
do iel = 1, ncel
  ! Compute the time scale
  xk = cvara_k(iel)
  xe = cvara_ep(iel)
  xrom = cromo(iel)
  xnu  = cpro_pcvlo(iel)/xrom
  ttke = xk / xe
  if (iturb.eq.50) then
    ttmin = cv2fct*sqrt(xnu/xe)
    time_scale = max(ttke,ttmin)
  elseif (iturb.eq.51) then
    ttmin = cpalct*sqrt(xnu/xe)
    time_scale = sqrt(ttke**2 + ttmin**2)
  endif

  if (iturb.eq.50) then
    w5(iel) = - cell_f_vol(iel)*                                    &
         ( (cv2fc1-1.d0)*(cvara_phi(iel)-d2s3)/time_scale           &
               -cv2fc2*prdv2f(iel)/xrom/xk                          &
          -2.0d0*xnu/xe/time_scale * w1(iel) ) - xnu*w2(iel)
  elseif (iturb.eq.51) then
    w5(iel) = cell_f_vol(iel)
  endif
enddo

! If we extrapolate the source term:
if (istprv.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    c_st_a_p(iel) = c_st_a_p(iel) + w5(iel)
    smbr(iel) = smbr(iel) + thetp1*c_st_a_p(iel)
  enddo
! Otherwise: smbr
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w5(iel)
  enddo
endif

if (istprv.ge.0) then
  thetap = thetv
else
  thetap = 1.d0
endif

! Implicit term
do iel = 1, ncel

  ! Compute L^2
  xk = cvara_k(iel)
  xe = cvara_ep(iel)
  xnu  = viscl(iel)/crom(iel)

  llke = xk**1.5d0/xe
  if (iturb.eq.50) then
    llmin = cv2fet*(xnu**3/xe)**0.25d0
    l2 = ( cv2fcl*max(llke,llmin) )**2
  elseif (iturb.eq.51) then
    llmin = cpalet*(xnu**3/xe)**0.25d0
    l2 = cpalcl**2*(llke**2 + llmin**2)
  endif

  smbr(iel) = ( - cell_f_vol(iel)*cvara_var(iel) + smbr(iel) ) / l2
  ! ---> Matrix
  rovsdt(iel) = (rovsdt(iel) + cell_f_vol(iel)*thetap) / l2
enddo

!===============================================================================
! 3.3 Effective resolution in the equation of f_barre or alpha
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

iconvp = vcopt%iconv
idiffp = vcopt%idiff
ndircp = vcopt%ndircl
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
! all boundary convective flux with upwind
icvflb = 0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits &
 ( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , cvar_var        , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 4. Resolution of the equation of phi
!===============================================================================

ivar = iphi

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

call field_get_val_s(ivarfl(ivar), cvar_var)
cvara_var => cvara_phi

if (vcopt%iwarni.ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1100) label
endif

!     S as Source, V as Variable
thets  = thetst
thetv  = vcopt%thetav

call field_get_val_s(ivisct, cpro_pcvto)
if (istprv.ge.0) then
  call field_get_key_int(ivisct, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_prev_s(ivisct, cpro_pcvto)
  endif
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 4.1 User source terms
!===============================================================================

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(ivar)    ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   smbr   , rovsdt )

!     If we extrapolate the source terms
if (istprv.ge.0) then
  do iel = 1, ncel
    !       Save for exchange
    tuexpe = c_st_phi_p(iel)
    !       For the future and the next time step
    c_st_phi_p(iel) = smbr(iel)
    !       Second member of previous time step
    !       We suppose -rovsdt > 0: we implicit
    !       the user source term (the rest)
    smbr(iel) = rovsdt(iel)*cvara_var(iel) - thets*tuexpe
    !       Diagonal
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel)   = rovsdt(iel)*cvara_var(iel) + smbr(iel)
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

!===============================================================================
! 4.2 Mass source term
!===============================================================================


if (ncesmp.gt.0) then

  !       Integer equal to 1 (for navsto: nb of over-iter)
  iiun = 1

  !       We increment smbr by -Gamma.var_prev and rovsdt by Gamma
  call catsma                                                     &
 ( ncelet , ncel   , ncesmp , iiun   , isto2t ,                   &
   icetsm , itypsm(1,ivar) ,                                      &
   cell_f_vol , cvara_var  , smacel(1,ivar) , smacel(1,ipr) ,     &
   smbr   ,  rovsdt , w2 )

  ! If we extrapolate the source term we put Gamma Pinj in the prev. TS
  if (istprv.ge.0) then
    do iel = 1, ncel
      c_st_phi_p(iel) = c_st_phi_p(iel) + w2(iel)
    enddo
  !       Otherwise we put it directly in smbr
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w2(iel)
    enddo
  endif

endif


!===============================================================================
! 4.3 Mass accumulation term \f$ -\dfrad{dRO}{dt}VOLUME \f$
!    and unstable over time term
!===============================================================================

! ---> Adding the matrix diagonal

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
           + vcopt%istat*(crom(iel)/dt(iel))*cell_f_vol(iel)
enddo

!===============================================================================
! 4.4 Source term of phi
!     \f$ \phi_fbar\f$:
!     \f[ smbr = \rho f_barre - \dfrac{\phi}{k} P_k +\dfrac{2}{k}
!                         \dfrac{\mu_t}{\sigma_k} \grad{\phi} \cdot \grad{k} \f]
!     BL-V2/K:
!     \f[ smbr = \rho \alpha f_h + \rho (1-\alpha^p) f_w - \dfrac{\phi}{k} P_k
!          +\dfrac{2}{k} \dfrac{\mu_t}{\sigma_k} \grad{\phi} \cdot \grad{k} \f]
!        with \f$ f_w=-\dfrac{\epsilon}{2} \cdot \dfrac{\phi}{k} \f$ and
!             \f$ f_h = \dfrac{1}{T} \cdot
!                                (C1-1+C2 \dfrac{P_k}{\epsilon \rho} (2/3-\phi)
!===============================================================================

!     Exmplicit term, store temporarily in W2

do iel = 1, ncel
  xk = cvara_k(iel)
  xe = cvara_ep(iel)
  xnu  = viscl(iel)/crom(iel)
  if (iturb.eq.50) then
    ! The term in f_bar is taken at the current and not previous time step
    ! ... a priori better
    ! Remark: if we keep this choice, we have to modify the case
    !         of the second-order (which need the previous value time step
    !         for extrapolation).
    w2(iel)   =  cell_f_vol(iel)*                                   &
         ( xrom*cvar_fb(iel)                                        &
           +2.d0/xk*cpro_pcvto(iel)/sigmak*w1(iel) )
  elseif (iturb.eq.51) then
    ttke = xk / xe
    ttmin = cpalct*sqrt(xnu/xe)
    tt = sqrt(ttke**2 + ttmin**2)
    fhomog = -1.d0/tt*(cpalc1-1.d0+cpalc2*prdv2f(iel)/xe/xrom)*     &
             (cvara_phi(iel)-d2s3)
    w2(iel)   = cell_f_vol(iel)*                                    &
         ( cvara_al(iel)**3*fhomog*xrom                             &
           +2.d0/xk*cpro_pcvto(iel)/sigmak*w1(iel) )
         !FIXME implit negativ w1 and fhomog
  endif

enddo

! If we extrapolate the source term: prev. TS
if (istprv.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    c_st_phi_p(iel) = c_st_phi_p(iel) + w2(iel)
    smbr(iel) = smbr(iel) + thetp1*c_st_phi_p(iel)
  enddo
!     Otherwise: smbr
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w2(iel)
  enddo
endif

!     Implicit term
do iel = 1, ncel
  xrom = cromo(iel)
  if (iturb.eq.50) then
    smbr(iel) = smbr(iel)                                         &
         - cell_f_vol(iel)*prdv2f(iel)*cvara_phi(iel)/cvara_k(iel)
  elseif (iturb.eq.51) then
    smbr(iel) = smbr(iel)                                         &
              - cell_f_vol(iel)*(prdv2f(iel)+xrom*cvara_ep(iel)/2 &
                                    *(1.d0-cvara_al(iel)**3))     &
         *cvara_phi(iel)/cvara_k(iel)
  endif
enddo

! ---> Matrix

if (istprv.ge.0) then
  thetap = thetv
else
  thetap = 1.d0
endif
do iel = 1, ncel
  xrom = cromo(iel)
  if (iturb.eq.50) then
    rovsdt(iel) = rovsdt(iel)                                     &
         + cell_f_vol(iel)*max(prdv2f(iel),0.d0)/cvara_k(iel)*thetap
  elseif (iturb.eq.51) then
    rovsdt(iel) = rovsdt(iel)                                     &
         + cell_f_vol(iel)*(max(prdv2f(iel),0.d0)+xrom*cvara_ep(iel)/2    &
                                    *(1.d0-cvara_al(iel)**3))         &
           /cvara_k(iel)*thetap
  endif
enddo

!===============================================================================
! 4.5 Diffusion terms
!===============================================================================
! ---> Viscosity
! Normally, in the phi-model equations, only turbulent viscosity
!  turbulente takes place in phi diffusion (the term with mu disappeared
!  passing from \f$f\f$ to \f$ \overline{f})\f$. But as it stands,
!  it makes the calculation unstable (because \f$\mu_t\f$ tends towards 0
!  at the wall what decouples \f$ \phi \f$ of its boundary condition and
!  the molecular diffusion term is integred in \f$ \overline{f} \f$, it is as if it
!  was treated as explicit)
!  -> we add artificially diffusion (knowing that as k=0, the phi value
!  does not matter)

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

if (vcopt%idiff.ge.1) then
  do iel = 1, ncel
    if (iturb.eq.50) then
      w2(iel) = viscl(iel)      + visct(iel)/sigmak
    elseif (iturb.eq.51) then
      w2(iel) = viscl(iel)/2.d0 + visct(iel)/sigmak !FIXME
    endif
  enddo

  call viscfa &
( imvisf ,                                                       &
  w2     ,                                                       &
  viscf  , viscb  )

  ! Translate coefa into cofaf and coefb into cofbf
  do ifac = 1, nfabor

    iel = ifabor(ifac)

    hint = w2(iel)/distb(ifac)

    ! Translate coefa into cofaf and coefb into cofbf
    cofafp(ifac) = -hint*coefap(ifac)
    cofbfp(ifac) = hint*(1.d0-coefbp(ifac))

  enddo

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0

    ! Translate coefa into cofaf and coefb into cofbf
    cofafp(ifac) = 0.d0
    cofbfp(ifac) = 0.d0
  enddo

endif

!===============================================================================
! 4.6 Effective resolution of the phi equation
!===============================================================================

if (istprv.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*c_st_phi_p(iel)
  enddo
endif

iconvp = vcopt%iconv
idiffp = vcopt%idiff
ndircp = vcopt%ndircl
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
! all boundary convective flux with upwind
icvflb = 0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits &
 ( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , cvar_var        , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 5. Clipping
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(iphi), vcopt)

call clpv2f(ncel, vcopt%iwarni)

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, rovsdt)
deallocate(w1, w2, w3)
deallocate(w5)
deallocate(dpvar)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                         &
'   ** Resolution pour V2F (phi et f_bar/alpha)               ',/,&
'      ----------------------------------------        ',/)
 1100    format(/,'           Resolution pour la variable ',A8,/)

#else

 1000    format(/,                                         &
'   ** Solving V2F (phi and f_bar/alpha)'               ,/,&
'      ---------------------------------'               ,/)
 1100    format(/,'           Solving variable ',A8                  ,/)

#endif

!----
! End
!----

return

end subroutine

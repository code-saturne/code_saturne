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

!> \file turbkw.f90
!>
!> \brief Solving the \f$ k - \omega \f$ SST for incompressible flows
!> or slightly compressible flows for one time step.
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

subroutine turbkw &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   tslagr , ckupdc , smacel )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstnum
use cstphy
use optcal
use lagran
use parall
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
double precision tslagr(ncelet,*)
double precision ckupdc(6,ncepdp), smacel(ncesmp,nvar)

! Local variables

character(len=80) :: chaine
integer          iel   , ifac  , inc   , iprev,  iccocg, ivar
integer          ii, f_id , iiun
integer          iclipk(1), iclipw, iclpkmx(1)
integer          nswrgp, imligp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          iflmas, iflmab
integer          iwarnp
integer          istprv
integer          init
integer          imucpp, idftnp, iswdyp
integer          key_t_ext_id
integer          iroext
integer          iviext
integer          kclipp, clip_w_id, clip_k_id

integer          icvflb, imasac
integer          ivoid(1)

double precision rnorm , d2s3, divp23, epz2
double precision deltk , deltw, a11, a12, a22, a21
double precision unsdet, romvsd
double precision prdtur, xk, xw, xeps, xnu
double precision visct , rho, visclc, visctc, hint
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision thetp1, thetak, thetaw, thets, thetap, epsrsp
double precision tuexpk, tuexpw
double precision cdkw, xarg1, xxf1, xgamma, xbeta, sigma, produc
double precision xlt, xdelta, xrd, xfd, xs2pw2, xdist, xdiff, fddes
double precision var, vrmin(2), vrmax(2)
double precision utaurf,ut2,ypa,ya,xunorm, limiter, nu0
double precision turb_schmidt
double precision normp

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrk, smbrw
double precision, allocatable, dimension(:) :: tinstk, tinstw, xf1
double precision, allocatable, dimension(:,:) :: gradk, grado, grad
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:) :: s2pw2
double precision, allocatable, dimension(:) :: w1, w2
double precision, allocatable, dimension(:) :: gdkgdw
double precision, allocatable, dimension(:) :: w5, w6
double precision, allocatable, dimension(:) :: prodk, prodw
double precision, allocatable, dimension(:) :: gamk, gamw
double precision, allocatable, dimension(:) :: usimpk, usimpw
double precision, allocatable, dimension(:) :: w7
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: rotfct
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: bromo, cromo
double precision, dimension(:), pointer :: coefa_k, coefb_k, coefaf_k, coefbf_k
double precision, dimension(:), pointer :: coefa_o, coefb_o, coefaf_o, coefbf_o
double precision, dimension(:), pointer :: cvar_k, cvara_k, cvar_omg, cvara_omg
double precision, dimension(:), pointer :: cvar_var
double precision, dimension(:), pointer :: cpro_pcvto, cpro_pcvlo
double precision, dimension(:), pointer :: viscl, cvisct
double precision, dimension(:), pointer :: c_st_k_p, c_st_omg_p
double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cpro_divukw, cpro_s2kw
double precision, dimension(:), pointer :: ddes_fd_coeff
double precision, dimension(:), pointer :: w_dist
double precision, dimension(:), pointer :: cpro_k_clipped
double precision, dimension(:), pointer :: cpro_w_clipped

type(var_cal_opt) :: vcopt_w, vcopt_k

!===============================================================================

!===============================================================================
! 1.Initialization
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrk(ncelet), smbrw(ncelet))
allocate(tinstk(ncelet), tinstw(ncelet), xf1(ncelet))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet))
allocate(dpvar(ncelet))
allocate(gdkgdw(ncelet))
allocate(prodk(ncelet), prodw(ncelet))
if (iddes.eq.1) then
  allocate(gradv(3, 3, ncelet))
  allocate(s2pw2(ncelet))
endif

epz2 = epzero**2
iclpkmx(1) = 0

call field_get_val_s(ivisct, cvisct)
call field_get_val_s(iviscl, viscl)

call field_get_key_int(ivarfl(ik), kimasf, iflmas)
call field_get_key_int(ivarfl(ik), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, cromo)
call field_get_val_s(ibrom, bromo)

call field_get_val_s(ivisct, cpro_pcvto)
call field_get_val_s(iviscl, cpro_pcvlo)

call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_prev_s(ivarfl(ik), cvara_k)
call field_get_val_s(ivarfl(iomg), cvar_omg)
call field_get_val_prev_s(ivarfl(iomg), cvara_omg)

call field_get_coefa_s(ivarfl(ik), coefa_k)
call field_get_coefb_s(ivarfl(ik), coefb_k)
call field_get_coefaf_s(ivarfl(ik), coefaf_k)
call field_get_coefbf_s(ivarfl(ik), coefbf_k)

call field_get_coefa_s(ivarfl(iomg), coefa_o)
call field_get_coefb_s(ivarfl(iomg), coefb_o)
call field_get_coefaf_s(ivarfl(iomg), coefaf_o)
call field_get_coefbf_s(ivarfl(iomg), coefbf_o)

call field_get_key_struct_var_cal_opt(ivarfl(iomg), vcopt_w)
call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt_k)

thets  = thetst

call field_get_val_s(is2kw, cpro_s2kw)
call field_get_val_s(idivukw, cpro_divukw)

call field_get_key_int(ivarfl(ik), kstprv, istprv)
if (istprv.ge.0) then
  call field_get_val_s(istprv, c_st_k_p)
  call field_get_key_int(ivarfl(iep), kstprv, istprv)
  if (istprv.ge.0) then
    call field_get_val_s(istprv, c_st_omg_p)
  endif
  if (istprv.ge.0) istprv = 1
endif

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

if (istprv.ge.0) then
  call field_get_key_int(icrom, key_t_ext_id, iroext)
  if (iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
    call field_get_val_prev_s(ibrom, bromo)
  endif
  call field_get_key_int(iviscl, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_prev_s(iviscl, cpro_pcvlo)
  endif
  call field_get_key_int(ivisct, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_prev_s(ivisct, cpro_pcvto)
  endif
endif

call field_get_id("wall_distance", f_id)
call field_get_val_s(f_id, w_dist)

if (vcopt_k%iwarni.ge.1) then
  write(nfecra,1000)
endif

d2s3 = 2.d0/3.d0

!===============================================================================
! 2.1 Compute dk/dxj.dw/dxj
!     stored in gdkgdw
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradk(3,ncelet), grado(3,ncelet))

iccocg = 1
inc = 1
iprev = 1

call field_gradient_scalar(ivarfl(ik), iprev, imrgra, inc,          &
                           iccocg,                                  &
                           gradk)

call field_gradient_scalar(ivarfl(iomg), iprev, imrgra, inc,        &
                           iccocg,                                  &
                           grado)

if (iddes.eq.1) then
  call field_get_val_s_by_name("hybrid_blend", ddes_fd_coeff)

  call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,        &
                             gradv)
  do iel = 1, ncel
    s2pw2(iel) = gradv(1,1,iel)**2.d0 + gradv(2,1,iel)**2.d0 + gradv(3,1,iel)**2.d0 &
               + gradv(1,2,iel)**2.d0 + gradv(2,2,iel)**2.d0 + gradv(3,2,iel)**2.d0 &
               + gradv(1,3,iel)**2.d0 + gradv(2,3,iel)**2.d0 + gradv(3,3,iel)**2.d0
  enddo
end if

do iel = 1, ncel
  gdkgdw(iel) = gradk(1,iel)*grado(1,iel) &
              + gradk(2,iel)*grado(2,iel) &
              + gradk(3,iel)*grado(3,iel)
enddo

! Free memory
deallocate(gradk, grado)

!===============================================================================
! 2.2. Compute the weight f1 (stored in xf1)
!===============================================================================

do iel = 1, ncel
  w2(iel) = max(w_dist(iel),epzero)
enddo

! En cas d'ordre 2 on utilise les valeurs en n car le terme en (1-f1)*gdkgdw
! sera une propriété. Du coup, on aura quand meme certaines "constantes"
! intervenant dans des termes en n+1/2 (ex sigma_k pour la diffusion) calcules
! a partir de f1 en n -> mais l'effet sur les "constantes" est faible
! -> a garder en tete si on fait vraiment de l'ordre 2 en temps en k-omega
do iel = 1, ncel
  rho = cromo(iel)
  xnu = cpro_pcvlo(iel)/rho
  xk = cvara_k(iel)
  xw  = cvara_omg(iel)
  cdkw = 2*rho/ckwsw2/xw*gdkgdw(iel)
  cdkw = max(cdkw,1.d-20)
  xarg1 = max(sqrt(xk)/cmu/xw/w2(iel), 500.d0*xnu/xw/w2(iel)**2)
  xarg1 = min(xarg1, 4.d0*rho*xk/ckwsw2/cdkw/w2(iel)**2)
  xf1(iel) = tanh(xarg1**4)
enddo

!===============================================================================
! 3. Unsteady terms (stored in tinstk and tinstw)
!===============================================================================

do iel = 1, ncel
  rho = crom(iel)
  romvsd = rho*volume(iel)/dt(iel)
  tinstk(iel) = vcopt_k%istat*romvsd
  tinstw(iel) = vcopt_w%istat*romvsd
enddo

!===============================================================================
! 4. Compute production terms
!     stored in: prodk,prodw
!      En sortie de l'etape on conserve gdkgdw,xf1,prodk,tinstW
!===============================================================================

do iel = 1, ncel
  xk   = cvara_k(iel)
  xw   = cvara_omg(iel)
  xeps = cmu*xw*xk
  visct = cpro_pcvto(iel)
  rho = cromo(iel)
  prodw(iel) = visct*cpro_s2kw(iel)                    &
             - d2s3*rho*xk*cpro_divukw(iel)

  ! The negative part is implicited
  xxf1   = xf1(iel)
  xgamma = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
  tinstw(iel) = tinstw(iel)                                         &
              + max( d2s3*rho*volume(iel)                           &
                   *(rho*xgamma*xk/(visct*xw))*cpro_divukw(iel), 0.d0)

  ! Take the min between prodw and the low Reynold one
  if (prodw(iel).gt.ckwc1*rho*xeps) then
    prodk(iel) = ckwc1*rho*xeps
  else
    prodk(iel) = prodw(iel)
    tinstk(iel) = tinstk(iel) + max(d2s3*volume(iel)*rho*cpro_divukw(iel), 0.d0)
  endif
enddo

!===============================================================================
! 5. Take into account rotation/curvature correction, if necessary
!===============================================================================

! Spalart-Shur correction: the production terms are multiplied by a
! 'rotation function'

if (irccor.eq.1) then

  ! Allocate an array for the rotation function
  allocate(rotfct(ncel))

  ! Compute the rotation function (gdkgdw array not used)
  call rotcor(dt, rotfct, gdkgdw)

  do iel = 1, ncel
    prodk(iel) = prodk(iel)*rotfct(iel)
    prodw(iel) = prodw(iel)*rotfct(iel)
  enddo

  ! rotfct array is used later in case of renforced coupling (ikecou = 1).
  ! The array is deallocated at the end of the subroutine.

endif

!===============================================================================
! 6. Compute buoyancy terms
!     stored in: prodk, prodw, w2
!===============================================================================

if (igrake.eq.1) then

  ! Allocate a temporary array for the gradient calculation
  allocate(grad(3,ncelet))

  ! --- Buoyant term:     G = Beta*g*GRAD(T)/PrT/rho
  !     Here is computed: G =-g*GRAD(rho)/PrT/rho

  iccocg = 1
  inc = 1

  nswrgp = vcopt_k%nswrgr
  epsrgp = vcopt_k%epsrgr
  imligp = vcopt_k%imligr
  iwarnp = vcopt_k%iwarni
  climgp = vcopt_k%climgr
  extrap = vcopt_k%extrag

  ! BCs on rho: Dirichlet ROMB
  !  NB: viscb is used as COEFB

  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  f_id = -1

  call gradient_s                                                 &
 ( f_id   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , epsrgp , climgp , extrap ,                            &
   cromo  , bromo  , viscb  ,                                     &
   grad   )


  ! Buoyancy production
  !   prodk=min(P,c1*eps)+G
  !   prodw=P+(1-ce3)*G
  if (iscalt.gt.0.and.nscal.ge.iscalt) then
    call field_get_key_double(ivarfl(isca(iscalt)), ksigmas, turb_schmidt)
    prdtur = turb_schmidt
  else
    prdtur = 1.d0
  endif

  do iel = 1, ncel
    rho = cromo(iel)
    visct = cpro_pcvto(iel)

    w2(iel) = -(grad(1,iel)*gx + grad(2,iel)*gy + grad(3,iel)*gz) / &
               (rho*prdtur)

    prodw(iel) = prodw(iel)+visct*max(w2(iel),zero)
    prodk(iel) = prodk(iel)+visct*w2(iel)

    ! Implicit Buoyant terms when negativ
    tinstk(iel) = tinstk(iel)                                        &
                + max(-volume(iel)*visct/cvara_k(iel)*w2(iel), 0.d0)
  enddo

  ! Free memory
  deallocate(grad)

endif

!===============================================================================
! 7. Take user source terms into account
!     explicit parts stored in: smbrk, smbrw
!     implicit parts stored in: usimpk, usimpw
!===============================================================================

allocate(usimpk(ncelet), usimpw(ncelet))

do iel = 1, ncel
  smbrk(iel) = 0.d0
  smbrw(iel) = 0.d0
  usimpk(iel) = 0.d0
  usimpw(iel) = 0.d0
enddo

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(ik)      ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   smbrk  , usimpk )

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(iomg)    ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   smbrw  , usimpw )

! If source terms are extrapolated over time
if (istprv.ge.0) then

  thetak = vcopt_k%thetav
  thetaw = vcopt_w%thetav

  do iel = 1, ncel

    ! Recover the value at time (n-1)
    tuexpk = c_st_k_p(iel)
    tuexpw = c_st_omg_p(iel)

    ! Save the values for the next time-step
    c_st_k_p(iel) = smbrk(iel)
    c_st_omg_p(iel) = smbrw(iel)

    ! Explicit Part
    smbrk(iel) = - thets*tuexpk
    smbrw(iel) = - thets*tuexpw
    ! It is assumed that (-usimpk > 0) and though this term is implicit
    smbrk(iel) = usimpk(iel)*cvara_k(iel) + smbrk(iel)
    smbrw(iel) = usimpw(iel)*cvara_omg(iel) + smbrw(iel)

    ! Implicit part
    tinstk(iel) = tinstk(iel) -usimpk(iel)*thetak
    tinstw(iel) = tinstw(iel) -usimpw(iel)*thetaw
  enddo

! If no extrapolation over time
else
  do iel = 1, ncel
    ! Explicit Part
    smbrk(iel) = smbrk(iel) + usimpk(iel)*cvara_k(iel)
    smbrw(iel) = smbrw(iel) + usimpw(iel)*cvara_omg(iel)

    ! Implicit part
    tinstk(iel) = tinstk(iel) + max(-usimpk(iel),zero)
    tinstw(iel) = tinstw(iel) + max(-usimpw(iel),zero)
  enddo
endif


!===============================================================================
! 8. Finalization of explicit and implicit source terms

!      Les termes sont stockes dans     smbrk, smbrw
!      En sortie de l'etape on conserve smbrk,smbrw,gdkgdw
!===============================================================================
! Standard k-w SST RANS model
if(iddes.ne.1) then
  do iel = 1, ncel

    visct  = cpro_pcvto(iel)
    rho    = cromo(iel)
    xk     = cvara_k(iel)
    xw     = cvara_omg(iel)
    xxf1   = xf1(iel)
    xgamma = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
    xbeta  = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2

    smbrk(iel) = smbrk(iel) + volume(iel)*(                         &
                                            prodk(iel)              &
                                          - cmu*rho*xw*xk )

    smbrw(iel) = smbrw(iel)                                                &
               + volume(iel)*(                                             &
                               rho*xgamma/visct*prodw(iel)                 &
                             - xbeta*rho*xw**2                             &
                             + 2.d0*rho/xw*(1.d0-xxf1)/ckwsw2*gdkgdw(iel)  &
                             )

    tinstw(iel) = tinstw(iel) + volume(iel)*max(-2.d0*rho/xw**2*(1.d0-xxf1) &
                                                /ckwsw2*gdkgdw(iel), 0.d0)
  enddo
! DDES mode for k-w SST
else
  do iel = 1, ncel

    visct  = cpro_pcvto(iel)
    xnu    = cpro_pcvlo(iel)
    rho    = cromo(iel)
    xk     = cvara_k(iel)
    xw     = cvara_omg(iel)
    xxf1   = xf1(iel)
    xgamma = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
    xbeta  = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2
    xs2pw2 = max(sqrt(s2pw2(iel)),epzero)
    xdist  = max(w_dist(iel),epzero)

    xlt    = sqrt(xk)/(cmu*xw)
    xdelta = volume(iel)**(1.d0/3.d0)
    xrd    = (visct + xnu)/ ( xs2pw2 *(xkappa**2) *(xdist**2))
    xfd    = 1.d0 -tanh((8.d0*xrd)**3)
    xdiff  = max((xlt-(cddes*xdelta)),0.d0)
    fddes  = xlt/(xlt-xfd*xdiff)

    ddes_fd_coeff(iel) = xfd

    ! Storage of the Fd coefficient and check locally if
    ! DDES is activated for post-propcessing
    ! NB: for RANS zones (L_RANS < L_LES) Fd is clipped to 0
    if ((cddes*xdelta).ge.xlt)then
      ddes_fd_coeff(iel) = 0.d0
    endif

    smbrk(iel) = smbrk(iel) + volume(iel)*(                         &
                                            prodk(iel)              &
                                          - cmu*rho*xw*xk*fddes)

    smbrw(iel) = smbrw(iel)                                                &
               + volume(iel)*(                                             &
                               rho*xgamma/visct*prodw(iel)                 &
                             - xbeta*rho*xw**2                             &
                             + 2.d0*rho/xw*(1.d0-xxf1)/ckwsw2*gdkgdw(iel)  &
                             )

    tinstw(iel) = tinstw(iel) + volume(iel)*max(-2.d0*rho/xw**2*(1.d0-xxf1) &
                                                /ckwsw2*gdkgdw(iel), 0.d0)
  enddo
endif

! If the solving of k-omega is uncoupled, negative source terms are implicited
if (ikecou.eq.0) then
  do iel=1,ncel
    xw    = cvara_omg(iel)
    xxf1  = xf1(iel)
    xbeta = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2
    rho = crom(iel)
    tinstk(iel) = tinstk(iel) + volume(iel)*cmu*rho*xw
    tinstw(iel) = tinstw(iel) + 2.d0*volume(iel)*xbeta*rho*xw
  enddo
endif

! Free memory
deallocate(gdkgdw)

!===============================================================================
! 9 Prise en compte des termes sources lagrangien
!   couplage retour
!===============================================================================

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1, ncel

    ! Termes sources explicte et implicte sur k
    smbrk(iel)  = smbrk(iel) + tslagr(iel,itske)

    ! Termes sources explicte sur omega : on reprend la constante CE4 directement
    !    du k-eps sans justification ... a creuser si necessaire
    smbrw(iel)  = smbrw(iel)                                      &
                + ce4 *tslagr(iel,itske) * cromo(iel)             &
                /cpro_pcvto(iel)

    ! Termes sources implicite sur k
    tinstk(iel) = tinstk(iel) + max(-tslagr(iel,itsli),zero)

    ! Termes sources implicte sur omega
    tinstw(iel) = tinstw(iel)                                     &
                + max( (-ce4*tslagr(iel,itske)/cvara_k(iel)) , zero)
  enddo

endif

!===============================================================================
! 10. Mass source terms (Implicit and explicit parts)

!===============================================================================

if (ncesmp.gt.0) then

  allocate(gamk(ncelet), gamw(ncelet))

  ! Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

  ! On incremente SMBRS par -Gamma.var_prev et ROVSDT par Gamma
  ivar = ik

  call catsma &
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t ,                                                       &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , cvara_k      , smacel(1,ivar) , smacel(1,ipr) ,       &
   smbrk  , tinstk , gamk )

  ivar = iomg

  call catsma &
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t ,                                                       &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , cvara_omg    , smacel(1,ivar) , smacel(1,ipr) ,       &
   smbrw  , tinstw , gamw )

  ! Si on extrapole les TS on met Gamma Pinj dans c_st_k_p, c_st_omg_p
  if (istprv.ge.0) then
    do iel = 1, ncel
      c_st_k_p(iel) = c_st_k_p(iel) + gamk(iel)
      c_st_omg_p(iel) = c_st_omg_p(iel) + gamw(iel)
    enddo
  !  Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbrk(iel) = smbrk(iel) + gamk(iel)
      smbrw(iel) = smbrw(iel) + gamw(iel)
    enddo
  endif

  !Free memory
  deallocate(gamk, gamw)

endif

! Finalisation des termes sources
if (istprv.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) + thetp1 * c_st_k_p(iel)
    smbrw(iel) = smbrw(iel) + thetp1 * c_st_omg_p(iel)
  enddo
endif


!===============================================================================
! 11.1 Re-set Boundary conditions flux coefficient for k and omega

!     The definition of cofaf requires hint=(mu+muT/sigma)/distb where sigma
!     is not constant in the k-omega model (and not directly accessible)
!===============================================================================

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! --- Physical Propreties
  visclc = viscl(iel)
  visctc = cvisct(iel)

  xxf1 = xf1(iel)

  ! k
  sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
  hint = (visclc+visctc/sigma)/distb(ifac)

  ! Translate coefa into cofaf and coefb into cofbf
  coefaf_k(ifac) = -hint*coefa_k(ifac)
  coefbf_k(ifac) = hint*(1.d0-coefb_k(ifac))

  ! Omega
  sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
  hint = (visclc+visctc/sigma)/distb(ifac)

  ! Translate coefa into cofaf and coefb into cofbf
  coefaf_o(ifac) = -hint*coefa_o(ifac)
  coefbf_o(ifac) = hint*(1.d0-coefb_o(ifac))

enddo

!===============================================================================
! 11.2 Prise en compte des termes de conv/diff dans le second membre pour le
!      couplage renforcé k-omega (ikecou == 1)

!      Tableaux de travail              w7
!      Les termes sont stockes dans     w5 ET w6, PUIS AJOUTES A smbrk, smbrw
!      En sortie de l'etape on conserve w2-6,smbrk,smbrw,usimpk
!===============================================================================

if (ikecou.eq.1) then

  allocate(w5(ncelet), w6(ncelet))
  allocate(w7(ncelet))

  do iel = 1, ncel
    w5 (iel) = 0.d0
    w6 (iel) = 0.d0
  enddo

  ! ---> Traitement de k
  ivar   = ik

  call field_get_label(ivarfl(ivar), chaine)

  if (vcopt_k%idiff.ge. 1) then

    do iel = 1, ncel
      xxf1 = xf1(iel)
      sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
      w7(iel) = viscl(iel)                                        &
              + vcopt_k%idifft*cvisct(iel)/sigma
    enddo
    call viscfa &
 ( imvisf ,                                                       &
   w7     ,                                                       &
   viscf  , viscb  )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

  iccocg = 1
  inc    = 1
  imucpp = 0
  idftnp = ISOTROPIC_DIFFUSION
  imasac = 1
  iconvp = vcopt_k%iconv
  idiffp = vcopt_k%idiff
  nswrgp = vcopt_k%nswrgr
  imligp = vcopt_k%imligr
  ircflp = vcopt_k%ircflu
  ischcp = vcopt_k%ischcv
  isstpp = vcopt_k%isstpc
  iwarnp = vcopt_k%iwarni
  blencp = vcopt_k%blencv
  epsrgp = vcopt_k%epsrgr
  climgp = vcopt_k%climgr
  extrap = vcopt_k%extrag
  relaxp = vcopt_k%relaxv
  thetap = vcopt_k%thetav
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
 ( idtvar , ivarfl(ivar)    , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                            &
   iwarnp , imucpp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,                   &
   cvara_k         , cvara_k         ,                                     &
   coefa_k , coefb_k , coefaf_k , coefbf_k ,                               &
   imasfl , bmasfl ,                                                       &
   viscf  , viscb  , rvoid  , rvoid  ,                                     &
   rvoid  , rvoid  ,                                                       &
   icvflb , ivoid  ,                                                       &
   w5     )

  if (vcopt_k%iwarni.ge.2) then
    rnorm = sqrt(cs_gdot(ncel,smbrk,smbrk))
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif

  ! ---> Traitement de omega
  ivar   = iomg

  call field_get_label(ivarfl(ivar), chaine)

  if (vcopt_w%idiff.ge. 1) then
    do iel = 1, ncel
      xxf1 = xf1(iel)
      sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
      w7(iel) = viscl(iel)                                        &
              + vcopt_w%idifft*cvisct(iel)/sigma
    enddo
    call viscfa &
 ( imvisf ,                                                       &
   w7     ,                                                       &
   viscf  , viscb  )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

  iccocg = 1
  inc    = 1
  imucpp = 0
  idftnp = ISOTROPIC_DIFFUSION
  iconvp = vcopt_w%iconv
  imasac = 1
  idiffp = vcopt_w%idiff
  nswrgp = vcopt_w%nswrgr
  imligp = vcopt_w%imligr
  ircflp = vcopt_w%ircflu
  ischcp = vcopt_w%ischcv
  isstpp = vcopt_w%isstpc
  iwarnp = vcopt_w%iwarni
  blencp = vcopt_w%blencv
  epsrgp = vcopt_w%epsrgr
  climgp = vcopt_w%climgr
  extrap = vcopt_w%extrag
  relaxp = vcopt_w%relaxv
  thetap = vcopt_w%thetav
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
 ( idtvar , ivarfl(ivar)    , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                            &
   iwarnp , imucpp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,                   &
   cvara_omg       , cvara_omg       ,                                     &
   coefa_o , coefb_o , coefaf_o , coefbf_o ,                               &
   imasfl , bmasfl ,                                                       &
   viscf  , viscb  , rvoid  , rvoid  ,                                     &
   rvoid  , rvoid  ,                                                       &
   icvflb , ivoid  ,                                                       &
   w6     )

  if (vcopt_w%iwarni.ge.2) then
    rnorm = sqrt(cs_gdot(ncel,smbrw,smbrw))
    write(nfecra,1100) chaine(1:8) ,rnorm
  endif

  do iel = 1,ncel
    smbrk(iel) = smbrk(iel) + w5(iel)
    smbrw(iel) = smbrw(iel) + w6(iel)
  enddo

endif

!===============================================================================
! 11.3 k-omega coupling (ikecou == 1)

!===============================================================================

!  Ordre 2 non pris en compte
if (ikecou.eq.1) then

  ! Take into account, if necessary, the Spalart-Shur rotation/curvature
  ! correction of the production term
  if (irccor.eq.2) then
    do iel = 1, ncel
      w1(iel) = rotfct(iel)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = 1.d0
    enddo
  endif

  do iel = 1, ncel

    rho = crom(iel)

    ! RESOLUTION COUPLEE

    romvsd     = 1.d0/(rho*volume(iel))
    smbrk(iel) = smbrk(iel)*romvsd
    smbrw(iel) = smbrw(iel)*romvsd
    divp23     = d2s3*max(cpro_divukw(iel),zero)
    produc     = w1(iel)*cpro_s2kw(iel)+w2(iel)
    xk         = cvara_k(iel)
    xw         = cvara_omg(iel)
    xxf1       = xf1(iel)
    xgamma     = xxf1*ckwgm1 + (1.d0-xxf1)*ckwgm2
    xbeta      = xxf1*ckwbt1 + (1.d0-xxf1)*ckwbt2

    a11 = 1.d0/dt(iel)                                            &
         - 1.d0/xw*min(produc,zero)+divp23+cmu*xw
    a12 = cmu*xk
    a21 = 0.d0
    a22 = 1.d0/dt(iel)+xgamma*divp23+2.d0*xbeta*xw

    unsdet = 1.d0/(a11*a22 -a12*a21)

    deltk = ( a22*smbrk(iel) -a12*smbrw(iel) )*unsdet
    deltw = (-a21*smbrk(iel) +a11*smbrw(iel) )*unsdet

    ! NOUVEAU TERME SOURCE POUR CODITS

    romvsd = rho*volume(iel)/dt(iel)

    smbrk(iel) = romvsd*deltk
    smbrw(iel) = romvsd*deltw

  enddo


  ! on enleve la convection/diffusion au temps n a SMBRK et SMBRW
  ! s'ils ont ete calcules
  do iel = 1, ncel
    smbrk(iel) = smbrk(iel) - w5(iel)
    smbrw(iel) = smbrw(iel) - w6(iel)
  enddo

  ! Free memory
  deallocate(w5, w6)
  deallocate(w7)

endif

!===============================================================================
! 14. Solving
!===============================================================================

! ---> turbulent kinetic (k) energy treatment
ivar = ik

! Face viscosity
if (vcopt_k%idiff.ge. 1) then

  do iel = 1, ncel
    xxf1 = xf1(iel)
    sigma = xxf1*ckwsk1 + (1.d0-xxf1)*ckwsk2
    w1(iel) = viscl(iel)                                          &
             + vcopt_k%idifft*cvisct(iel)/sigma
  enddo
  call viscfa &
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Solving k
iconvp = vcopt_k%iconv
idiffp = vcopt_k%idiff
ndircp = vcopt_k%ndircl
nswrsp = vcopt_k%nswrsm
nswrgp = vcopt_k%nswrgr
imligp = vcopt_k%imligr
ircflp = vcopt_k%ircflu
ischcp = vcopt_k%ischcv
isstpp = vcopt_k%isstpc
iescap = 0
imucpp = 0
idftnp = ISOTROPIC_DIFFUSION
iswdyp = vcopt_k%iswdyn
iwarnp = vcopt_k%iwarni
blencp = vcopt_k%blencv
epsilp = vcopt_k%epsilo
epsrsp = vcopt_k%epsrsm
epsrgp = vcopt_k%epsrgr
climgp = vcopt_k%climgr
extrap = vcopt_k%extrag
relaxp = vcopt_k%relaxv
thetap = vcopt_k%thetav
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0
init   = 1

call codits &
 ( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp , normp  ,                                              &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_k         , cvara_k         ,                            &
   coefa_k , coefb_k , coefaf_k , coefbf_k ,                      &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinstk , smbrk  , cvar_k          , dpvar  ,                   &
   rvoid  , rvoid  )

! ---> Omega treatment
ivar = iomg

! Face viscosity
if (vcopt_w%idiff.ge. 1) then
  do iel = 1, ncel
    xxf1 = xf1(iel)
    sigma = xxf1*ckwsw1 + (1.d0-xxf1)*ckwsw2
    w1(iel) = viscl(iel)                                          &
                        + vcopt_w%idifft*cvisct(iel)/sigma
  enddo
  call viscfa &
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! Solving omega
iconvp = vcopt_w%iconv
idiffp = vcopt_w%idiff
ndircp = vcopt_w%ndircl
nswrsp = vcopt_w%nswrsm
nswrgp = vcopt_w%nswrgr
imligp = vcopt_w%imligr
ircflp = vcopt_w%ircflu
ischcp = vcopt_w%ischcv
isstpp = vcopt_w%isstpc
iescap = 0
imucpp = 0
idftnp = ISOTROPIC_DIFFUSION
iswdyp = vcopt_w%iswdyn
iwarnp = vcopt_w%iwarni
blencp = vcopt_w%blencv
epsilp = vcopt_w%epsilo
epsrsp = vcopt_w%epsrsm
epsrgp = vcopt_w%epsrgr
climgp = vcopt_w%climgr
extrap = vcopt_w%extrag
relaxp = vcopt_w%relaxv
thetap = vcopt_w%thetav
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0
init   = 1

call codits &
 ( idtvar , init   , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp , normp  ,                                              &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_omg       , cvara_omg       ,                            &
   coefa_o , coefb_o , coefaf_o , coefbf_o ,                      &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinstw , smbrw  , cvar_omg        , dpvar  ,                   &
   rvoid  , rvoid  )

!===============================================================================
! 15. Clipping
!===============================================================================

! Calcul des Min/Max avant clipping, pour affichage
do ii = 1, 2
  if (ii.eq.1) then
    cvar_var => cvar_k
  elseif (ii.eq.2) then
    cvar_var => cvar_omg
  endif

  vrmin(ii) =  grand
  vrmax(ii) = -grand
  do iel = 1, ncel
    var = cvar_var(iel)
    vrmin(ii) = min(vrmin(ii),var)
    vrmax(ii) = max(vrmax(ii),var)
  enddo
enddo

! On clippe simplement k et omega par valeur absolue
call field_get_key_id("clipping_id", kclipp)

call field_get_key_int(ivarfl(ik), kclipp, clip_k_id)
if (clip_k_id.ge.0) then
  call field_get_val_s(clip_k_id, cpro_k_clipped)
endif

call field_get_key_int(ivarfl(iomg), kclipp, clip_w_id)
if (clip_w_id.ge.0) then
  call field_get_val_s(clip_w_id, cpro_w_clipped)
endif

do iel = 1, ncel
  if (clip_k_id.ge.0) &
    cpro_k_clipped(iel) = 0.d0
  if (clip_w_id.ge.0) &
    cpro_w_clipped(iel) = 0.d0
enddo

iclipk(1) = 0
iclipw = 0
do iel = 1, ncel
  xk = cvar_k(iel)
  xw = cvar_omg(iel)
  if (abs(xk).le.epz2) then
    iclipk(1) = iclipk(1) + 1
    if (clip_k_id.ge.0) &
      cpro_k_clipped(iel) = epz2 - xk
    cvar_k(iel) = epz2
  elseif (xk.le.0.d0) then
    iclipk(1) = iclipk(1) + 1
    if (clip_k_id.ge.0) &
      cpro_k_clipped(iel) = - xk
    cvar_k(iel) = -xk
  endif
  if (abs(xw).le.epz2) then
    iclipw = iclipw + 1
    if (clip_w_id.ge.0) &
      cpro_w_clipped(iel) = epz2 - xw
    cvar_omg(iel) = epz2
  elseif (xw.le.0.d0) then
    iclipw = iclipw + 1
    if (clip_w_id.ge.0) &
      cpro_w_clipped(iel) = - xw
    cvar_omg(iel) = -xw
  endif
enddo

! ---  Stockage nb de clippings pour log

call log_iteration_clipping_field(ivarfl(ik), iclipk(1), 0,    &
                                  vrmin(1:1), vrmax(1:1),iclipk(1), iclpkmx(1))
call log_iteration_clipping_field(ivarfl(iomg), iclipw, 0,  &
                                  vrmin(2:2), vrmax(2:2),iclipk(1), iclpkmx(1))

!===============================================================================
! 16. Advanced reinit
!===============================================================================

! Automatic reinitialization at the end of the first iteration:
! wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa,
! then y so we have an idea of the wall distance in complex geometries.
! Then U is initialized with a Reichard layer,
! Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15
! k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cmu)
! The blending function is chosen so that the asymptotic behavior
! and give the correct peak of k (not the same blending than for the EBRSM
! because k profile is not the same for k-omega).
! For omega, far from the wall we take eps/Cmu/k, but next to the wall,
! omega solution is enforced.

!TODO FIXME: why not just before? Are the BC uncompatible?
if (ntcabs.eq.1.and.reinit_turb.eq.1) then
  call field_get_val_v(ivarfl(iu), vel)

  utaurf = 0.05d0*uref
  nu0 = viscl0 / ro0

  do iel = 1, ncel
    ! Compute the velocity magnitude
    xunorm = vel(1,iel)**2 + vel(2,iel)**2 + vel(3,iel)**2
    xunorm = sqrt(xunorm)

    ya = w_dist(iel)
    ypa = ya*utaurf/nu0
    ! Velocity magnitude is imposed (limitted only), the direction is
    ! conserved
    if (xunorm.le.1.d-12*uref) then
      limiter = 1.d0
    else
      limiter = min(utaurf/xunorm*(2.5d0*dlog(1.d0+0.4d0*ypa)            &
      +7.8d0*(1.d0-dexp(-ypa/11.d0)          &
      -(ypa/11.d0)*dexp(-0.33d0*ypa))),      &
      1.d0)
    endif

    vel(1,iel) = limiter*vel(1,iel)
    vel(2,iel) = limiter*vel(2,iel)
    vel(3,iel) = limiter*vel(3,iel)

    ut2 = 0.05d0*uref
    xeps = utaurf**3*min(1.d0/(xkappa*15.d0*nu0/utaurf), &
    1.d0/(xkappa*ya))
    cvar_k(iel) = xeps/2.d0/nu0*ya**2                    &
    * exp(-ypa/25.d0)**2                        &
    + ut2**2/sqrt(cmu)*(1.d0-exp(-ypa/25.d0))**2

    cvar_omg(iel) = ut2**3/(xkappa*15.d0*nu0/ut2)/(ut2**2/sqrt(cmu))/cmu

  enddo
endif

! Free memory
deallocate(viscf, viscb)
deallocate(smbrk, smbrw)
deallocate(tinstk, tinstw, xf1)
deallocate(w1, w2, usimpk, usimpw)
deallocate(dpvar)
deallocate(prodk, prodw)
if (allocated(gradv)) deallocate(gradv)
if (allocated(s2pw2)) deallocate(s2pw2)

if (allocated(rotfct))  deallocate(rotfct)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/, &
'   ** RESOLUTION DU K-OMEGA'                     ,/,&
'      ---------------------'                     ,/)
 1100 format(1X,A8,' : BILAN EXPLICITE = ',E14.5)

#else

 1000 format(/, &
'   ** SOLVING K-OMEGA'                           ,/,&
'      ---------------'                           ,/)
 1100 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)

#endif

!----
! End
!----

return

end subroutine

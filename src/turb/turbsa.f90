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
! Purpose:
! --------
!> \file turbsa.f90
!> \brief Solving op the equation of \f$ \tilde{\nu} \f$, which is the scalar
!>        quantity defined by the Spalart-Allmaras model for 1 time-step.
!>
!> Please refer to the
!> <a href="../../theory.pdf#spalart"><b>Spalart-Allmaras model</b></a>
!> section of the theory guide for more informations.
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
!> \param[in]     icepdc        number of ncepdp cells with head losses
!> \param[in]     icetsm        number of cells with mass source
!> \param[in]     itypsm        type of mass source for the variables
!>                               (cf. cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     ckupdc        work array for head losses
!> \param[in]     smacel        value of variables associated to the
!                               mass source
!                               for ivar = ipr, smacel = mass flux
!> \param[in]     itypfb        boundary face types
!______________________________________________________________________________!

subroutine turbsa &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel ,                                              &
   itypfb )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstnum
use cstphy
use optcal
use mesh
use field
use field_operator
use parall
use pointe, only: dispar
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          itypfb(nfabor)
integer          icvflb
integer          ivoid(1)

double precision dt(ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

integer          iel   , ifac  , inc   , iccocg, iprev, ivar
integer          iiun
integer          nswrgp, imligp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          iflmas, iflmab
integer          iwarnp
integer          istprv
integer          ipatrg
integer          imucpp, idftnp, iswdyp

double precision romvsd
double precision visct , rom
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision thetv, thetp1, thetap
double precision tuexpn
double precision cofbnu
double precision chi  , chi3, taussa, nusa, distbf, fw, fv1, fv2
double precision gsa , rsa , dsigma, cv13
double precision nu0, dsa0, hssa, omega, sbar, cst2, cst3

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: tsimp
double precision, allocatable, dimension(:) :: rhssa, tinssa, trgrdu
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: trgrdn, vort
double precision, allocatable, dimension(:) :: tsexp
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: csab1r, rotfct
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, cromo
double precision, dimension(:), pointer :: cvar_nusa, cvara_nusa
double precision, dimension(:), pointer :: cpro_pcvto
double precision, dimension(:), pointer :: viscl, cvisct
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: c_st_nusa_p

type(var_cal_opt) :: vcopt_nusa

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays for the turbulence resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(tsimp(ncelet))
allocate(trgrdn(ncelet), vort(ncelet))
allocate(rhssa(ncelet))
allocate(tinssa(ncelet), trgrdu(ncelet))

! Allocate work arrays
allocate(w1(ncelet))
allocate(tsexp(ncelet))
allocate(dpvar(ncelet))

call field_get_val_s(icrom, crom)
call field_get_val_s(ivisct, cvisct)
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ibrom, brom)
call field_get_val_s(ivarfl(inusa), cvar_nusa)
call field_get_val_prev_s(ivarfl(inusa), cvara_nusa)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

ivar   = inusa

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt_nusa)

thetv  = vcopt_nusa%thetav

call field_get_key_int(ivarfl(inusa), kstprv, istprv)
if (istprv.ge.0) then
  call field_get_val_s(istprv, c_st_nusa_p)
  istprv = 1
endif

call field_get_val_s(icrom, cromo)

if (istprv.ge.0) then
  if (iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
  endif
  if (iviext.gt.0) then
    call field_get_val_prev_s(ivisct, cpro_pcvto)
  endif
endif

if (vcopt_nusa%iwarni.ge.1) then
  write(nfecra,1000)
endif

! Calculation of some constants
dsigma = 1.d0 / csasig
cv13 = csav1**3

! To avoid numerical problem, constant used to prevent taussa from
! being negative (see Oliver TA 2008)
cst2 = 0.7d0
cst3 = 0.9d0

!===============================================================================
! 2. Compute the vorticity omega, the trace of the velocity gradient
!    and the gradient of nusa
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3, 3, ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,       &
                           gradv)

! vort = omega**2 = dudy**2 + dvdx**2 + dudz**2 + dwdx**2 + dvdz**2 + dwdy**2
!                - 2*dudy*dvdx - 2*dudz*dwdx - 2*dvdz*dwdy
!
!        = 2 Oij.Oij
! trgrdu = dudx + dvdy + dwdz

do iel = 1, ncel
  vort(iel) = (gradv(2, 1, iel) - gradv(1, 2, iel))**2   &
            + (gradv(3, 1, iel) - gradv(1, 3, iel))**2   &
            + (gradv(3, 2, iel) - gradv(2, 3, iel))**2
  trgrdu(iel) = gradv(1, 1, iel) + gradv(2, 2, iel) + gradv(3, 3, iel)
enddo

! Free memory
deallocate(gradv)

! Allocate a temporary array for the gradient calculation
allocate(grad(3,ncelet))

! Compute the gradient of nusa

iccocg = 1

call field_gradient_scalar(ivarfl(inusa), iprev, imrgra, inc,       &
                           iccocg,                                  &
                           grad)

! trgrdn = GRAD(nusa)**2
do iel = 1, ncel
  trgrdn(iel) = grad(1,iel)**2 + grad(2,iel)**2 + grad(3,iel)**2
enddo

! Free memory
deallocate(grad)

!===============================================================================
! 3. Compute the buoyant term
!===============================================================================

! Gravity is not taken into account at the moment

!===============================================================================
! 4. Source terms are finalized

!      stored in rhssa
!===============================================================================

! Herebelow, we only handle  the case where all the walls have
! the same roughness
! To extend it, we should be able to link every fluid cell to a boundary face
! (and then give it the appropriate roughness value)

ipatrg = 0
dsa0 = -999.d0
hssa = -999.d0

call field_get_coefb_s(ivarfl(inusa), coefbp)

do ifac = 1, nfabor
  if (itypfb(ifac).eq.iparug) then
    ipatrg = 1
    cofbnu = coefbp(ifac)
    ! Roughness of the wall
    dsa0   = distb(ifac) *cofbnu/(1.d0-cofbnu)
    hssa   = exp(8.5d0*xkappa)*dsa0
  endif
  if (ipatrg.ne.0) exit
enddo

if (irangp.ge.0) then
  call parcpt(ipatrg)
  if (ipatrg.ne.0) then
    call parsom(dsa0)
    dsa0=dsa0/ipatrg
  endif
endif

! Take into account the Spalart-Shur rotation/curvature correction, if necessary
! => variable production term coefficient (csab1)
allocate(csab1r(ncel))

if (irccor.eq.1) then

  ! Allocate temporary array for rotation function
  allocate(rotfct(ncel))

  ! Compute the rotation function (w1 array not used)
  call rotcor(dt, rotfct, w1)

  do iel = 1, ncel
    csab1r(iel) = csab1*rotfct(iel)
  enddo

  ! Free memory
  deallocate(rotfct)

else
  do iel = 1, ncel
    csab1r(iel) = csab1
  enddo
endif

! If source terms are extrapolated, rho is rho^n
!                                 visct is visct^n
do iel = 1, ncel

  if (istprv.ge.0 .and. iviext.gt.0) then
    visct = cpro_pcvto(iel)
  else
    visct = cvisct(iel)
  endif
  rom   = cromo(iel)
  ! kinematic viscosity
  nu0   = viscl(iel)/rom
  ! We have to know if there is any rough wall
  distbf= dispar(iel)
  ! viscosity of SA
  nusa  = cvara_nusa(iel)
  chi   = nusa/nu0
  ! If we have a rough wall
  if (ipatrg.ne.0) then
    distbf = distbf + dsa0
    chi  = chi + 0.5d0* hssa/distbf
  endif
  chi3  = chi**3
  fv1   = chi3/(chi3 + cv13 )
  fv2   = 1.d0 - nusa /(nu0 + nusa*fv1)

  ! Numerical fix to prevent taussa to be smaller than 0
  ! (reported in Oliver T.A. 2008)
  sbar = nusa/(xkappa*distbf)**2*fv2
  omega = sqrt(vort(iel))

  if (sbar.ge.-cst2*omega) then
    taussa = omega+sbar
  else
    taussa = omega*(1.d0 + &
                   (cst2**2*omega+cst3*sbar)/((cst3-2.d0*cst2)*omega-sbar))
  endif

  ! Computation of fw
  if (nusa.ge.10.d0*taussa*(xkappa*distbf)**2) then
    rsa = 10.d0
  else
    rsa   = nusa/(taussa*(xkappa*distbf)**2)
  endif
  gsa   = rsa + csaw2*(rsa**6-rsa)
  fw    = gsa*( (1.d0+csaw3**6)/(gsa**6+csaw3**6))**(1.d0/6.d0)

  rhssa(iel) = volume(iel)*rom*(                                 &
     dsigma * csab2*trgrdn(iel)+csab1r(iel)*taussa*nusa-csaw1*fw*(nusa/distbf)**2)

  ! Implicitation of the negative source terms of the SA equation.
  ! NB: this term could be negative, and if so, then we explicit it.
  tinssa(iel) = (max(csaw1*fw*nusa/distbf**2-csab1r(iel)*taussa,0.d0)         &
                      )*rom*volume(iel)

enddo

! Free memory
deallocate(csab1r)

!===============================================================================
! 5. Take user source terms into account

!      omega**2 = vort and the trace of the velocity gradient = trgrdu
!        are available
!      The explicit part is stored in    tsexp
!      The implicit part is stored in    tsimp
!===============================================================================
do iel = 1, ncel
  tsimp(iel) = 0.d0
  tsexp (iel) = 0.d0
enddo

call cs_user_turbulence_source_terms &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(inusa)   ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   tsexp  , tsimp )

!===============================================================================
! 6. User source terms and d/dt(rho) and div(rho u) are taken into account
!      stored in rhssa
!===============================================================================

! If source terms are extrapolated
if (istprv.ge.0) then

  do iel = 1, ncel

     ! Ts^(n-1) (Term User EXPlicit Nusa)
     tuexpn =c_st_nusa_p(iel)

    ! The explicit user source terms are stored for the next time step
    ! We store the explicit source terms of time n (model source terms + user
    ! source terms)
    c_st_nusa_p(iel) = rhssa(iel) + tsexp(iel)


    ! --- Extrapolated explicit source terms
    rhssa(iel) = - thetst*tuexpn

    rhssa(iel) = tsimp(iel)*cvara_nusa(iel) + rhssa(iel)

    ! --- Implicit user source terms
    ! Here it is assumed that -tsimp > 0. That is why it is implicited
    tinssa(iel) = tinssa(iel) - tsimp(iel)*thetv

  enddo

! If source terms are not extrapolated, then they are directly added to the RHS
else
  do iel = 1, ncel
    rhssa(iel) = rhssa(iel) + tsimp(iel)*cvara_nusa(iel) + tsexp(iel)

    ! --- Implicit user source terms
    tinssa(iel) = tinssa(iel) + max(-tsimp(iel),zero)
  enddo
endif

! --- rho/dt and div(rho u)
!     Extrapolated or not in coherence with bilsc2
do iel = 1, ncel
  rom = crom(iel)
  romvsd = rom*volume(iel)/dt(iel)

  ! tinssa already contains the negativ implicited source term
  tinssa(iel) = tinssa(iel)                                        &
               +vcopt_nusa%istat*romvsd
enddo


!===============================================================================
! 7. Lagrangian source terms (Explicit part)
!===============================================================================

! Not accounted for at the moment.

!===============================================================================
! 8. Explicit mass source terms

!    Gamma.var_prev is stored in w1
!===============================================================================

if (ncesmp.gt.0) then

  ! Integer equal to 1. (in navsto: nb of sub-iteration)
  iiun = 1

  ! --- Explicit and Implicit part
  !     -Gamma.var_prev is added to the RHS and Gamma to tinssa
  ivar = inusa

  call catsma &
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
   isto2t ,                                                       &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , cvara_nusa     , smacel(1,ivar) , smacel(1,ipr) ,     &
   rhssa  , tinssa , w1 )

  ! --- Explicit part: Gamma Pinj
  !     (if we extrapolate source terms, Gamma.var_prev is stored in prev. TS)
  if (istprv.ge.0) then
    do iel = 1, ncel
      c_st_nusa_p(iel) = c_st_nusa_p(iel) + w1(iel)
    enddo
  else
    do iel = 1, ncel
      rhssa(iel) = rhssa(iel) + w1(iel)
    enddo
  endif

endif

! Finalization of the extrapolated explicit source terms
if (istprv.ge.0) then
  thetp1 = 1.d0 + thetst
  do iel = 1, ncel
    rhssa(iel) = rhssa(iel) + thetp1 * c_st_nusa_p(iel)
  enddo
endif

!===============================================================================
! 9. Solving of the transport equation on nusa
!===============================================================================

ivar = inusa

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Face viscosity

if (vcopt_nusa%idiff.ge.1) then

  do iel = 1, ncel
    rom = crom(iel)

    ! diffusibility: 1/sigma*(mu_laminaire+ rho*nusa)
    w1(iel) = dsigma *( viscl(iel)                                &
                        + vcopt_nusa%idifft*cvara_nusa(iel)*rom )
  enddo

  call viscfa                                                     &
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

! --- Solving

iconvp = vcopt_nusa%iconv
idiffp = vcopt_nusa%idiff
ndircp = ndircl(ivar)
nswrsp = vcopt_nusa%nswrsm
nswrgp = vcopt_nusa%nswrgr
imligp = vcopt_nusa%imligr
ircflp = vcopt_nusa%ircflu
ischcp = vcopt_nusa%ischcv
isstpp = vcopt_nusa%isstpc
iescap = 0
imucpp = 0
idftnp = vcopt_nusa%idften
iswdyp = vcopt_nusa%iswdyn
iwarnp = vcopt_nusa%iwarni
blencp = vcopt_nusa%blencv
epsilp = vcopt_nusa%epsilo
epsrsp = vcopt_nusa%epsrsm
epsrgp = vcopt_nusa%epsrgr
climgp = vcopt_nusa%climgr
extrap = vcopt_nusa%extrag
relaxp = vcopt_nusa%relaxv
thetap = vcopt_nusa%thetav
! all boundary convective flux with upwind
icvflb = 0

call codits &
 ( idtvar , ivarfl(ivar)    , iconvp , idiffp , ndircp ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   cvara_nusa      , cvara_nusa      ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinssa , rhssa  , cvar_nusa       , dpvar ,                    &
   rvoid  , rvoid  )

!===============================================================================
! 10. Clipping
!===============================================================================

call clipsa(ncel)

! Free memory
deallocate(viscf, viscb)
deallocate(tsimp)
deallocate(rhssa)
deallocate(tinssa, trgrdu)
deallocate(trgrdn, vort)
deallocate(w1)
deallocate(tsexp)
deallocate(dpvar)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** Resolution de Spalart-Allmaras                         ',/,&
'      ------------------------------------                   ',/)
#else

 1000 format(/,                                                   &
'   ** Solving Spalart-Allmaras      '                         ,/,&
'      ------------------------------'                         ,/)
#endif

!----
! End
!----

return

end subroutine

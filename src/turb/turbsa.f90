!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine turbsa &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   ckupdc , smacel ,                                              &
   itypfb )

!===============================================================================
! Purpose:
! --------

! Solving op the equation of nusa, which is the scalar quantity defined by
! the Spalart-Allmaras model for 1 time-step.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use mesh
use field
use field_operator
use parall
use pointe, only: dispar

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

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

integer          iel   , ifac  , inc   , iccocg, iprev, ivar
integer          iiun
integer          nswrgp, imligp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ipcvst, ipcvis, iflmas, iflmab
integer          iwarnp, ipp
integer          iptsta
integer          ipcvto, ipcvlo
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
double precision surfn, nu0, dsa0, hssa, omega, sbar, cst2, cst3

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
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp

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
ipcvst = ipproc(ivisct)
ipcvis = ipproc(iviscl)
call field_get_val_s(ibrom, brom)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

ivar   = inusa
thetv  = thetav(ivar)

call field_get_val_s(icrom, cromo)
ipcvto = ipcvst
ipcvlo = ipcvis

if(isto2t.gt.0) then
  if (iroext.gt.0) then
    call field_get_val_prev_s(icrom, cromo)
  endif
  if(iviext.gt.0) then
    ipcvto = ipproc(ivista)
    ipcvlo = ipproc(ivisla)
  endif
endif

! If source terms are extrapolated
if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
else
  iptsta = 0
endif

if(iwarni(inusa).ge.1) then
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

iccocg = 1
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

! Herebelow, we only handle  the case where all the walls have the same roughness
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

if(irangp.ge.0) then
  call parcpt(ipatrg)
  if(ipatrg.ne.0) then
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
  call rotcor(dt, rtpa, rotfct, w1)
  !==========

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

  visct = propce(iel,ipcvto)
  rom   = cromo(iel)
  ! Kinematic viscosity
  nu0   = propce(iel,ipcvis)/rom
  ! We have to know if there is any rough wall
  distbf= dispar(iel)
  ! viscosity of SA
  nusa  = rtpa(iel,inusa)
  chi   = nusa/nu0
  ! If we have a rough wall
  if(ipatrg.ne.0) then
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
  ! NB : this term could be negative, and if so, then we explicit it.
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

call ustssa                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , vort   , trgrdu ,                            &
   tsexp  , tsimp )

!===============================================================================
! 6. User source terms and d/dt(rho) and div(rho u) are taken into account

!      stored in rhssa
!===============================================================================

! If source terms are extrapolated
if (isto2t.gt.0) then

  do iel = 1, ncel

     ! Ts^(n-1) (Term User EXPlicit Nusa)
     tuexpn =propce(iel,iptsta)

    ! The explicit user source terms are stored for the next time step
    ! On stoque les TS explicites du temps n (TS model + TS utilisateur)
    propce(iel,iptsta) = rhssa(iel) + tsexp(iel)


    ! --- Extrapolated explicit source terms
    rhssa(iel) = - thetst*tuexpn

    rhssa(iel) = tsimp(iel)*rtpa(iel,inusa) + rhssa(iel)

    ! --- Implicit user source terms
    ! Here it is assumed that -tsimp > 0. That is why it is implicited
    tinssa(iel) = tinssa(iel) - tsimp(iel)*thetv

  enddo

! If source terms are not extrapolated, then they are directly added to the RHS
else
  do iel = 1, ncel
    rhssa(iel) = rhssa(iel) + tsimp(iel)*rtpa(iel,inusa) + tsexp(iel)

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
               +istat(inusa)*romvsd
enddo


!===============================================================================
! 7. Lagrangian source terms (Explicit part)
!===============================================================================

! Not accounted for at the moment.

!===============================================================================
! 8. Explicit mass source terms

!    Gamma*RTPAi is stored in w1
!===============================================================================

if (ncesmp.gt.0) then

  ! Integer equal to 1. (in navsto: nb of sub-iteration)
  iiun = 1

  ! --- Explicit and Implicit part
  !     -Gamma RTPA is added to the RHS and Gamma*theta to tinssa
  ivar = inusa

  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   ,                            &
                              isto2t , thetv ,                    &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipr) ,       &
   rhssa  , tinssa , w1 )

  ! --- Explicit part: Gamma*RTPAi
  !     (if we extrapolate source terms, Gamma*RTPAi is stored in propce)
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta) = propce(iel,iptsta) + w1(iel)
    enddo
  else
    do iel = 1, ncel
      rhssa(iel) = rhssa(iel) + w1(iel)
    enddo
  endif

endif

! Finalization of the extrapolated explicit source terms
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thetst
  do iel = 1, ncel
    rhssa(iel) = rhssa(iel) + thetp1    * propce(iel,iptsta)
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

ipp    = ipprtp(ivar)

! Face viscosity

if (idiff(ivar).ge.1) then

  do iel = 1, ncel
    rom = crom(iel)

    ! diffusibility: 1/sigma*(mu_laminaire+ rho*nusa)
    w1(iel) = dsigma *( propce(iel,ipcvis)                        &
                        + idifft(ivar)*rtpa(iel,inusa)*rom )
  enddo

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

  ! Be carefull with the walls:
  !  If we have a smooth wall then nusa is zero at the wall
  !  If we have a rough wall then nusa_wall*(1- IprF/d0)=Vipr

  do ifac = 1, nfabor

    iel   = ifabor(ifac)
    surfn = surfbn(ifac)

    ! Smooth wall
    if (itypfb(ifac).eq.iparoi) then
      viscb(ifac) = dsigma * propce(iel,ipcvis)*surfn/distb(ifac)

    ! Rough wall
    elseif (itypfb(ifac).eq.iparug) then

      rom = crom(iel)

      ! dsa0 is recomputed in case of many different roughness
      cofbnu = coefbp(ifac)

      ! Roughness of the wall
      dsa0   = distb(ifac) *cofbnu/(1.d0-cofbnu)
      hssa   = exp(8.5d0*xkappa)*dsa0

      ! For rough walls: nusa_F*(IprF/d0+1) = nusa_Ipr
      viscb(ifac) = dsigma * ( propce(iel,ipcvis)                    &
                   + idifft(ivar)*rtpa(iel,inusa)*rom                &
                   * dsa0/(distb(ifac)+dsa0)            )*surfn/distb(ifac)

    endif

  enddo

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

! --- Solving

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
ndircp = ndircl(ivar)
nitmap = nitmax(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
iescap = 0
imucpp = 0
idftnp = idften(ivar)
iswdyp = iswdyn(ivar)
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetap = thetav(ivar)
! all boundary convective flux with upwind
icvflb = 0

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   tinssa , rhssa  , rtp(1,ivar)     , dpvar ,                    &
   rvoid  , rvoid  )

!===============================================================================
! 10. Clipping
!===============================================================================

call clipsa(ncelet, ncel, nvar, rtp)
!==========

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
'   ** RESOLUTION DE SPALART-ALLMARAS                         ',/,&
'      ------------------------------------                   ',/)
#else

 1000 format(/,                                                   &
'   ** SOLVING SPALART-ALLMARAS      '                         ,/,&
'      ------------------------------'                         ,/)
#endif

!----
! End
!----

return

end subroutine

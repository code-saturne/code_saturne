!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
! Function :
! --------

!> \file driflu.f90
!>
!> \brief Compute the modified convective flux for scalars with a drift.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iflid         index of the current drift scalar field
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp           calculated variables at cell centers
!>                               (at current time step)
!> \param[in]     rtpa          calculated variables at cell centers
!>                               (at previous time step)
!> \param[in]     propce        physical properties at cell centers
!> \param[in,out] imasfl        scalar mass flux at interior face centers
!> \param[in,out] bmasfl        scalar mass flux at boundary face centers
!> \param[in,out] rovsdt        Non stationnary term and mass aggregation term
!> \param[in,out] smbrs         right hand side for the scalar iscal
!______________________________________________________________________________

subroutine driflu &
( iflid  ,                                                       &
  dt     , rtp    , rtpa   , propce ,                            &
  imasfl , bmasfl ,                                              &
  rovsdt , smbrs  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use pointe
use field
use mesh
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          iflid

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision imasfl(nfac), bmasfl(ndimfb)
double precision rovsdt(ncelet), smbrs(ncelet)

! Local variables

integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg
integer          ipcvst, ipcvsl, iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp
integer          ircflp, ischcp, isstpp
integer          ippu  , ippv  , ippw
integer          isou  , jsou
integer          f_id
integer          iflmb0, idftnp, iphydp, ivisep, itypfl
integer          keysca, iscal, keydri, iscdri, icvflb
integer          ivoid(1)

double precision epsrgp, climgp, extrap, blencp
double precision thetap
double precision rhovdt
double precision omega
double precision rho
double precision relaxp

double precision rvoid(1)

character*80     fname

double precision, dimension(:), allocatable :: w1, viscce
double precision, dimension(:), allocatable :: coefap, coefbp
double precision, dimension(:), allocatable :: cofafp, cofbfp
double precision, dimension(:,:), allocatable :: coefa1
double precision, dimension(:,:,:), allocatable :: coefb1
double precision, dimension(:,:), allocatable :: drift, vel, dudt
double precision, dimension(:), allocatable :: viscf, viscb
double precision, dimension(:), allocatable :: flumas, flumab

double precision, dimension(:), pointer :: taup
double precision, dimension(:), pointer :: taufpt
double precision, dimension(:,:), pointer :: coefav, cofafv
double precision, dimension(:,:,:), pointer :: coefbv, cofbfv
double precision, dimension(:), pointer :: imasfl_mix, bmasfl_mix
double precision, dimension(:), pointer :: brom, crom

!===============================================================================

!===============================================================================
! 0. Key words for fields
!===============================================================================

! Key id for scalar id
call field_get_key_id("scalar_id", keysca)
call field_get_key_int(iflid, keysca, iscal)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)
call field_get_key_int(iflid, keydri, iscdri)

! Pointers to the mass fluxes of the mix (based on mix velocity)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl_mix)
call field_get_val_s(iflmab, bmasfl_mix)

!===============================================================================
! 1. Initialization
!===============================================================================

ivar = isca(iscal)

! --- Physical properties
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
ipcvst = ipproc(ivisct)

! --- Brownian diffusivity
if (ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! Name of the scalar
call field_get_name(ivarfl(ivar), fname)

! Index of the corresponding relaxation time (taup)
call field_get_id('drift_tau_'//trim(fname), f_id)
call field_get_val_s(f_id, taup)

! Index of the corresponding interaction time particle--eddies (taufpt)
if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS)) then
  call field_get_id('drift_turb_tau_'//trim(fname), f_id)
  call field_get_val_s(f_id, taufpt)
endif

! Vector containing all the additional convective terms
allocate(drift(3, ncelet))
allocate(vel(3, ncelet))
allocate(dudt(3, ncelet))
allocate(w1(ncelet), viscce(ncelet))
allocate(coefap(ndimfb), coefbp(ndimfb))
allocate(cofafp(ndimfb), cofbfp(ndimfb))
allocate(coefa1(3, ndimfb), coefb1(3, 3, ndimfb))
allocate(viscf(nfac), viscb(nfabor))
allocate(flumas(nfac), flumab(nfabor))

do ifac = 1, nfac
  viscf(ifac) = 0.d0
  flumas(ifac) = 0.d0
enddo

do ifac = 1, nfabor
  viscb(ifac) = 0.d0
  flumab(ifac) = 0.d0
enddo

!===============================================================================
! 2. initialize the additional convective flux with the gravity term
!===============================================================================

do iel = 1, ncel
  rho = crom(iel)
  drift(1, iel) = rho*taup(iel)*gx
  drift(2, iel) = rho*taup(iel)*gy
  drift(3, iel) = rho*taup(iel)*gz
enddo

!===============================================================================
! 3. Computation of the turbophoresis and the thermophoresis terms
!===============================================================================

! Initialized to 0
do iel = 1, ncel
  viscce(iel) = 0.d0
enddo

if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS).and.iturb.ne.0) then

  ! The diagonal part is easy to implicit (Grad (K) . n = (K_j - K_i)/IJ)

  ! Compute the K=1/3*trace(K) coefficient (diffusion of Zaichik)

  if (itytur.eq.3) then

    do iel = 1, ncel

      ! Correction by Omega
      omega = taup(iel)/taufpt(iel)
      ! FIXME: use idifft or not?
      viscce(iel) = 1.d0/3.d0                               &
                  * taup(iel)/(1.d0+omega)*( rtp(iel,ir11)  &
                                           + rtp(iel,ir22)  &
                                           + rtp(iel,ir33) )
    enddo

  elseif (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then

    do iel = 1, ncel

      ! Correction by Omega
      omega = taup(iel)/taufpt(iel)

      viscce(iel) = 2.d0*taup(iel)/(1.d0+omega)*rtp(iel,ik)
    enddo

  else

  endif

endif

if (btest(iscdri, DRIFT_SCALAR_THERMOPHORESIS)) then

  ! propce(iel,ipcvsl): contains the Brownian motion
  !------------------------------------------------

  if (ipcvsl.gt.0) then

    do iel = 1, ncel
      viscce(iel) = viscce(iel) + propce(iel,ipcvsl)
    enddo

  else

    do iel = 1, ncel
      viscce(iel) = viscce(iel) + visls0(iscal)
    enddo

  endif

endif

if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS).or.    &
    btest(iscdri, DRIFT_SCALAR_THERMOPHORESIS)) then

  iphydp = 0
  inc    = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  ! Face diffusivity of 1. to compute (Grad K . n)_face
  do iel = 1, ncelet
    w1(iel) = 1.d0
  enddo

  call viscfa &
  ( imvisf ,            &
    w1     ,            &
    viscf  , viscb  )

  ! Homogeneous Neumann BC
  do ifac = 1, nfabor
    ! BCs for gradients
    coefap(ifac) = 0.d0
    coefbp(ifac) = 1.d0
    ! BCs for fluxes
    cofafp(ifac) = 0.d0
    cofbfp(ifac) = 0.d0
  enddo

  init   = 0

  call itrmas &
  !==========
 ( init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rvoid  ,                                                       &
   viscce ,                                                       &
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   viscf  , viscb  ,                                              &
   w1     , w1     , w1     ,                                     &
   flumas , flumab )

! TODO add extradiagonal part
endif

!===============================================================================
! 4. Centrifugal force (particular derivative Du/Dt)
!===============================================================================

if (btest(iscdri, DRIFT_SCALAR_CENTRIFUGALFORCE)) then

  do iel = 1, ncel
    vel(1, iel) = rtp(iel, iu)
    vel(2, iel) = rtp(iel, iv)
    vel(3, iel) = rtp(iel, iw)

    rhovdt = crom(iel)*volume(iel)/dt(iel)

    dudt(1,iel) = -rhovdt*(rtp(iel,iu)-rtpa(iel,iu))
    dudt(2,iel) = -rhovdt*(rtp(iel,iv)-rtpa(iel,iv))
    dudt(3,iel) = -rhovdt*(rtp(iel,iw)-rtpa(iel,iw))
  enddo

  iconvp = 1
  idiffp = 0
  nswrgp = nswrgr(iu)
  imligp = imligr(iu)
  ircflp = ircflu(iu)
  ischcp = ischcv(iu)
  isstpp = isstpc(iu)
  inc    = 1
  ivisep = 0
  ippu   = ipprtp(iu)
  ippv   = ipprtp(iv)
  ippw   = ipprtp(iw)
  iwarnp = iwarni(iu)
  idftnp = idften(iu)
  blencp = blencv(iu)
  epsrgp = epsrgr(iu)
  climgp = climgr(iu)
  thetap = thetav(iu)
  relaxp = relaxv(iu)
  icvflb = 0

  ! Get Boundary conditions of the velocity
  call field_get_coefa_v (ivarfl(iu), coefav)
  call field_get_coefb_v (ivarfl(iu), coefbv)
  call field_get_coefaf_v(ivarfl(iu), cofafv)
  call field_get_coefbf_v(ivarfl(iu), cofbfv)

  ! Warning: bilsc adds "-( grad(u) . u)"
  call bilscv &
  !==========
 ( idtvar , iu     , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   ippu   , iwarnp , idftnp ,                                     &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   vel    , vel    ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   imasfl_mix , bmasfl_mix ,                                      &
   viscf  , viscb  , rvoid  , rvoid  ,                            &
   icvflb , ivoid  ,                                              &
   dudt   )

  do iel = 1, ncel
    drift(1, iel) = drift(1, iel) + taup(iel)*dudt(1, iel)/volume(iel)
    drift(2, iel) = drift(2, iel) + taup(iel)*dudt(2, iel)/volume(iel)
    drift(3, iel) = drift(3, iel) + taup(iel)*dudt(3, iel)/volume(iel)
  enddo

endif

!===============================================================================
! 5. Electrophoresis term
!===============================================================================

if (btest(iscdri, DRIFT_SCALAR_ELECTROPHORESIS)) then

  !TODO
  call csexit(1)

endif

!===============================================================================
! 6. Finalization of the mass flux
!===============================================================================

! Zero additional flux at the boundary
do ifac = 1, nfabor

  do isou = 1, 3
    coefa1(isou, ifac) = 0.d0
    do jsou = 1, 3
      if (isou.eq.jsou) then
        coefb1(isou, jsou, ifac) = 0.d0 !FIXME Dirichlet or Neumann ?
      else
        coefb1(isou, jsou, ifac) = 0.d0
      endif
    enddo
  enddo

enddo

init   = 0
inc    = 1
iflmb0 = 0
itypfl = 0 ! drift has already been multiplied by rho
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

call inimav &
!==========
 ( ivar   , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   drift  ,                                                       &
   coefa1 , coefb1 ,                                              &
   flumas , flumab )

do ifac = 1, nfac
  imasfl(ifac) = imasfl_mix(ifac) + flumas(ifac)
enddo

do ifac = 1, nfabor
  bmasfl(ifac) = bmasfl_mix(ifac) + flumab(ifac)
enddo

!===============================================================================
! 7. Mass aggregation term of the additional part "div(rho(u_p-u_f))"
!===============================================================================

init = 1
iconvp = iconv(ivar)
thetap = thetav(ivar)

call divmas &
 ( ncelet , ncel   , nfac   , nfabor , init   , nfecra ,    &
   ifacel , ifabor ,                                        &
   flumas , flumab ,                                        &
   w1     )

! NB: if the porosity module is swiched on, the the porosity is already
! taken into account in w1

! --> mass aggregation term
do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) + iconvp*thetap*w1(iel)
  smbrs(iel) = smbrs(iel) - iconvp*w1(iel)*rtpa(iel,ivar)
enddo

! Free memory
deallocate(viscce)
deallocate(drift)
deallocate(vel)
deallocate(dudt)
deallocate(w1)
deallocate(viscf, viscb)
deallocate(flumas, flumab)
deallocate(coefap, coefbp)
deallocate(cofafp, cofbfp)
deallocate(coefa1, coefb1)

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine

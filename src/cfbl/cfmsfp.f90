!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine cfmsfp &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   flumas , flumab ,                                              &
   trflms , trflmb )

!===============================================================================
! FONCTION :
! ----------

!  "MASS FLUX" AT THE FACES CALCULATION FOR THE CFL RESTRICTION CALCULATION
!   AND THE SOLVING OF THE PRESSURE

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itspdv           ! e  ! <-- ! calcul termes sources prod et dissip           !
!                  !    !     !  (0 : non , 1 : oui)                           !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! flumas(nfac)     ! tr ! --> ! flux de masse aux faces internes               !
! flumab(nfabor    ! tr ! --> ! flux de masse aux faces de bord                !
! trflms(nfac)     ! tr ! --- ! tableau de travail                             !
! trflmb(nfabor    ! tr ! --- ! tableau de travail                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
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
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision flumas(nfac), flumab(nfabor)
double precision trflms(nfac), trflmb(nfabor)

! Local variables

integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, ii, jj
integer          nswrgp, imligp, iwarnp

integer          ivar0 , isou
integer          imaspe, iflmb0, itypfl
integer          icliup, iclivp, icliwp, iclvar
integer          itsqdm, iiun  , iextts

double precision epsrgp, climgp, extrap
double precision flui  , fluj  , pfac  , thetv, rom

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10, w11, w12
double precision, dimension(:), pointer :: crom, brom

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(w10(ncelet), w11(ncelet), w12(ncelet))

! --- Number of the computational variable
!     Pressure
ivar     = ipr

! Density

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

! ---> Mass flux initialization

do ifac = 1, nfac
  flumas(ifac) = 0.d0
enddo
do ifac = 1, nfabor
  flumab(ifac) = 0.d0
enddo

!===============================================================================
! 2. MASS FLUX AT THE FACES
!===============================================================================

!     2.1 SOURCE TERMS OF THE MOMENTUM EQUATIONS
!     ==========================================

!     FX -> W5 , FY -> W6 , FZ -> W7

!     Some first tests (double expansion waves in a shock tube)
!       has shown that taking into account all the
!       momentum equation terms in the mass equation seems to be a
!       bad idea (in particular the convective term, but the diffusive
!       term, the transposed gradient, the mass and user source terms
!       were all null in the considered tests).
!       However, it may be due to a bug at that early stage of implementation
!       of the algorithm (but we didn't find it).
!     We thus recommand not to take into account the momentum source terms,
!       except the gravity term (because it is in balance with the pressure
!       gradient and because its effect is visible at equilibrium).
!     However, we keep here the implementation of the preliminary tests
!       (1.1.0.h version) with an overall test so that the correction is not
!       active (thus, there is no user question and there is always the
!       possibility to perform other tests in the future).
!     Note that, with these terms, the thoeretical analysis is harder
!     (Without these terms we are in the configuration Euler + gravity)

! --- Initialization
do iel = 1, ncel
  w10(iel) = 0.d0
enddo
do iel = 1, ncel
  w11(iel) = 0.d0
enddo
do iel = 1, ncel
  w12(iel) = 0.d0
enddo


!     Test on momentum source terms
itsqdm = 0
if (itsqdm.ne.0) then

  ! --- User soure term
  call ustsnv &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iu  ,                                                          &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smacel , w10    , w9     ) !FIXME


  ! --- Convective term of the momentum equation
  if(iconv(iu).ge.1) then

    icliup = iclrtp(iu ,icoef)
    iclivp = iclrtp(iv ,icoef)
    icliwp = iclrtp(iw ,icoef)

    init   = 1
    inc    = 1
    iccocg = 1
    iflmb0 = 1
    nswrgp = nswrgr(iu)
    imligp = imligr(iu)
    iwarnp = iwarni(iu)
    epsrgp = epsrgr(iu)
    climgp = climgr(iu)
    extrap = extrag(iu)

    imaspe = 1
    itypfl = 1

!     Mass flux calculation
    call inimas                                                   &
    !==========
 ( iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   crom   , brom   ,                                              &
   rtpa (1,iu)  , rtpa (1,iv)  , rtpa (1,iw)  ,                   &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   coefb(1,icliup) , coefb(1,iclivp) , coefb(1,icliwp) ,          &
   flumas , flumab )

!     Calculation of the convective term along the three spatial directions
!       without reconstruction
    do isou = 1, 3
      if(isou.eq.1) ivar0  = iu
      if(isou.eq.1) iclvar = icliup
      if(isou.eq.2) ivar0  = iv
      if(isou.eq.2) iclvar = iclivp
      if(isou.eq.3) ivar0  = iw
      if(isou.eq.3) iclvar = icliwp

      do ifac = 1, nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        flui = 0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
        fluj = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
        trflms(ifac) = -(flui*rtpa(ii,ivar0)+fluj*rtpa(jj,ivar0))
      enddo

      do ifac = 1, nfabor
        ii = ifabor(ifac)
        flui = 0.5d0*( flumab(ifac) +abs(flumab(ifac)) )
        fluj = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )
        pfac = coefa(ifac,iclvar)                                 &
             + coefb(ifac,iclvar)*rtpa(ii,ivar0)
        trflmb(ifac) = - ( flui*rtpa(ii,ivar0) + fluj*pfac )
      enddo

      init = 0
      if(isou.eq.1) then
        call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,          &
             ifacel,ifabor,trflms,trflmb,w10)
      elseif(isou.eq.2) then
        call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,          &
             ifacel,ifabor,trflms,trflmb,w11)
      elseif(isou.eq.3) then
        call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,          &
             ifacel,ifabor,trflms,trflmb,w12)
      endif

    enddo

  endif


! --- Viscous term

  if( idiff(iu).ge.1 ) then

    do iel = 1, ncelet
      w8(iel) = 1.d0
      w9(iel) = 0.d0
    enddo

    call cfdivs                                                   &
    !==========
 ( rtpa   , propce ,                                              &
   w10    , w8     , w9     , w9     )
   !        ------

    call cfdivs                                                   &
    !==========
 ( rtpa   , propce ,                                              &
   w11    , w9     , w8     , w9     )
   !        ------

    call cfdivs                                                   &
    !==========
 ( rtpa   , propce ,                                              &
   w12    , w9     , w9     , w8     )
   !        ------

  endif

!FIXME
! --- Mass sourceterm
!     All is explicit for the moment... to be changed when this
!      term is tested (see remark at the begining of this file)

  iiun   = 1
  iextts = 0
  thetv  = 0.d0
  do isou = 1, 3
    if(isou.eq.1) then
      ivar0  = iu
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr)   ,                             &
   w10   , w1     , w2 )
      do iel = 1, ncel
        w10(iel) = w10(iel) + w2(iel)
      enddo

    elseif(isou.eq.2) then
      ivar0  = iv
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr)   ,                             &
   w11   , w1     , w2 )
      do iel = 1, ncel
        w11(iel) = w11(iel) + w2(iel)
      enddo

    elseif(isou.eq.3) then
      ivar0  = iw
      call catsma                                                 &
      !==========
 ( ncelet, ncel   , ncesmp , iiun   , iextts , thetv  ,           &
   icetsm, itypsm(1,ivar0) , volume , rtpa(1,ivar0)   ,           &
   smacel(1,ivar0), smacel(1,ipr)   ,                             &
   w12   , w1     , w2 )
      do iel = 1, ncel
        w12(iel) = w12(iel) + w2(iel)
      enddo

    endif

  enddo

endif
!     End of the test on momentum source terms


! --- Volumic forces term (gravity)
do iel = 1, ncel
  rom = crom(iel)
  w5(iel) = gx + w10(iel)/rom
  w6(iel) = gy + w11(iel)/rom
  w7(iel) = gz + w12(iel)/rom
enddo

! 2.2 MASS FLUX CALCULATION AT THE FACES
! ======================================

! --- Calculation of the convective "velocities at the cell centers
!     (Calculation of u^n+dt*f^n)

do iel = 1, ncel
  w10(iel) = rtpa(iel,iu) + dt(iel)*w5(iel)
  w11(iel) = rtpa(iel,iv) + dt(iel)*w6(iel)
  w12(iel) = rtpa(iel,iw) + dt(iel)*w7(iel)
enddo

! --- Calculation of the flux with an INIMAS call

!     In order to avoid a misfit boundary condition, we impose a homogeneous
!       Neumann condition. Note that it is only usefull for gradient
!       reconstruction. The boundary value does not matter since the flux
!       is updated afterwards.

!     To avoid the declaration of other arrays, we build a zero flux with
!       nul avec :
!       TRFLMB = 1 pour COEFB = 1 et
!       INC    = 0 pour COEFA = 0

!     We take ROM = W1 = 1 and ROMB = TRFLMB = 1
!       (ROMB is also usefull for COEFB=1)

do iel = 1, ncel
  w1(iel) = 1.d0
enddo
do ifac = 1, nfabor
  trflmb(ifac) = 1.d0
enddo

icliup = iclrtp(iu ,icoef)
iclivp = iclrtp(iv ,icoef)
icliwp = iclrtp(iw ,icoef)

init   = 1
!              ^ Initialization of the mass flux
inc    = 0
!              ^ As mentioned above, for a zero flux
iccocg = 1
ivar0  = 0
iflmb0 = 1
nswrgp = 0
!              ^ Reconstruction is useless here
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

imaspe = 1
itypfl = 1

call inimas                                                       &
!==========
 ( ivar0  , ivar0  , ivar0  , imaspe , itypfl ,                   &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrgp , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   w1     , trflmb ,                                              &
   w10    , w11    , w12    ,                                     &
   coefa(1,icliup) , coefa(1,iclivp) , coefa(1,icliwp) ,          &
   trflmb          , trflmb          , trflmb          ,          &
   flumas , flumab )

! Free memory
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)
deallocate(w10, w11, w12)

!--------
! FORMATS
!--------


!----
! FIN
!----

return

end subroutine

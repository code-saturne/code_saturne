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

subroutine atprke &
!================

 ( nscal  ,                                                       &
   ipcvto,                                                        &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   tinstk , tinste )

!===============================================================================
! FONCTION :
! --------
!  SPECIFIQUE AU CAS ATMOSPHERIQUE :
!  CALCUL DU TERME DE PRODUCTION LIEE A LA FLOTTABILITE :
!  G = G*GRAD(THETA)/PRDTUR/THETA
!-------------------------------------------------------------------------------
! ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscal            ! i  ! <-- ! total number of scalars                        !
! ipcvto           ! i  ! <-- ! pointer for turbulent viscosity
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !              faces de bord                     !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !              faces de bord                     !
! tinstk(ncelet)   ! tr ! --> ! prod et terme de gravite pour eq k             !
! tinste(ncelet)   ! tr ! --> ! prod et terme de gravite pour eq eps           !
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
use dimens, only: ndimfb
use numvar
use entsor
use cstnum
use cstphy
use optcal
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use atincl

!===============================================================================

implicit none

! Arguments

integer          nscal
integer          ipcvto



double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision rtp (ncelet,*), rtpa (ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tinstk(ncelet), tinste(ncelet)

! Local variables
integer         iel
integer         itpp , iqw, icltpp,iclqw
integer         ipcliq
integer         iccocg, inc
integer         iivar
integer         nswrgp, imligp
integer         iwarnp

double precision gravke, prdtur
double precision theta_virt
double precision qldia,qw
double precision epsrgp, climgp, extrap

double precision xent,yent,zent,dum
double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1. Initialisation
!===============================================================================

! Allocate work arrays
allocate(grad(ncelet,3))


!===============================================================================
! 2. Calcul des derivees de la temperature potentielle
!===============================================================================

if (ippmod(iatmos).eq.1) then
  call dry_atmosphere()
elseif (ippmod(iatmos).eq.2) then
  call humid_atmosphere()
endif

! Deallocate work arrays
deallocate(grad)

return
contains

!**************************************************************************
!* repository of tactical functions/subroutines
!**************************************************************************
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------

subroutine dry_atmosphere()

! computes the production term in case of dry atmosphere
! ie. when ippmod(iatmos) eq 1

! Computation of the gradient of the potential temperature

itpp = isca(iscalt)
icltpp = iclrtp(itpp,icoef)

! computational options:

iccocg = 1
inc = 1

nswrgp = nswrgr(itpp)
epsrgp = epsrgr(itpp)
imligp = imligr(itpp)
iwarnp = iwarni(itpp)
climgp = climgr(itpp)
extrap = extrag(itpp)

iivar = itpp

! computes the turbulent production/destruction terms:
! dry atmo: (1/sigmas*theta)*(dtheta/dz)*gz

call grdcel                                                     &
     !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp ,imligp,            &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,itpp), coefa(1,icltpp) , coefb(1,icltpp) ,              &
   grad   )

! Production and gravity terms
! TINSTK=P+G et TINSTE = P + (1-CE3)*G

if (iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

! En production lineaire, on multiplie tout de suite le terme
! de gravite par VISCT, car le reste est deja multiplie.
! Dans les autres cas, la multiplication est faite plus tard.

if (iturb.eq.21) then
  do iel = 1, ncel
    gravke = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) &
           / (rtpa(iel,itpp)*prdtur)
    tinste(iel) = tinstk(iel) + propce(iel,ipcvto)*max(gravke,zero)
    tinstk(iel) = tinstk(iel) + propce(iel,ipcvto)*gravke
  enddo
else
  do iel = 1, ncel
    gravke = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) &
           / (rtpa(iel,itpp)*prdtur)
    tinste(iel) = tinstk(iel) + max(gravke,zero)
    tinstk(iel) = tinstk(iel) + gravke
  enddo
endif
end subroutine dry_atmosphere

!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------

subroutine humid_atmosphere()

! computes the production term in case of humid atmosphere
! ie. when ippmod(iatmos) eq 2

implicit none
double precision pphy

double precision, dimension(:), allocatable :: etheta
double precision, dimension(:), allocatable :: eq
double precision, dimension(:), allocatable :: gravke_theta
double precision, dimension(:), allocatable :: gravke_qw

allocate(etheta(ncelet))
allocate(eq(ncelet))
allocate(gravke_theta(ncelet))
allocate(gravke_qw(ncelet))

! Computation of the gradient of the potential temperature

itpp = isca(iscalt)
icltpp = iclrtp(itpp,icoef)
iqw = isca(itotwt)
iclqw = iclrtp(iqw,icoef)
ipcliq = ipproc(iliqwt)

! compute the coefficients etheta,eq

do iel = 1, ncel
  ! calculate the physical pressure 'pphy'
  if (imeteo.eq.0) then
    call atmstd(xyzcen(3,iel),pphy,dum,dum)
  else
    call intprf (                                                 &
         nbmett, nbmetm,                                          &
         ztmet, tmmet, phmet, xyzcen(3,iel), ttcabs, pphy )
  endif
  qw = rtpa(iel,iqw) ! total water content
  qldia = propce(iel,ipcliq) ! liquid water content
  call etheq(pphy,rtpa(iel,itpp),qw,qldia,                      &
             nebdia(iel),nn(iel),etheta(iel),eq(iel))
enddo
! options for gradient calculation

iccocg = 1
inc = 1

nswrgp = nswrgr(itpp)
epsrgp = epsrgr(itpp)
imligp = imligr(itpp)
iwarnp = iwarni(itpp)
climgp = climgr(itpp)
extrap = extrag(itpp)

iivar = itpp

! computes the turbulent production/destruction terms:
! humid atmo: (1/sigmas*theta_v)*(dtheta_l/dz)*gz

call grdcel                                                       &
     !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp ,imligp,            &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,itpp), coefa(1,icltpp) , coefb(1,icltpp) ,              &
   grad   )

! Production and gravity terms
! TINSTK = P + G et TINSTE = P + (1-CE3)*G

if(iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

! store now the production term due to theta_liq in gravke_theta

if (iturb.eq.21) then
  ! For linear production multiply immediately by the turbulent
  ! viscosity VISCT.
  do iel = 1, ncel
    qw = rtpa(iel,iqw) ! total water content
    qldia = propce(iel,ipcliq) ! liquid water content
    theta_virt = rtpa(iel,itpp)*(1.d0 + (rvsra-1)*qw - rvsra*qldia)
    gravke = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz)            &
           / (theta_virt*prdtur)
    gravke_theta(iel) = propce(iel,ipcvto)*gravke*etheta(iel)
  enddo
else
  ! For the other cases no multiplication done
  do iel = 1, ncel
    qw = rtpa(iel,iqw) ! total water content
    qldia = propce(iel,ipcliq) ! liquid water content
    theta_virt = rtpa(iel,itpp)*(1.d0 + (rvsra - 1)*qw - rvsra*qldia)
    gravke = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz)            &
           / (theta_virt*prdtur)
    gravke_theta(iel) = gravke*etheta(iel)
  enddo
endif

! ----------------------------------------------------------------
! now gradient of humidity and it's associated production term
! ----------------------------------------------------------------

nswrgp = nswrgr(iqw)
epsrgp = epsrgr(iqw)
imligp = imligr(iqw)
iwarnp = iwarni(iqw)
climgp = climgr(iqw)
extrap = extrag(iqw)

iivar = iqw

! computes the turbulent production/destruction terms:
! humid atmo: (1/sigmas*theta_v)*(dtheta_l/dz)*gz

call grdcel                                                       &
  !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp ,imligp,            &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,iqw), coefa(1,iclqw) , coefb(1,iclqw) ,                 &
   grad   )

! Production and gravity terms
! TINSTK = P + G et TINSTE = P + (1-CE3)*G

if (iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

! store the production term due to qw in gravke_qw

if (iturb.eq.21) then
  ! For linear production multiply immediately by the turbulent
  ! viscosity VISCT.
  do iel = 1, ncel
    qw = rtpa(iel,iqw) ! total water content
    qldia = propce(iel,ipcliq) !liquid water content
    theta_virt = rtpa(iel,itpp)*(1.d0 + (rvsra - 1.d0)*qw - rvsra*qldia)
    gravke = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz)                 &
           / (theta_virt*prdtur)
    gravke_qw(iel) = propce(iel,ipcvto)*gravke*eq(iel)
  enddo
else
  ! For the other cases no multiplication done
  do iel = 1, ncel
    qw = rtpa(iel,iqw) ! total water content
    qldia = propce(iel,ipcliq) ! liquid water content
    theta_virt = rtpa(iel,itpp)*(1.d0 + (rvsra - 1.d0)*qw - rvsra*qldia)
    gravke = (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz)                 &
           / (theta_virt*prdtur)
    gravke_qw(iel) = gravke*eq(iel)
  enddo
endif

! termination

do iel = 1, ncel
  gravke = gravke_theta(iel) + gravke_qw(iel)
  tinste(iel) = tinstk(iel) + max(gravke,zero)
  tinstk(iel) = tinstk(iel) + gravke
enddo
end subroutine humid_atmosphere

end subroutine atprke

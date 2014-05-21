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

subroutine cfmspr &
!================

 ( nvar   , nscal  , iterns , ncepdp , ncesmp ,                   &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , vela   ,                   &
   ckupdc , smacel )

!===============================================================================
! FONCTION :
! ----------

! SOLVING OF A CONVECTION-DIFFUSION EQUATION WITH SOURCE TERMS
!   FOR PRESSURE ON ONE TIME-STEP
!   (COMPRESSIBLE ALGORITHM IN P,U,E)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iterns           ! i  ! <-- ! Navier-Stokes iteration number                 !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet)   ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp)   ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! vela             ! ra ! <-- ! variable value at time step beginning          !
! ckupdc           ! tr ! <-- ! work array for the head loss                   !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! variable value associated to the mass source   !
! (ncesmp,*   )    !    !     ! term (for ivar=ipr, smacel is the mass flux    !
!                  !    !     ! \f$ \Gamma^n \f$)                              !
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
use pointe, only: itypfb
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal, iterns
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision vela  (3  ,ncelet)

! Local variables

character*80     chaine
integer          ivar
integer          ifac  , iel
integer          init  , inc   , iccocg, isqrt , ii, jj
integer          iflmas, iflmab
integer          ippvar, ipp   , iphydp, icvflb
integer          nswrgp, imligp, iwarnp
integer          istatp, iconvp, idiffp, ireslp, ndircp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          ivoid(1)
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision epsrsp
double precision sclnor

integer          imucpp, idftnp, iswdyp
integer          imvis1

double precision thetv, relaxp, hint

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: smbrs, rovsdt
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, allocatable, dimension(:) :: w10
double precision, allocatable, dimension(:) :: wflmas, wflmab
double precision, allocatable, dimension(:) :: wbfa, wbfb
double precision, allocatable, dimension(:) :: dpvar

double precision, allocatable, dimension(:) :: c2

double precision, dimension(:), pointer :: coefaf_p, coefbf_p
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: rhopre, crom
double precision, dimension(:), pointer :: cvar_pr

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate temporary arrays for the mass resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbrs(ncelet), rovsdt(ncelet))
allocate(wflmas(nfac), wflmab(nfabor))

allocate(wbfa(nfabor), wbfb(nfabor))

! Allocate work arrays
allocate(w1(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(w10(ncelet))
allocate(dpvar(ncelet))

! --- Number of computational variable and post for pressure

ivar   = ipr
ippvar = ipprtp(ivar)

! --- Mass flux associated to energy
call field_get_key_int(ivarfl(isca(ienerg)), kimasf, iflmas)
call field_get_key_int(ivarfl(isca(ienerg)), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_s(icrom, crom)
call field_get_val_prev_s(icrom, rhopre)

call field_get_val_s(ivarfl(ipr), cvar_pr)

call field_get_label(ivarfl(ivar), chaine)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

call field_get_coefaf_s(ivarfl(ipr), coefaf_p)
call field_get_coefbf_s(ivarfl(ipr), coefbf_p)

! Computation of the boundary coefficients for the pressure gradient
! recontruction in accordance with the diffusion boundary coefficients (coefaf_p,
! coefbf_p). Always a homogeneous Neumann except at walls where the hydrostatic
! pressure is taken into account (icfgrp option).
do ifac = 1, nfabor
  iel = ifabor(ifac)
  hint = dt(iel) / distb(ifac)
  wbfa(ifac) = -coefaf_p(ifac) / hint
  wbfb(ifac) = 1.d0
enddo

!===============================================================================
! 2. SOURCE TERMS
!===============================================================================

! --> Initialization

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo


!     MASS SOURCE TERM
!     ================

if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    smbrs(iel) = smbrs(iel) + smacel(iel,ipr)*volume(iel)
  enddo
endif


!     UNSTEADY TERM
!     =============

! --- Calculation of the square of sound velocity c2.
!     Pressure is an unsteady variable in this algorithm
!     Varpos has been modified for that.

allocate(c2(ncelet))

call cf_thermo_c_square( cvar_pr, crom, c2, ncel)
!======================

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel) + istat(ivar)*(volume(iel)/(dt(iel)*c2(iel)))
enddo

!===============================================================================
! 3. "MASS FLUX" AND FACE "VISCOSITY" CALCULATION
!===============================================================================

! Computation of the "convective flux" for the density:
! (u + dt f) is computed at internal faces and stored in wflmas,
! the values at boundary faces in wflmab won't be taken into account.

call cfmsfp                                                                     &
!==========
( nvar   , nscal  , iterns , ncepdp , ncesmp ,                                  &
  icepdc , icetsm , itypsm ,                                                    &
  dt     , propce , vela   ,                                                    &
  ckupdc , smacel ,                                                             &
  wflmas , wflmab )

! Mass flux at internal faces (upwind scheme for the density).
do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  wflmas(ifac) = -0.5d0*                                                        &
                 ( crom(ii)*(wflmas(ifac)+abs(wflmas(ifac)))                    &
                 + crom(jj)*(wflmas(ifac)-abs(wflmas(ifac))))
enddo

! Mass flux at boundary faces.
if (icfgrp.eq.1) then
  ! The hydrostatic pressure gradient contribution has to be added to the mass
  ! flux if it has been taken into account in the pressure B.C. at walls.
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    if (itypfb(ifac).eq.iparoi) then
      wflmab(ifac) = -bmasfl(ifac)                                              &
                     -dt(iel)*crom(iel)*(  gx*surfbo(1,ifac)                    &
                                         + gy*surfbo(2,ifac)                    &
                                         + gz*surfbo(3,ifac))
    else
      wflmab(ifac) = -bmasfl(ifac)
    endif
  enddo
else
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    wflmab(ifac) = -bmasfl(ifac)
  enddo
endif

init = 0
call divmas(init,wflmas,wflmab,smbrs)

! (Delta t)_ij is calculated as the "viscocity" associated to the pressure
imvis1 = 1

call viscfa                                                                     &
!==========
( imvis1 ,                                                                      &
  dt     ,                                                       &
  viscf  , viscb  )

!===============================================================================
! 4. SOLVING
!===============================================================================

istatp = istat (ivar)
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
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetv  = thetav(ivar)
icvflb = 0

call codits                                                                     &
!==========
( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap ,                &
  imrgra , nswrsp , nswrgp , imligp , ircflp ,                                  &
  ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,                         &
  imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                                  &
  blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,                         &
  relaxp , thetv  ,                                                             &
  rtpa(1,ivar)    , rtpa(1,ivar)    ,                                           &
  wbfa   , wbfb   ,                                                             &
  coefaf_p        , coefbf_p        ,                                           &
  wflmas          , wflmab          ,                                           &
  viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,                         &
  rvoid  , rvoid  ,                                                             &
  icvflb , ivoid  ,                                                             &
  rovsdt , smbrs  , rtp(1,ivar)     , dpvar  ,                                  &
  rvoid  , rvoid  )

!===============================================================================
! 5. PRINTINGS AND CLIPPINGS
!===============================================================================

! --- User intervention for a finer management of the bounds and possible
!       corrective treatement.

call cf_check_pressure(cvar_pr, ncel)
!=====================

! --- Explicit balance (see codits : the increment is withdrawn)

if (iwarni(ivar).ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                     &
                 - istat(ivar)*(volume(iel)/dt(iel))                            &
                   *(rtp(iel,ivar)-rtpa(iel,ivar))                              &
                   * max(0,min(nswrsm(ivar)-2,1))
  enddo
  isqrt = 1
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200) chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. COMMUNICATION OF P
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ivar))
endif

!===============================================================================
! 7. ACOUSTIC MASS FLUX CALCULATION AT THE FACES
!===============================================================================

! mass flux = [dt (grad P).n] + [rho (u + dt f)]

! Computation of [dt (grad P).n] by itrmas
init   = 1
inc    = 1
iccocg = 1
iphydp = 0
! festx,y,z    = rvoid
! viscf, viscb = arithmetic mean at faces
! viscelx,y,z  = dt
! This flux is stored as the mass flux of the energy

call itrmas                                                                     &
!==========
( init   , inc    , imrgra , iccocg , nswrgp , imligp ,                         &
  iphydp , iwarnp ,                                                             &
  epsrgp , climgp , extrap ,                                                    &
  rvoid  ,                                                                      &
  rtp(1,ivar)     ,                                                             &
  wbfa   , wbfb   ,                                                             &
  coefaf_p        , coefbf_p        ,                                           &
  viscf  , viscb  ,                                                             &
  dt     , dt     , dt     ,                                                    &
  imasfl, bmasfl)

! Incrementation of the flux with [rho (u + dt f)].n = wflmas
! (added with a negative sign since wflmas,wflmab was used above
! in the right hand side).
do ifac = 1, nfac
  imasfl(ifac) = imasfl(ifac) - wflmas(ifac)
enddo
do ifac = 1, nfabor
  bmasfl(ifac) = bmasfl(ifac) - wflmab(ifac)
enddo

!===============================================================================
! 8. UPDATING OF THE DENSITY
!===============================================================================

if (igrdpp.gt.0) then

  do iel = 1, ncel
    ! Backup of the current density values
    rhopre(iel) = crom(iel)
    ! Update of density values
    crom(iel) = crom(iel)+(rtp(iel,ivar)-rtpa(iel,ivar))/c2(iel)
  enddo

!===============================================================================
! 9. DENSITY COMMUNICATION
!===============================================================================

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(crom)
    call synsca(rhopre)
  endif

endif

! There are no clippings on density because we consider that pressure
! and energy have been checked and are correct so that the density
! is also correct through the state law or the linearized law.

deallocate(c2)
deallocate(viscf, viscb)
deallocate(smbrs, rovsdt)
deallocate(wflmas, wflmab)
deallocate(w1)
deallocate(w7, w8, w9)
deallocate(w10)
deallocate(dpvar)
deallocate(wbfa, wbfb)

!--------
! FORMATS
!--------

 1000 format(/,                                                                 &
'   ** RESOLUTION FOR THE VARIABLE ',A8                        ,/,              &
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)

!----
! END
!----

return
end subroutine

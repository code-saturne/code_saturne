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

subroutine cfdivs &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtp    , propce , propfb ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   diverg , ux     , uy     , uz     )

!===============================================================================
! FONCTION :
! ----------
!                                   v
! CALCULE DIVERG = DIVERG + DIV(SIGMA .U)

!          v               t
! AVEC SIGMA = MU (GRAD(U) + GRAD(U)) + (KAPPA - 2/3 MU) DIV(U) Id

! ET MU = MU_LAMINAIRE + MU_TURBULENT

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! diverg(ncelet    ! tr ! --> ! div(sigma.u)                                   !
! ux,y,z(ncelet    ! tr ! <-- ! composantes du vecteur u                       !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use entsor
use numvar
use optcal
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use pointe, only:coefau, coefbu

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision rtp(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision diverg(ncelet)
double precision ux(ncelet), uy(ncelet), uz(ncelet)

! Local variables

integer          inc, iel, ifac, ii, jj
integer          iclvar
integer          nswrgp, imligp, iwarnp
integer          ipcvis, ipcvst, ipcvsv

logical          ilved

double precision epsrgp, climgp, extrap
double precision vecfac, kappa, mu, trgdru
double precision sigma(3,3)

double precision, allocatable, dimension(:) :: vistot
double precision, allocatable, dimension(:,:,:) :: gradv
double precision, allocatable, dimension(:,:) :: tempv

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays
allocate(vistot(ncelet))
allocate(gradv(ncelet,3,3))
allocate(tempv(3, ncelet))

ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
if(iviscv.gt.0) then
  ipcvsv = ipproc(iviscv)
else
  ipcvsv = 0
endif

! --- Calcul de la viscosite totale

if (itytur.eq.3 ) then
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis)
  enddo
else
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)
  enddo
endif

! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(vistot)
  if (ipcvsv.gt.0) then
    call synsca(propce(1,ipcvsv))
  endif
endif

!===============================================================================
! 2. Compute the divegence of (sigma.u)
!===============================================================================

! --- Velocity gradient
inc    = 1
nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(iu)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

ilved = .false.

! WARNING: gradv(iel, xyz, uvw)
call grdvec &
!==========
( iu     , imrgra , inc    , nswrgp , imligp ,                   &
  iwarnp , nfecra ,                                              &
  epsrgp , climgp , extrap ,                                     &
  ilved  ,                                                       &
  rtp(1,iu) ,  coefau , coefbu,                                  &
  gradv  )

! --- Compute the vector \tens{\sigma}.\vect{v}
!     i.e. sigma_ij v_j e_i

! Variable kappa in space
if (ipcvsv.gt.0) then
  do iel = 1, ncel
    kappa = propce(iel,ipcvsv)
    mu = vistot(iel)
    trgdru = gradv(iel, 1, 1)+gradv(iel, 2, 2)+gradv(iel, 3, 3)

    sigma(1, 1) = mu * 2.d0*gradv(iel, 1, 1)  &
                - (kappa-2.d0/3.d0)*trgdru

    sigma(2, 2) = mu * 2.d0*gradv(iel, 2, 2)  &
                - (kappa-2.d0/3.d0)*trgdru

    sigma(3, 3) = mu * 2.d0*gradv(iel, 3, 3)  &
                - (kappa-2.d0/3.d0)*trgdru

    sigma(1, 2) = mu * (gradv(iel, 1, 2) + gradv(iel, 2, 1))
    sigma(2, 1) = sigma(1, 2)

    sigma(2, 3) = mu * (gradv(iel, 2, 3) + gradv(iel, 3, 2))
    sigma(3, 2) = sigma(2, 3)

    sigma(1, 3) = mu * (gradv(iel, 1, 3) + gradv(iel, 3, 1))
    sigma(3, 1) = sigma(1, 3)

    tempv(1, iel) = sigma(1, 1)*ux(iel) &
                  + sigma(1, 2)*uy(iel) &
                  + sigma(1, 3)*uz(iel)
    tempv(2, iel) = sigma(2, 1)*ux(iel) &
                  + sigma(2, 2)*uy(iel) &
                  + sigma(2, 3)*uz(iel)
    tempv(3, iel) = sigma(3, 1)*ux(iel) &
                  + sigma(3, 2)*uy(iel) &
                  + sigma(3, 3)*uz(iel)
  enddo

else

  do iel = 1, ncel
    kappa = viscv0
    mu = vistot(iel)
    trgdru = gradv(iel, 1, 1)+gradv(iel, 2, 2)+gradv(iel, 3, 3)

    sigma(1, 1) = mu * 2.d0*gradv(iel, 1, 1)  &
                - (kappa-2.d0/3.d0)*trgdru

    sigma(2, 2) = mu * 2.d0*gradv(iel, 2, 2)  &
                - (kappa-2.d0/3.d0)*trgdru

    sigma(3, 3) = mu * 2.d0*gradv(iel, 3, 3)  &
                - (kappa-2.d0/3.d0)*trgdru

    sigma(1, 2) = mu * (gradv(iel, 1, 2) + gradv(iel, 2, 1))
    sigma(2, 1) = sigma(1, 2)

    sigma(2, 3) = mu * (gradv(iel, 2, 3) + gradv(iel, 3, 2))
    sigma(3, 2) = sigma(2, 3)

    sigma(1, 3) = mu * (gradv(iel, 1, 3) + gradv(iel, 3, 1))
    sigma(3, 1) = sigma(1, 3)

    tempv(1, iel) = sigma(1, 1)*ux(iel) &
                  + sigma(1, 2)*uy(iel) &
                  + sigma(1, 3)*uz(iel)
    tempv(2, iel) = sigma(2, 1)*ux(iel) &
                  + sigma(2, 2)*uy(iel) &
                  + sigma(2, 3)*uz(iel)
    tempv(3, iel) = sigma(3, 1)*ux(iel) &
                  + sigma(3, 2)*uy(iel) &
                  + sigma(3, 3)*uz(iel)
  enddo

endif

! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(tempv)
endif

! Initialize diverg(ncel+1, ncelet)
!  (unused value, but need to be initialized to avoid Nan values)
if (ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    diverg(iel) = 0.d0
  enddo
endif

! --- Interior faces contribution

do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  vecfac = surfac(1,ifac)*(tempv(1, ii)+tempv(1, jj))*0.5d0               &
         + surfac(2,ifac)*(tempv(2, ii)+tempv(2, jj))*0.5d0               &
         + surfac(3,ifac)*(tempv(3, ii)+tempv(3, jj))*0.5d0
  diverg(ii) = diverg(ii) + vecfac
  diverg(jj) = diverg(jj) - vecfac
enddo

! --- Boundary faces contribution

do ifac = 1, nfabor
  ii = ifabor(ifac)
  vecfac = surfbo(1,ifac)*tempv(1, ii)                                    &
         + surfbo(2,ifac)*tempv(2, ii)                                    &
         + surfbo(3,ifac)*tempv(3, ii)
  diverg(ii) = diverg(ii) + vecfac
enddo

return

end subroutine
